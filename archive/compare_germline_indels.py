#!/usr/bin/env python3
"""
Step 2: Compare Indels Called Germline By Our Pipeline to Those Called Germline By DeepVariant

Simple approach: Process chromosome by chromosome, position by position
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from datetime import datetime
import time

# ====================
# CONFIGURATION
# ====================

# Team B indel caller data
TEAM_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/results"
TEAM_SAMPLE_FILE = "sample1_indels_only.parquet"

# DeepVariant data (already converted to single-indel format)
DV_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/deep_variant"
DV_SAMPLE_FILE = "GTEx-sample1.indels.parquet"

# Output directory
OUTPUT_DIR = "/home/ubuntu/data/indel_comparison_results"

# ====================
# UTILITY FUNCTIONS
# ====================

def convert_to_single_indel_format(ref, alt):
    """
    Convert from combined VCF format to single-indel format for comparison.
    """
    anchor = ref[0]  # First base is always the anchor

    if len(alt) > len(ref):
        # INSERTION: bases inserted after anchor
        inserted_length = len(alt) - len(ref)
        inserted_bases = alt[1:1+inserted_length]
        new_ref = anchor
        new_alt = anchor + inserted_bases
        return new_ref, new_alt

    elif len(alt) < len(ref):
        # DELETION: bases deleted after anchor
        deleted_length = len(ref) - len(alt)
        deleted_bases = ref[1:1+deleted_length]
        new_ref = anchor + deleted_bases
        new_alt = anchor
        return new_ref, new_alt

    else:
        # Same length - shouldn't happen for indels
        return ref, alt


def main():
    """
    Main function to run the comparison.
    """
    print("=" * 80)
    print("STEP 2: COMPARE GERMLINE INDEL CALLS")
    print("Team Pipeline vs DeepVariant")
    print("=" * 80)

    # Load Team B data
    print("\nLoading Team B data...")
    team_path = Path(TEAM_DATA_DIR) / TEAM_SAMPLE_FILE
    if not team_path.exists():
        print(f"ERROR: Team file not found: {team_path}")
        sys.exit(1)
    team_df = pd.read_parquet(team_path)
    print(f"  Loaded {len(team_df):,} rows from Team B parquet")

    # Load DeepVariant data
    print("\nLoading DeepVariant data...")
    dv_path = Path(DV_DATA_DIR) / DV_SAMPLE_FILE
    if not dv_path.exists():
        print(f"ERROR: DV file not found: {dv_path}")
        sys.exit(1)
    dv_df = pd.read_parquet(dv_path)
    print(f"  Loaded {len(dv_df):,} rows from DeepVariant parquet")

    # Filter Team B for germline ALT indels in not-difficult regions
    print("\nFiltering Team B data...")
    team_germline_alts = team_df[
        (team_df['status'] == 'GERMLINE') &
        (team_df['allele_type'] == 'ALT') &
        (team_df['in_notdifficult'] == True)
    ].copy()
    print(f"  Team B germline ALT indels in not-difficult regions: {len(team_germline_alts):,}")

    # Get unique chromosomes
    chromosomes = sorted(set(team_germline_alts['chrom'].unique()) | set(dv_df['chrom'].unique()))
    print(f"\nChromosomes to process: {', '.join(chromosomes)}")

    # Initialize tracking lists for all matches
    all_matches = []
    total_team_alts = 0
    total_dv_alts = 0

    print("\n" + "=" * 80)
    print("PROCESSING BY CHROMOSOME")
    print("=" * 80)

    for chrom in chromosomes:
        start_chr = time.time()
        print(f"\nProcessing {chrom}...")

        # Get ALL team data for this chromosome (for position lookups)
        team_chrom_all = team_df[team_df['chrom'] == chrom]

        # Filter Team B germline ALTs for counting
        team_germline_chrom = team_germline_alts[team_germline_alts['chrom'] == chrom]

        # Filter DV data for this chromosome
        dv_chrom = dv_df[dv_df['chrom'] == chrom]

        # Apply DV filtering on the smaller chromosome subset (much faster!)
        dv_germline = dv_chrom[
            (dv_chrom['filter'] == 'PASS') &
            #(dv_chrom['gt'] != '0/0') &  # Exclude homozygous reference
            #(dv_chrom['gt'] != './.') &  # Exclude no-calls
            (dv_chrom['in_notdifficult'] == True)
        ]

        print(f"  Team B: {len(team_germline_chrom):,} germline ALT indels")
        print(f"  DeepVariant: {len(dv_germline):,} PASS indels with variant genotypes")

        total_team_alts += len(team_germline_chrom)
        total_dv_alts += len(dv_germline)

        if len(team_germline_chrom) == 0 or len(dv_germline) == 0:
            print(f"  No overlap possible, skipping...")
            continue

        # Get unique positions with Team B germline indels
        team_positions = set(team_germline_chrom['pos'].unique())
        dv_positions = set(dv_germline['pos'].unique())

        # Find overlapping positions
        overlapping_positions = team_positions & dv_positions
        print(f"  Overlapping positions: {len(overlapping_positions):,}")

        # Process each overlapping position
        chr_matches = 0
        for pos in overlapping_positions:
            # Single lookup: Get ALL Team B data at this position
            team_pos_all = team_chrom_all[team_chrom_all['pos'] == pos]

            # Extract REF from these rows
            ref_rows = team_pos_all[team_pos_all['allele_type'] == 'REF']
            if len(ref_rows) == 0:
                continue
            ref_seq = ref_rows.iloc[0]['allele_seq']

            # Extract germline ALT rows from these rows
            team_pos_alts = team_pos_all[
                (team_pos_all['allele_type'] == 'ALT') &
                (team_pos_all['status'] == 'GERMLINE') &
                (team_pos_all['in_notdifficult'] == True)
            ]

            if len(team_pos_alts) == 0:
                continue

            # Check if multi-allelic site
            has_multiple_alts = len(team_pos_alts) > 1

            # Get DeepVariant data at this position
            dv_pos_data = dv_germline[dv_germline['pos'] == pos]

            # Compare each Team B ALT with each DV entry
            for _, team_alt in team_pos_alts.iterrows():
                alt_seq = team_alt['allele_seq']

                # Apply conversion if multi-allelic site
                if has_multiple_alts:
                    converted_ref, converted_alt = convert_to_single_indel_format(ref_seq, alt_seq)
                else:
                    converted_ref, converted_alt = ref_seq, alt_seq

                # Check for matches in DeepVariant
                for _, dv_row in dv_pos_data.iterrows():
                    if converted_ref == dv_row['ref'] and converted_alt == dv_row['alt']:
                        # Found a match!
                        chr_matches += 1
                        all_matches.append(f"{chrom}:{pos}:{converted_ref}:{converted_alt}")
                        break  # Found match for this Team B alt, move to next

        print(f"  Matches found: {chr_matches:,}")
        elapsed = time.time() - start_chr
        print(f"  Completed in {elapsed:.1f} seconds")

    # Calculate summary statistics
    print("\n" + "=" * 80)
    print("OVERALL STATISTICS")
    print("=" * 80)

    total_matches = len(all_matches)

    print(f"\n### BASIC COUNTS ###")
    print(f"Team B germline ALT indels: {total_team_alts:,}")
    print(f"DeepVariant germline indels: {total_dv_alts:,}")
    print(f"Exact matches: {total_matches:,}")

    if total_team_alts > 0:
        print(f"\n### MATCH RATES ###")
        print(f"Team B indels found in DV: {total_matches/total_team_alts*100:.1f}% ({total_matches:,}/{total_team_alts:,})")

    if total_dv_alts > 0:
        print(f"DV indels found in Team B: {total_matches/total_dv_alts*100:.1f}% ({total_matches:,}/{total_dv_alts:,})")

    # Save summary results
    if total_matches > 0:
        Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save list of matched variants
        matches_file = f"{OUTPUT_DIR}/matched_germline_indels_{timestamp}.txt"
        with open(matches_file, 'w') as f:
            f.write("# Matched germline indels (chrom:pos:ref:alt)\n")
            for match in sorted(all_matches):
                f.write(f"{match}\n")
        print(f"\n### OUTPUT FILES ###")
        print(f"Matched indels saved to: {matches_file}")

        # Save summary statistics
        summary_file = f"{OUTPUT_DIR}/summary_stats_{timestamp}.txt"
        with open(summary_file, 'w') as f:
            f.write("GERMLINE INDEL COMPARISON SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Team B germline ALT indels: {total_team_alts:,}\n")
            f.write(f"DeepVariant germline indels: {total_dv_alts:,}\n")
            f.write(f"Exact matches: {total_matches:,}\n\n")
            f.write(f"Team B match rate: {total_matches/total_team_alts*100:.1f}%\n")
            f.write(f"DV match rate: {total_matches/total_dv_alts*100:.1f}%\n")

        print(f"Summary saved to: {summary_file}")

    print("\n" + "=" * 80)
    print("Analysis complete!")


if __name__ == "__main__":
    main()