#!/usr/bin/env python3
"""
Compare Team B indel calls with DeepVariant calls
Focus on germline/somatic classification agreement
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from datetime import datetime

class Tee:
    """Class to duplicate stdout to both console and log file"""
    def __init__(self, filename):
        self.file = open(filename, 'w')
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
        self.file.flush()

    def flush(self):
        self.file.flush()

# Configuration
TEAM_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/results"
TEAM_SAMPLE_FILE = "sample1_indels_only.parquet"
DV_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/deep_variant"
DV_SAMPLE_FILE = "GTEx-sample1.indels.parquet"
OUTPUT_DIR = "/home/ubuntu/data/indel_comparison_results"

# Testing mode - set to None to process all chromosomes, or list specific ones
TESTING = None #['chr1']  # Examples: None, ['chr1'], ['chr1', 'chr2', 'chr21']

def convert_to_single_indel_format(ref, alt):
    """Convert from combined VCF format to single-indel format"""
    anchor = ref[0]
    if len(alt) > len(ref):
        # INSERTION
        inserted_length = len(alt) - len(ref)
        inserted_bases = alt[1:1+inserted_length]
        return anchor, anchor + inserted_bases
    elif len(alt) < len(ref):
        # DELETION
        deleted_length = len(ref) - len(alt)
        deleted_bases = ref[1:1+deleted_length]
        return anchor + deleted_bases, anchor
    else:
        return ref, alt

def main():
    # Set up logging
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{OUTPUT_DIR}/compare_DV_log_{timestamp}.txt"
    tee = Tee(log_file)

    print(f"Logging to: {log_file}")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    print("\nLoading data...")

    # Load Team B data
    team_path = Path(TEAM_DATA_DIR) / TEAM_SAMPLE_FILE
    team_df = pd.read_parquet(team_path)
    print(f"  Loaded {len(team_df):,} Team B rows")

    # Load DeepVariant data
    dv_path = Path(DV_DATA_DIR) / DV_SAMPLE_FILE
    dv_df = pd.read_parquet(dv_path)
    print(f"  Loaded {len(dv_df):,} DeepVariant rows")

    # Filter datasets
    print("\nFiltering datasets...")
    team_df = team_df[team_df['in_notdifficult'] == True].copy()
    print(f"  Team B after filtering: {len(team_df):,} rows")

    # Note that there are some with PASS but genotype 0/0; we accept those.
    dv_df = dv_df[(dv_df['in_notdifficult'] == True) & (dv_df['filter'] == 'PASS')].copy()
    print(f"  DeepVariant after filtering: {len(dv_df):,} rows")

    # Initialize counters
    false_positives_count = 0
    true_negative_count = 0
    false_negatives_pos_in_DV = 0
    false_negatives_pos_not_in_DV = 0
    true_positives_pos_in_DV = 0
    true_positives_pos_not_in_DV = 0

    # Initialize registers
    false_positive_register = []
    false_negative_register = []

    # Get chromosomes
    chromosomes = sorted(set(team_df['chrom'].unique()) | set(dv_df['chrom'].unique()))

    # Apply testing filter if enabled
    if TESTING is not None:
        chromosomes = [chr for chr in TESTING if chr in chromosomes]
        if not chromosomes:
            print(f"ERROR: None of the test chromosomes {TESTING} found in data")
            sys.exit(1)
        print(f"\nTESTING MODE: Processing {len(chromosomes)} chromosome(s): {', '.join(chromosomes)}")
    
    print(f"\nProcessing {len(chromosomes)} chromosomes...")

    for chrom in chromosomes:
        print(f"\n{chrom}:")
        # Filter for chromosome
        team_chrom = team_df[team_df['chrom'] == chrom]
        dv_chrom = dv_df[dv_df['chrom'] == chrom]

        # Get Team B positions with ALT rows (somatic or germline)
        team_alt_positions = set(team_chrom[
            (team_chrom['allele_type'] == 'ALT') &
            ((team_chrom['status'] == 'SOMATIC') | (team_chrom['status'] == 'GERMLINE'))
        ]['pos'].unique())

        # Get DV positions
        dv_positions = set(dv_chrom['pos'].unique())

        # Split Team B positions into those in DV and not in DV
        team_pos_in_dv = team_alt_positions & dv_positions
        team_pos_not_in_dv = team_alt_positions - dv_positions

        print(f"  Team ALT positions: {len(team_alt_positions):,}")
        print(f"  DV positions: {len(dv_positions):,}")
        print(f"  Overlapping: {len(team_pos_in_dv):,}")
        print(f"  Team-only: {len(team_pos_not_in_dv):,}")

        # Count positions not in DV (vectorized approach)
        # Filter for ALT first (cheaper), then check positions
        team_alts = team_chrom[team_chrom['allele_type'] == 'ALT']
        team_alts_not_in_dv = team_alts[team_alts['pos'].isin(team_pos_not_in_dv)]
        status_counts = team_alts_not_in_dv['status'].value_counts()

        false_negatives_pos_not_in_DV += status_counts.get('GERMLINE', 0)
        true_positives_pos_not_in_DV += status_counts.get('SOMATIC', 0)

        # Process overlapping positions
        for pos in team_pos_in_dv:
            # Get all Team B data at this position
            team_pos_all = team_chrom[team_chrom['pos'] == pos]

            # Get REF
            ref_rows = team_pos_all[team_pos_all['allele_type'] == 'REF']
            if len(ref_rows) == 1:
                ref_seq = ref_rows.iloc[0]['allele_seq']
            else:
                raise ValueError(f"Expected 1 REF but found {len(ref_rows)} for {chrom}:{pos}")
            
            # Get Team B ALTs
            team_alts = team_pos_all[(team_pos_all['allele_type'] == 'ALT')]

            # Check if multi-allelic
            has_multiple_alts = len(team_alts) > 1

            # Get DV data at this position (already filtered for PASS, we accept genotype 0/0)
            dv_germline_at_pos = dv_chrom[dv_chrom['pos'] == pos]

            # Process each Team B ALT
            for _, team_alt in team_alts.iterrows():
                alt_seq = team_alt['allele_seq']
                team_status = team_alt['status']

                # Apply conversion if multi-indel allelic resulting in complex REF/ALT
                if has_multiple_alts:
                    converted_ref, converted_alt = convert_to_single_indel_format(ref_seq, alt_seq)
                else:
                    converted_ref, converted_alt = ref_seq, alt_seq

                # Check for match in DV germline calls
                found_in_dv_germline = False
                for _, dv_row in dv_germline_at_pos.iterrows():
                    if converted_ref == dv_row['ref'] and converted_alt == dv_row['alt']:
                        found_in_dv_germline = True
                        break

                # Categorize
                if team_status == 'GERMLINE':
                    if found_in_dv_germline:
                        true_negative_count += 1  # Both agree germline
                    else:
                        false_negatives_pos_in_DV += 1  # Team says germline, DV doesn't
                        false_negative_register.append({
                            'chrom': chrom,
                            'pos': pos,
                            'ref': converted_ref,
                            'alt': converted_alt
                        })
                elif team_status == 'SOMATIC':
                    if found_in_dv_germline:
                        false_positives_count += 1  # Team says somatic, DV says germline
                        false_positive_register.append({
                            'chrom': chrom,
                            'pos': pos,
                            'ref': converted_ref,
                            'alt': converted_alt
                        })
                    else:
                        true_positives_pos_in_DV += 1  # Both agree not germline

    # Print results
    print("\nProcessing complete!")
    print("\n=== RESULTS ===")
    print("Note: Counts are of ALT alleles and do not take into account allele depth")
    print("-" * 60)
    print(f"False Positives (Team somatic, DV germline): {false_positives_count}")
    print(f"True Negatives (Both germline): {true_negative_count}")
    print(f"False Negatives in DV (Team germline, DV not): {false_negatives_pos_in_DV}")
    print(f"False Negatives not in DV: {false_negatives_pos_not_in_DV}")
    print(f"True Positives in DV (Both not germline): {true_positives_pos_in_DV}")
    print(f"True Positives not in DV: {true_positives_pos_not_in_DV}")

    total_false_negatives = false_negatives_pos_in_DV + false_negatives_pos_not_in_DV
    total_true_positives = true_positives_pos_in_DV + true_positives_pos_not_in_DV

    print(f"\nTotal False Negatives: {total_false_negatives}")
    print(f"Total True Positives: {total_true_positives}")

    # Save registers (use same timestamp as log)
    if false_positive_register:
        fp_df = pd.DataFrame(false_positive_register)
        fp_file = f"{OUTPUT_DIR}/false_positives_{timestamp}.csv"
        fp_df.to_csv(fp_file, index=False)
        print(f"\nFalse positives saved to: {fp_file}")

    if false_negative_register:
        fn_df = pd.DataFrame(false_negative_register)
        fn_file = f"{OUTPUT_DIR}/false_negatives_{timestamp}.csv"
        fn_df.to_csv(fn_file, index=False)
        print(f"False negatives saved to: {fn_file}")

    print(f"\nCompleted at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()