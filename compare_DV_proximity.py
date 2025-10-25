#!/usr/bin/env python3
"""
Compare Team B indel calls with DeepVariant calls using PROXIMITY matching
Focus on germline/somatic classification agreement using position proximity
rather than exact allele matching
Note that there are some efficency optimisations that could be made to this code now that the actual sequences are not being compared, but it is pretty fast and making the changes would change the alogo compared to the exact matching version making them harder to compare directly.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from datetime import datetime
import bisect

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

# Output file prefix - set to descriptive name like "s1_all", "s1_coding", "s2_all", etc.
OUTPUT_PREFIX = "s1_coding"  # Example: "s1_all" will create files like "s1_all_compare_DV_proximity_log_20251024_123456.txt"

# Proximity window for matching (in base pairs)
PROXIMITY_WINDOW = 10  # Will check for DV indels within +/- 10bp of Team B position

# Testing mode - set to None to process all chromosomes, or list specific ones
TESTING = None #, ['chr1'], ['chr1', 'chr2', 'chr21']

# Optional filter - set to True to only include coding regions
FILTER_CODING_ONLY = True  # If True, only includes positions where in_coding == True

def has_nearby_dv(pos, dv_positions_sorted, proximity):
    """Check if any DV position exists within proximity window using binary search.
    O(log n) time complexity instead of O(window_size)."""
    lo = pos - proximity
    hi = pos + proximity
    i = bisect.bisect_left(dv_positions_sorted, lo)
    return i < len(dv_positions_sorted) and dv_positions_sorted[i] <= hi

def convert_to_single_indel_format(ref, alt):
    """Convert from combined VCF format to single-indel format"""
    anchor = ref[0]
    if len(alt) > len(ref):
        # INSERTION: bases inserted after anchor
        # Example: REF=AGATCTG, ALT=ATGATCTG → REF=A, ALT=AT
        inserted_length = len(alt) - len(ref)
        inserted_bases = alt[1:1+inserted_length]
        return anchor, anchor + inserted_bases
    elif len(alt) < len(ref):
        # DELETION: bases deleted after anchor
        # Example: REF=AGATCTG, ALT=ACTG → REF=AGAT, ALT=A
        deleted_length = len(ref) - len(alt)
        deleted_bases = ref[1:1+deleted_length]
        return anchor + deleted_bases, anchor
    else:
        return ref, alt

def main():
    # Set up logging
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # Add prefix before "compare_" so files sort by run
    prefix_part = f"{OUTPUT_PREFIX}_" if OUTPUT_PREFIX else ""
    log_file = f"{OUTPUT_DIR}/{prefix_part}compare_DV_proximity_log_{timestamp}.txt"
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
    print("  Team B filters: in_notdifficult == True")
    print("  AND removing: somatic ALTs where near_germ == True (keeping all REF rows)")

    # Note: REF rows can have status='SOMATIC' when the reference allele is nearly absent
    # (e.g., 1 re1ad out of 30), indicating a germline variant has replaced the reference.
    # We must keep ALL REF rows regardless of status since we need them for variant comparison.
    # Filter: Keep all in not-difficult regions, but remove only somatic ALT alleles near germline sites
    # Keep ALL REF rows and germline ALTs regardless of near_germ status
    team_filter = (team_df['in_notdifficult'] == True) & \
                  ~((team_df['allele_type'] == 'ALT') & (team_df['status'] == 'SOMATIC') & (team_df['near_germ'] == True))
    if FILTER_CODING_ONLY:
        team_filter = team_filter & (team_df['in_coding'] == True)
        print("  Additional filter: in_coding == True (coding regions only)")

    team_df = team_df[team_filter].copy()
    print(f"  Team B after filtering: {len(team_df):,} rows")

    # Note that there are some with PASS but genotype 0/0; we accept those.
    print("  DeepVariant filters: in_notdifficult == True AND filter == 'PASS'")
    dv_filter = (dv_df['in_notdifficult'] == True) & (dv_df['filter'] == 'PASS')
    if FILTER_CODING_ONLY:
        dv_filter = dv_filter & (dv_df['in_coding'] == True)
        print("  Additional filter: in_coding == True (coding regions only)")

    dv_df = dv_df[dv_filter].copy()
    print(f"  DeepVariant after filtering: {len(dv_df):,} rows")

    # Initialize counters
    false_positives_count = 0
    true_negative_count = 0
    false_negatives_count = 0
    true_positives_count = 0

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

        # Get DV positions as a sorted list for efficient binary search
        dv_positions_sorted = sorted(dv_chrom['pos'].unique())

        # For proximity matching, check each Team B position for nearby DV indels
        team_pos_near_dv = set()
        team_pos_not_near_dv = set()

        for pos in team_alt_positions:
            # Check if any DV position exists within proximity window using binary search
            if has_nearby_dv(pos, dv_positions_sorted, PROXIMITY_WINDOW):
                team_pos_near_dv.add(pos)
            else:
                team_pos_not_near_dv.add(pos)

        print(f"  Team ALT positions: {len(team_alt_positions):,}")
        print(f"  DV positions: {len(dv_positions_sorted):,}")
        print(f"  Team positions with DV within ±{PROXIMITY_WINDOW}bp: {len(team_pos_near_dv):,}")
        print(f"  Team positions with no DV nearby: {len(team_pos_not_near_dv):,}")

        # Process positions with no DV nearby (no DV indel within proximity window)
        # These are positions where Team B calls an indel but DeepVariant doesn't call anything nearby
        # - GERMLINE calls here are false negatives (Team thinks germline, DV doesn't call)
        # - SOMATIC calls here are true positives (both agree not germline)
        # We collect false negatives for analysis but just count true positives (less problematic)
        for pos in team_pos_not_near_dv:
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

            # Process each Team B ALT
            for _, team_alt in team_alts.iterrows():
                alt_seq = team_alt['allele_seq']
                team_status = team_alt['status']

                # Apply conversion if neither ref nor alt is 1bp (anchor)
                # Single-indel format should always have either ref=1bp or alt=1bp
                if len(ref_seq) > 1 and len(alt_seq) > 1:
                    converted_ref, converted_alt = convert_to_single_indel_format(ref_seq, alt_seq)
                else:
                    converted_ref, converted_alt = ref_seq, alt_seq

                # Categorize (no DV nearby means DV doesn't call germline in the region)
                if team_status == 'GERMLINE':
                    false_negatives_count += 1
                    false_negative_register.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': converted_ref,
                        'alt': converted_alt,
                        'dv_nearby': False,
                        'allele_depth': team_alt['allele_depth'],
                        'total_depth': team_alt['total_depth']
                    })
                elif team_status == 'SOMATIC':
                    true_positives_count += 1

        # Process positions with DV nearby (DV indel exists within proximity window)
        # Using proximity matching: if ANY DV indel exists within ±PROXIMITY_WINDOW bp, we consider it a match
        # - Team GERMLINE + DV nearby = true negative (both agree germline, within proximity)
        # - Team SOMATIC + DV nearby = false positive (Team somatic, DV germline nearby)
        # We collect false positives for analysis
        for pos in team_pos_near_dv:
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

            # Process each Team B ALT
            # Since we're doing proximity matching, ANY DV indel nearby means DV calls germline in the region
            for _, team_alt in team_alts.iterrows():
                alt_seq = team_alt['allele_seq']
                team_status = team_alt['status']

                # Apply conversion if neither ref nor alt is 1bp (anchor)
                # Single-indel format should always have either ref=1bp or alt=1bp
                if len(ref_seq) > 1 and len(alt_seq) > 1:
                    converted_ref, converted_alt = convert_to_single_indel_format(ref_seq, alt_seq)
                else:
                    converted_ref, converted_alt = ref_seq, alt_seq

                # Categorize based on proximity (DV indel exists within window = germline region)
                if team_status == 'GERMLINE':
                    true_negative_count += 1  # Both agree germline (Team germline, DV has indel nearby)
                elif team_status == 'SOMATIC':
                    false_positives_count += 1  # Team says somatic, DV has germline indel nearby
                    false_positive_register.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': converted_ref,
                        'alt': converted_alt,
                        'allele_depth': team_alt['allele_depth'],
                        'total_depth': team_alt['total_depth']
                    })

    # Print results
    print("\nProcessing complete!")
    print("\n=== RESULTS (PROXIMITY MATCHING) ===")
    print(f"Note: Using proximity window of ±{PROXIMITY_WINDOW}bp for matching")
    print("Note: Counts are of ALT alleles and do not take into account allele depth")
    print("Note: Percentages use Team B calls as denominator to measure Team B accuracy")
    print("Note: DV indel within proximity window = germline region")
    print("-" * 60)

    # Total Team somatic and germline calls
    # These are used as denominators to measure Team B's accuracy
    # - For FP%: out of all Team somatic calls, how many were wrong?
    # - For TN%: out of all Team germline calls, how many were correct?
    total_team_somatic = false_positives_count + true_positives_count
    total_team_germline = true_negative_count + false_negatives_count
    total_team_calls = total_team_somatic + total_team_germline

    # Print with counts and percentages
    print(f"False Positives (Team somatic, DV germline nearby): {false_positives_count}", end="")
    if total_team_somatic > 0:
        print(f" ({false_positives_count/total_team_somatic*100:.1f}% of Team somatic)", end="")
    print()

    print(f"True Negatives (Both germline): {true_negative_count}", end="")
    if total_team_germline > 0:
        print(f" ({true_negative_count/total_team_germline*100:.1f}% of Team germline)", end="")
    print()

    print(f"False Negatives (Team germline, no DV nearby): {false_negatives_count}", end="")
    if total_team_germline > 0:
        print(f" ({false_negatives_count/total_team_germline*100:.1f}% of Team germline)", end="")
    print()

    print(f"True Positives (Team somatic, no DV nearby): {true_positives_count}", end="")
    if total_team_somatic > 0:
        print(f" ({true_positives_count/total_team_somatic*100:.1f}% of Team somatic)", end="")
    print()

    # Overall accuracy metrics
    if total_team_calls > 0:
        correct = true_negative_count + true_positives_count
        print(f"\n--- OVERALL METRICS ---")
        print(f"Total Team ALT calls analyzed: {total_team_calls:,}")
        print(f"Team somatic calls: {total_team_somatic:,} ({total_team_somatic/total_team_calls*100:.1f}%)")
        print(f"Team germline calls: {total_team_germline:,} ({total_team_germline/total_team_calls*100:.1f}%)")
        print(f"Agreement with DV (proximity): {correct:,}/{total_team_calls:,} ({correct/total_team_calls*100:.1f}%)")

        if total_team_somatic > 0:
            print(f"Somatic precision: {true_positives_count/total_team_somatic*100:.1f}%")
        if total_team_germline > 0:
            print(f"Germline precision: {true_negative_count/total_team_germline*100:.1f}%")

    # Save registers (use same timestamp as log)
    if false_positive_register:
        fp_df = pd.DataFrame(false_positive_register)
        fp_file = f"{OUTPUT_DIR}/{prefix_part}false_positives_proximity_{timestamp}.csv"
        fp_df.to_csv(fp_file, index=False)
        print(f"\nFalse positives saved to: {fp_file}")

    if false_negative_register:
        fn_df = pd.DataFrame(false_negative_register)
        fn_file = f"{OUTPUT_DIR}/{prefix_part}false_negatives_proximity_{timestamp}.csv"
        fn_df.to_csv(fn_file, index=False)
        print(f"False negatives saved to: {fn_file}")

    print(f"\nCompleted at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()