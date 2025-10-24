#!/usr/bin/env python3
"""
Utility to lookup a specific position in Team B and DeepVariant parquet files
Also shows DeepVariant indels within ±10bp for proximity matching analysis
Usage: python lookup_position.py <chromosome> <position>
   or: python lookup_position.py <chromosome>,<position>
Example: python lookup_position.py chr1 12345
Example: python lookup_position.py chr1,12345
Example: python util/lookup_position.py chr1,12345
"""

import pandas as pd
import sys
from pathlib import Path

# Configuration (same as compare_DV.py)
TEAM_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/results"
TEAM_SAMPLE_FILE = "sample1_indels_only.parquet"
DV_DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/deep_variant"
DV_SAMPLE_FILE = "GTEx-sample1.indels.parquet"

# Proximity window for searching nearby DV indels (same as compare_DV_proximity.py)
PROXIMITY_WINDOW = 25

def main():
    # Parse arguments - support both space and comma separated
    if len(sys.argv) == 2 and ',' in sys.argv[1]:
        # Handle comma-separated format: chr1,12345
        try:
            chrom, pos_str = sys.argv[1].split(',')
            pos = int(pos_str)
        except (ValueError, IndexError):
            print(f"Error: Invalid format '{sys.argv[1]}'")
            print("Usage: python lookup_position.py <chromosome>,<position>")
            print("   or: python lookup_position.py <chromosome> <position>")
            sys.exit(1)
    elif len(sys.argv) == 3:
        # Handle space-separated format: chr1 12345
        chrom = sys.argv[1]
        try:
            pos = int(sys.argv[2])
        except ValueError:
            print(f"Error: Position must be an integer, got '{sys.argv[2]}'")
            sys.exit(1)
    else:
        print("Usage: python lookup_position.py <chromosome>,<position>")
        print("   or: python lookup_position.py <chromosome> <position>")
        print("Example: python lookup_position.py chr1,12345")
        print("Example: python lookup_position.py chr1 12345")
        sys.exit(1)

    print(f"Looking up {chrom}:{pos}")
    print("=" * 80)

    # Load Team B data
    team_path = Path(TEAM_DATA_DIR) / TEAM_SAMPLE_FILE
    if not team_path.exists():
        print(f"ERROR: Team B file not found: {team_path}")
        sys.exit(1)

    print("\nLoading Team B data...")
    team_df = pd.read_parquet(team_path)

    # Filter for position
    team_rows = team_df[(team_df['chrom'] == chrom) & (team_df['pos'] == pos)]

    print(f"\n### TEAM B ROWS ({len(team_rows)} rows) ###")
    if len(team_rows) > 0:
        # Display columns in a useful order
        cols = ['chrom', 'pos', 'allele_type', 'allele_seq', 'status',
                'allele_depth', 'total_depth', 'near_germ', 'in_coding', 'in_notdifficult']
        # Only show columns that exist
        cols = [c for c in cols if c in team_rows.columns]
        print(team_rows[cols].to_string(index=False))

        # Show REF and ALT alleles
        ref_rows = team_rows[team_rows['allele_type'] == 'REF']
        alt_rows = team_rows[team_rows['allele_type'] == 'ALT']

        if len(ref_rows) > 0:
            print(f"\nREF allele: {ref_rows.iloc[0]['allele_seq']}")

        if len(alt_rows) > 0:
            print(f"ALT alleles ({len(alt_rows)}):")
            for _, alt_row in alt_rows.iterrows():
                status = alt_row['status']
                seq = alt_row['allele_seq']
                depth = alt_row.get('allele_depth', 'N/A')
                total = alt_row.get('total_depth', 'N/A')
                vaf = depth/total if depth != 'N/A' and total != 'N/A' and total > 0 else 'N/A'
                vaf_str = f"{vaf:.3f}" if vaf != 'N/A' else 'N/A'
                print(f"  {seq} ({status}, depth={depth}/{total}, VAF={vaf_str})")
    else:
        print("No rows found at this position")

    # Load DeepVariant data
    dv_path = Path(DV_DATA_DIR) / DV_SAMPLE_FILE
    if not dv_path.exists():
        print(f"ERROR: DeepVariant file not found: {dv_path}")
        sys.exit(1)

    print("\n" + "-" * 80)
    print("\nLoading DeepVariant data...")
    dv_df = pd.read_parquet(dv_path)

    # Filter for position
    dv_rows = dv_df[(dv_df['chrom'] == chrom) & (dv_df['pos'] == pos)]

    print(f"\n### DEEPVARIANT ROWS ({len(dv_rows)} rows) ###")
    if len(dv_rows) > 0:
        # Display key columns
        cols = ['chrom', 'pos', 'ref', 'alt', 'filter', 'gt', 'qual', 'gq',
                'dp', 'ad', 'vaf', 'in_notdifficult', 'original_ref', 'original_alt']
        # Only show columns that exist
        cols = [c for c in cols if c in dv_rows.columns]
        print(dv_rows[cols].to_string(index=False))

        # Summarize
        for _, dv_row in dv_rows.iterrows():
            ref = dv_row['ref']
            alt = dv_row['alt']
            filt = dv_row['filter']
            gt = dv_row.get('gt', 'N/A')
            vaf = dv_row.get('vaf', 'N/A')

            print(f"\n{ref} -> {alt}:")
            print(f"  Filter: {filt}")
            print(f"  Genotype: {gt}")
            vaf_str = f"{vaf:.3f}" if vaf != 'N/A' and not pd.isna(vaf) else 'N/A'
            print(f"  VAF: {vaf_str}")

            if 'original_ref' in dv_row and pd.notna(dv_row['original_ref']):
                print(f"  Original: {dv_row['original_ref']} -> {dv_row['original_alt']}")
    else:
        print("No rows found at this position")

    # Search for DeepVariant indels within proximity window
    print("\n" + "-" * 80)
    print(f"\n### DEEPVARIANT INDELS WITHIN ±{PROXIMITY_WINDOW}bp (proximity matching) ###")

    # Get positions within the window, excluding the exact position we already showed
    low_pos = pos - PROXIMITY_WINDOW
    high_pos = pos + PROXIMITY_WINDOW

    dv_nearby = dv_df[
        (dv_df['chrom'] == chrom) &
        (dv_df['pos'] >= low_pos) &
        (dv_df['pos'] <= high_pos) &
        (dv_df['pos'] != pos)  # Exclude exact position already shown above
    ].sort_values('pos')

    if len(dv_nearby) > 0:
        print(f"Found {len(dv_nearby)} DV indels in range {chrom}:{low_pos}-{high_pos} (excluding exact position):\n")

        # Group by position for clearer display
        for position in sorted(dv_nearby['pos'].unique()):
            pos_rows = dv_nearby[dv_nearby['pos'] == position]
            distance = position - pos
            distance_str = f"+{distance}" if distance > 0 else str(distance)

            print(f"Position {chrom}:{position} (offset {distance_str}bp):")
            for _, row in pos_rows.iterrows():
                ref = row['ref']
                alt = row['alt']
                filt = row['filter']
                vaf = row.get('vaf', 'N/A')
                vaf_str = f"{vaf:.3f}" if vaf != 'N/A' and not pd.isna(vaf) else 'N/A'

                deletion_size = len(ref) - len(alt) if len(ref) > len(alt) else 0
                insertion_size = len(alt) - len(ref) if len(alt) > len(ref) else 0

                if deletion_size > 0:
                    var_type = f"{deletion_size}bp deletion"
                elif insertion_size > 0:
                    var_type = f"{insertion_size}bp insertion"
                else:
                    var_type = "substitution"

                print(f"  {ref:10s} -> {alt:10s} ({var_type}, filter={filt}, VAF={vaf_str})")
    else:
        print(f"No DV indels found within ±{PROXIMITY_WINDOW}bp of position {pos}")

    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()