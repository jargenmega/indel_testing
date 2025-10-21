#!/usr/bin/env python3
"""
Compare germline indels between DeepVariant VCF and parquet files
"""

import pandas as pd
import gzip
import os
from pathlib import Path
from intervaltree import IntervalTree
from collections import defaultdict

# ====================
# FILE PATHS AND CONFIGURATION
# ====================

# Location of called indels parquet files
DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/results"

# Sample parquet file to analyze
SAMPLE_FILE = "sample1_indels_only.parquet"

# DeepVariant VCF file
VCF_FILE = "/home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.vcf.gz"

# BED file for not-difficult regions (complement to difficult regions)
BED_FILE = "/home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz"


def load_bed_intervals(bed_file):
    """
    Load BED regions into interval trees for fast O(log n) lookup.
    Returns: {chrom: IntervalTree} for not-difficult regions.
    """
    trees = defaultdict(IntervalTree)
    print(f"  Loading BED file: {bed_file}")

    open_func = gzip.open if bed_file.endswith('.gz') else open
    region_count = 0

    with open_func(bed_file, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            # Add 'chr' prefix if not present (for matching VCF format)
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
            start = int(parts[1])
            end = int(parts[2])
            trees[chrom][start:end] = True
            region_count += 1

    print(f"  Loaded {region_count:,} not-difficult regions")
    return trees


def position_in_bed(chrom, pos, bed_trees):
    """
    Check if position falls within any BED region.
    Returns True if in not-difficult region.
    """
    if chrom not in bed_trees:
        return False
    # Query returns overlapping intervals
    return len(bed_trees[chrom][pos]) > 0


def count_parquet_germline_indels(parquet_path):
    """
    Load parquet file and count germline indels.

    Since this is an indels_only.parquet file, all entries are indels.
    Germline indels are rows where:
    - status == 'GERMLINE'
    - allele_type == 'ALT'

    Note: We don't filter on near_germ because all germline ALTs have
    near_germ=True by definition (100% of germline ALTs are near germline sites).
    """
    print(f"Loading parquet file: {parquet_path}")

    if not os.path.exists(parquet_path):
        print(f"  ERROR: File not found: {parquet_path}")
        return None

    # Load the parquet file
    df = pd.read_parquet(parquet_path)
    print(f"  Total rows loaded: {len(df):,}")

    # Show data structure
    print(f"  Columns: {list(df.columns)}")

    # Count germline ALT indels in not-difficult regions
    germline_indels_all = df[
        (df['status'] == 'GERMLINE') &
        (df['allele_type'] == 'ALT')
    ]

    # Filter for not-difficult regions
    germline_indels = germline_indels_all[
        germline_indels_all['in_notdifficult'] == True
    ]

    print(f"  Germline ALT indels found: {len(germline_indels_all):,}")
    print(f"  After filtering for not-difficult regions: {len(germline_indels):,}")

    germline_indel_count = len(germline_indels)

    # Additional statistics
    if germline_indel_count > 0:
        # Count by chromosome
        chrom_counts = germline_indels['chrom'].value_counts().sort_index()
        print("\n  Germline indels by chromosome (top 10):")
        for chrom, count in chrom_counts.head(10).items():
            print(f"    {chrom}: {count:,}")

        # Get unique positions
        germline_positions = germline_indels[['chrom', 'pos']].drop_duplicates()
        print(f"\n  Unique positions with germline indels: {len(germline_positions):,}")

        # Classify by indel size (based on allele_seq length)
        single_base = germline_indels[germline_indels['allele_seq'].str.len() == 1]
        multi_base = germline_indels[germline_indels['allele_seq'].str.len() > 1]

        print(f"\n  Indel size breakdown (by ALT sequence length):")
        print(f"    Single-base indels (1bp): {len(single_base):,}")
        print(f"    Multi-base indels (>1bp): {len(multi_base):,}")

    return germline_indel_count


def count_vcf_indels(vcf_path, bed_trees):
    """
    Load VCF file and count indel variants in not-difficult regions.

    Indels are variants where:
    - REF and ALT have different lengths
    - Position falls within not-difficult regions (BED file)
    """
    print(f"\nLoading VCF file: {vcf_path}")

    if not os.path.exists(vcf_path):
        print(f"  ERROR: File not found: {vcf_path}")
        return None

    indel_count = 0
    indel_count_unfiltered = 0
    snv_count = 0
    total_variants = 0
    multi_allelic_count = 0
    filtered_out = 0

    # Open the gzipped VCF file
    with gzip.open(vcf_path, 'rt') as vcf:
        for line in vcf:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Parse VCF line
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])  # Convert to int for BED lookup
            ref = fields[3]
            alt = fields[4]

            total_variants += 1

            # Check for multi-allelic sites (comma in ALT field)
            if ',' in alt:
                multi_allelic_count += 1
                # For multi-allelic sites, count each alt allele separately
                alt_alleles = alt.split(',')
                for alt_allele in alt_alleles:
                    if len(ref) != len(alt_allele):
                        # This is an indel
                        indel_count_unfiltered += 1
                        # Check if in not-difficult region
                        if position_in_bed(chrom, pos, bed_trees):
                            indel_count += 1
                        else:
                            filtered_out += 1
                    else:
                        snv_count += 1
            else:
                # Single alternate allele
                if len(ref) != len(alt):
                    # This is an indel
                    indel_count_unfiltered += 1
                    # Check if in not-difficult region
                    if position_in_bed(chrom, pos, bed_trees):
                        indel_count += 1
                    else:
                        filtered_out += 1
                else:
                    # This is a SNV
                    snv_count += 1

    print(f"  Total variants loaded: {total_variants:,}")
    print(f"  Multi-allelic sites: {multi_allelic_count:,}")
    print(f"  Total indels found: {indel_count_unfiltered:,}")
    print(f"  After filtering for not-difficult regions: {indel_count:,}")
    print(f"  Filtered out (in difficult regions): {filtered_out:,}")

    return indel_count


def main():
    """Main function to compare germline indels in not-difficult regions."""

    print("=" * 80)
    print("Germline Indel Comparison: Parquet vs DeepVariant VCF")
    print("(Filtered for not-difficult regions)")
    print("=" * 80)

    # Construct full file paths
    parquet_path = os.path.join(DATA_DIR, SAMPLE_FILE)
    vcf_path = VCF_FILE

    print(f"\nConfiguration:")
    print(f"  Parquet file: {parquet_path}")
    print(f"  VCF file: {vcf_path}")
    print(f"  BED file (not-difficult regions): {BED_FILE}")

    # Load BED intervals for VCF filtering
    print("\nLoading not-difficult regions from BED file...")
    bed_trees = load_bed_intervals(BED_FILE)

    print("\n" + "=" * 80)
    print("PART 1: Counting Germline Indels in Parquet File")
    print("(Filtering on in_notdifficult column)")
    print("=" * 80)

    parquet_germline_count = count_parquet_germline_indels(parquet_path)

    print("\n" + "=" * 80)
    print("PART 2: Counting Indels in DeepVariant VCF")
    print("(Filtering using BED regions)")
    print("=" * 80)

    vcf_indel_count = count_vcf_indels(vcf_path, bed_trees)

    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY (NOT-DIFFICULT REGIONS ONLY)")
    print("=" * 80)

    if parquet_germline_count is not None:
        print(f"Parquet - Germline indels (in_notdifficult=True): {parquet_germline_count:,}")

    if vcf_indel_count is not None:
        print(f"DeepVariant - Total indels (in BED regions): {vcf_indel_count:,}")

    if parquet_germline_count is not None and vcf_indel_count is not None:
        difference = abs(parquet_germline_count - vcf_indel_count)
        print(f"\nDifference: {difference:,} indels")

        if parquet_germline_count > 0:
            ratio = vcf_indel_count / parquet_germline_count
            print(f"Ratio (VCF/Parquet): {ratio:.2f}")

    print("\n" + "=" * 80)
    print("Analysis complete!")


if __name__ == "__main__":
    main()