#!/usr/bin/env python3
"""
Convert DeepVariant VCF to parquet format with only indels.
Optionally filters for not-difficult regions.
"""

import gzip
import pandas as pd
import argparse
from pathlib import Path
from intervaltree import IntervalTree
from collections import defaultdict
import sys

def load_bed_intervals(bed_file):
    """Load BED regions into interval trees for fast O(log n) lookup."""
    trees = defaultdict(IntervalTree)
    print(f"Loading BED file: {bed_file}")

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
            # Add 'chr' prefix if not present
            if not chrom.startswith('chr'):
                chrom = f'chr{chrom}'
            start = int(parts[1])
            end = int(parts[2])
            trees[chrom][start:end] = True
            region_count += 1

    print(f"Loaded {region_count:,} regions")
    return trees


def position_in_bed(chrom, pos, bed_trees):
    """Check if position falls within any BED region."""
    if bed_trees is None:
        return True  # No filtering if no BED file provided
    if chrom not in bed_trees:
        return False
    return len(bed_trees[chrom][pos]) > 0


def convert_to_single_indel_format(ref, alt):
    """
    Convert from combined VCF format to single-indel format for comparison.

    In combined format: All ALTs share one REF (e.g., REF=AGATCTG)
    In single format: Each indel has minimal REF/ALT pair

    Rules:
    - First base is the anchor (all indels happen after it)
    - Insertions: REF=anchor, ALT=anchor+inserted_bases
    - Deletions: REF=anchor+deleted_bases, ALT=anchor
    """
    anchor = ref[0]  # First base is always the anchor

    if len(alt) > len(ref):
        # INSERTION: bases inserted after anchor
        # Example: REF=AGATCTG, ALT=ATGATCTG → REF=A, ALT=AT
        inserted_length = len(alt) - len(ref)
        inserted_bases = alt[1:1+inserted_length]
        new_ref = anchor
        new_alt = anchor + inserted_bases
        return new_ref, new_alt

    elif len(alt) < len(ref):
        # DELETION: bases deleted after anchor
        # Example: REF=AGATCTG, ALT=ACTG → REF=AGAT, ALT=A
        deleted_length = len(ref) - len(alt)
        deleted_bases = ref[1:1+deleted_length]
        new_ref = anchor + deleted_bases
        new_alt = anchor
        return new_ref, new_alt

    else:
        # Same length - no conversion needed (shouldn't happen for indels)
        return ref, alt


def parse_vcf_to_dataframe(vcf_path, bed_trees=None, max_variants=None):
    """
    Parse VCF file and extract indels into a DataFrame.

    Parameters:
    - vcf_path: Path to VCF file (can be gzipped)
    - bed_trees: Optional BED interval trees for region annotation
    - max_variants: Optional limit on number of variants to process

    Returns: DataFrame with indel information including in_notdifficult column
    """

    # Storage for parsed data
    data = []

    # Counters
    total_variants = 0
    total_indels = 0
    indels_in_notdifficult = 0
    multi_allelic = 0
    mixed_type = 0

    print(f"Processing VCF: {vcf_path}")

    with gzip.open(vcf_path, 'rt') as vcf:
        for line in vcf:
            # Skip header lines
            if line.startswith('#'):
                continue

            total_variants += 1

            # Progress indicator
            if total_variants % 100000 == 0:
                print(f"  Processed {total_variants:,} variants, found {total_indels:,} indels...")

            # Optional limit
            if max_variants and total_variants > max_variants:
                print(f"  Reached limit of {max_variants:,} variants")
                break

            # Parse VCF fields
            fields = line.strip().split('\t')
            if len(fields) < 10:  # Need at least FORMAT and one sample
                continue

            chrom = fields[0]
            pos = int(fields[1])
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_val = fields[6]
            info = fields[7]
            format_keys = fields[8].split(':')
            sample_vals = fields[9].split(':')

            # Create FORMAT dictionary for easy access
            format_dict = dict(zip(format_keys, sample_vals))

            # Handle multi-allelic sites
            alt_alleles = alt.split(',') if alt != '.' else []

            if len(alt_alleles) > 1:
                multi_allelic += 1

            # Check if position has any indels and mixed types
            has_snv = False
            has_insertion = False
            has_deletion = False
            indel_alleles = []  # Store indel alleles to process
            indel_types = []    # Track whether each is insertion or deletion

            for alt_allele in alt_alleles:
                if alt_allele != '.':
                    if len(ref) == len(alt_allele):
                        has_snv = True
                    else:
                        if len(alt_allele) > len(ref):
                            has_insertion = True
                            indel_types.append('insertion')
                        else:
                            has_deletion = True
                            indel_types.append('deletion')
                        indel_alleles.append(alt_allele)

            # Track mixed variant types
            if has_snv and (has_insertion or has_deletion):
                mixed_type += 1

            # Skip if no indels at this position
            if not (has_insertion or has_deletion):
                continue

            # Determine if conversion is needed
            has_multi_indels = len(indel_alleles) > 1
            has_mixed_indel_types = has_insertion and has_deletion

            # Conversion is needed if:
            # 1. Multiple indels at same position (especially mixed ins/del)
            # 2. Single deletion with extended REF (check if REF is longer than needed)
            needs_conversion = has_multi_indels

            # Check if position is in not-difficult region (only once per position)
            in_notdifficult = position_in_bed(chrom, pos, bed_trees) if bed_trees else True

            # Process each indel allele
            for i, alt_allele in enumerate(alt_alleles):
                # Skip non-indels
                if alt_allele == '.' or len(ref) == len(alt_allele):
                    continue

                total_indels += 1

                if in_notdifficult:
                    indels_in_notdifficult += 1

                # Extract genotype information
                gt = format_dict.get('GT', './.')
                gq = format_dict.get('GQ', '.')
                dp = format_dict.get('DP', '.')
                ad = format_dict.get('AD', '.')
                vaf = format_dict.get('VAF', '.')
                pl = format_dict.get('PL', '.')

                # Parse VAF for multi-allelic sites
                if vaf != '.' and ',' in vaf:
                    vaf_values = vaf.split(',')
                    vaf = vaf_values[i] if i < len(vaf_values) else '.'

                # Parse AD for multi-allelic sites
                if ad != '.' and ',' in ad:
                    ad_values = ad.split(',')
                    # AD format: ref_depth, alt1_depth, alt2_depth, ...
                    ref_ad = ad_values[0] if len(ad_values) > 0 else '.'
                    alt_ad = ad_values[i+1] if i+1 < len(ad_values) else '.'
                    ad = f"{ref_ad},{alt_ad}"

                # Determine indel type
                if len(alt_allele) > len(ref):
                    indel_type = 'insertion'
                    indel_size = len(alt_allele) - len(ref)
                else:
                    indel_type = 'deletion'
                    indel_size = len(ref) - len(alt_allele)

                # Apply conversion if needed
                original_ref = ref
                original_alt = alt_allele
                converted_ref = ref
                converted_alt = alt_allele
                conversion_applied = False

                if needs_conversion:
                    converted_ref, converted_alt = convert_to_single_indel_format(ref, alt_allele)
                    # Check if conversion actually changed anything
                    if converted_ref != ref or converted_alt != alt_allele:
                        conversion_applied = True

                # Store data
                data.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': converted_ref,
                    'alt': converted_alt,
                    'original_ref': original_ref if conversion_applied else None,
                    'original_alt': original_alt if conversion_applied else None,
                    'conversion_applied': conversion_applied,
                    'qual': qual,
                    'filter': filter_val,
                    'indel_type': indel_type,
                    'indel_size': indel_size,
                    'gt': gt,
                    'gq': gq,
                    'dp': dp,
                    'ad': ad,
                    'vaf': vaf,
                    'pl': pl,
                    'in_notdifficult': in_notdifficult,
                    'is_multi_allelic': len(alt_alleles) > 1,
                    'has_multi_indels': has_multi_indels,
                    'has_mixed_indel_types': has_mixed_indel_types
                })

    print(f"\nProcessing complete:")
    print(f"  Total variants: {total_variants:,}")
    print(f"  Total indels: {total_indels:,}")
    print(f"  Multi-allelic sites: {multi_allelic:,}")
    print(f"  Mixed SNV/indel sites: {mixed_type:,}")

    if bed_trees:
        print(f"  Indels in not-difficult regions: {indels_in_notdifficult:,} ({indels_in_notdifficult/total_indels*100:.1f}%)")
        print(f"  Indels in difficult regions: {total_indels - indels_in_notdifficult:,} ({(total_indels - indels_in_notdifficult)/total_indels*100:.1f}%)")

    print(f"  Total indels in parquet: {len(data):,}")

    return pd.DataFrame(data)


def main():
    parser = argparse.ArgumentParser(description='Convert VCF indels to parquet format with region annotation')
    parser.add_argument('--input', '-i', required=True, help='Input VCF file (can be .vcf.gz)')
    parser.add_argument('--output', '-o', required=True, help='Output parquet file')
    parser.add_argument('--bed', '-b', help='BED file for annotating not-difficult regions (adds in_notdifficult column)')
    parser.add_argument('--max-variants', '-m', type=int, help='Maximum variants to process (for testing)')

    args = parser.parse_args()

    # Load BED file if provided
    bed_trees = None
    if args.bed:
        if not Path(args.bed).exists():
            print(f"ERROR: BED file not found: {args.bed}", file=sys.stderr)
            return 1
        bed_trees = load_bed_intervals(args.bed)

    # Check input file
    if not Path(args.input).exists():
        print(f"ERROR: Input VCF not found: {args.input}", file=sys.stderr)
        return 1

    # Process VCF
    df = parse_vcf_to_dataframe(args.input, bed_trees, args.max_variants)

    if df.empty:
        print("WARNING: No indels found or all filtered out", file=sys.stderr)
        return 1

    # Convert numeric columns
    numeric_cols = ['pos', 'qual', 'indel_size', 'gq', 'dp', 'vaf']
    for col in numeric_cols:
        if col in df.columns:
            # Convert to numeric, keeping '.' as NaN
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Save to parquet
    print(f"\nSaving to parquet: {args.output}")
    df.to_parquet(args.output, compression='snappy')

    # Print summary statistics
    print("\nDataFrame summary:")
    print(f"  Shape: {df.shape}")
    print(f"  Chromosomes: {df['chrom'].nunique()}")
    print(f"  Unique positions: {df.groupby(['chrom', 'pos']).size().shape[0]}")

    print("\nIndel type breakdown:")
    print(df['indel_type'].value_counts())

    print("\nIndel size distribution:")
    print(df['indel_size'].describe())

    print("\nFilter values:")
    print(df['filter'].value_counts().head(10))

    if 'in_notdifficult' in df.columns and bed_trees:
        print("\nRegion annotation:")
        print(f"  In not-difficult regions: {df['in_notdifficult'].sum():,} ({df['in_notdifficult'].mean()*100:.1f}%)")
        print(f"  In difficult regions: {(~df['in_notdifficult']).sum():,} ({(~df['in_notdifficult']).mean()*100:.1f}%)")

    # Report on conversion statistics
    if 'conversion_applied' in df.columns:
        print("\nConversion statistics:")
        conversions = df['conversion_applied'].sum()
        print(f"  Conversions applied: {conversions:,} ({conversions/len(df)*100:.1f}%)")

        if conversions > 0:
            print(f"\n  Breakdown of converted positions:")
            print(f"    With multiple indels: {df[df['conversion_applied']]['has_multi_indels'].sum():,}")
            print(f"    With mixed indel types: {df[df['conversion_applied']]['has_mixed_indel_types'].sum():,}")

            # Show examples of converted positions
            print(f"\n  Example conversions (first 5):")
            converted_df = df[df['conversion_applied']].head()
            for idx, row in converted_df.iterrows():
                print(f"    {row['chrom']}:{row['pos']}")
                print(f"      Original: REF={row['original_ref']}, ALT={row['original_alt']}")
                print(f"      Converted: REF={row['ref']}, ALT={row['alt']}")

    # Report on multi-indel positions
    if 'has_multi_indels' in df.columns:
        print("\nMulti-indel positions:")
        multi_indel_positions = df[df['has_multi_indels']].groupby(['chrom', 'pos']).size()
        print(f"  Positions with multiple indels: {len(multi_indel_positions):,}")

        if 'has_mixed_indel_types' in df.columns:
            mixed_type_positions = df[df['has_mixed_indel_types']].groupby(['chrom', 'pos']).size()
            print(f"  Positions with mixed ins/del: {len(mixed_type_positions):,}")

    print(f"\nOutput saved to: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())