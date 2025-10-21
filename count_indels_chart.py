#!/usr/bin/env python3
"""
Fast indel analysis script for CRAM/BAM files with visualization.
Counts:
1. Number of reads with each indel count
2. For each indel count >1: how many sum to zero
"""

import pysam
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
from intervaltree import IntervalTree
import gzip
from datetime import datetime
import os

# Global configuration
MAX_READS = 10000000  # Set to an integer to cap the number of reads processed, None for no limit
BED_FILE = "/home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz"  # Filter to reads starting in these regions
MAX_INDEL_DISTANCE = 12  # Maximum distance (bp) between indels to be considered "close"

def load_bed_intervals(bed_file):
    """
    Load BED regions into interval trees for fast O(log n) lookup.
    Returns: {chrom: IntervalTree}
    """
    trees = defaultdict(IntervalTree)

    print(f"Loading BED file: {bed_file}...", file=sys.stderr)

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
            start = int(parts[1])
            end = int(parts[2])
            trees[chrom][start:end] = True
            region_count += 1

    print(f"Loaded {region_count:,} regions into interval trees", file=sys.stderr)
    return trees

def read_starts_in_bed(read, bed_trees):
    """
    Check if the read start position falls within any BED region using interval tree.
    O(log n) lookup time.
    """
    if read.reference_name not in bed_trees:
        return False

    read_start = read.reference_start
    # Query returns overlapping intervals - we check if read start falls in any
    return len(bed_trees[read.reference_name][read_start]) > 0

def parse_cigar_indels(cigar_tuples):
    """
    Extract indel information from CIGAR tuples.
    Returns list of indel sizes (positive=insertion, negative=deletion)
    """
    indels = []
    for op, length in cigar_tuples:
        if op == 1:  # Insertion
            indels.append(length)
        elif op == 2:  # Deletion
            indels.append(-length)
    return indels

def parse_cigar_close_indels(cigar_tuples):
    """
    Extract only "close" indels from CIGAR tuples.
    An indel is only included if it's within MAX_INDEL_DISTANCE bp of another indel.
    Returns list of indel sizes (positive=insertion, negative=deletion)
    """
    if cigar_tuples is None:
        return []

    close_indels = []
    last_indel_size = None
    last_indel_added = False
    distance_since_indel = 0

    for op, length in cigar_tuples:
        if op == 1 or op == 2:  # Insertion or Deletion
            indel_size = length if op == 1 else -length

            if last_indel_size is not None and distance_since_indel <= MAX_INDEL_DISTANCE:
                # Close to previous indel
                if not last_indel_added:
                    close_indels.append(last_indel_size)
                    last_indel_added = True
                close_indels.append(indel_size)
                last_indel_size = indel_size
                last_indel_added = True
            else:
                # First indel or too far from previous
                last_indel_size = indel_size
                last_indel_added = False

            distance_since_indel = 0

        elif op in [0, 7, 8, 4]:  # Match, =, X, or Soft clip (consume query)
            distance_since_indel += length

    return close_indels

def main():
    if len(sys.argv) < 2:
        print("Usage: python count_indels_chart.py <cram_file> [reference.fasta] [output_chart.png]")
        sys.exit(1)

    cram_file = sys.argv[1]
    reference = sys.argv[2] if len(sys.argv) > 2 else None

    # Create descriptive output filename with timestamp
    if len(sys.argv) > 3:
        output_chart = sys.argv[3]
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sample_name = os.path.basename(cram_file).replace('.cram', '').replace('.bam', '')
        output_chart = f"plots/indel_distribution_frameshift_analysis_{sample_name}_{timestamp}.png"

    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)

    # Load BED regions into interval trees
    bed_trees = load_bed_intervals(BED_FILE)

    # Statistics counters
    indel_count_distribution = defaultdict(int)  # {num_indels: count}
    sum_zero_by_count = defaultdict(int)  # {num_indels: count_that_sum_to_zero}
    total_reads = 0
    reads_in_bed = 0
    reads_outside_bed = 0

    # Counter for printing first 10 non-zero sum examples
    nonzero_examples_printed = 0
    max_examples = 10

    # Counter for printing reads with >3 close indels
    many_indels_printed = 0
    max_many_indels_examples = 10

    print(f"Processing {cram_file}...", file=sys.stderr)

    # Open the alignment file
    save_kwargs = {}
    if reference:
        save_kwargs['reference_filename'] = reference

    with pysam.AlignmentFile(cram_file, "rc", **save_kwargs) as bam:
        for read in bam:
            total_reads += 1

            # Check if we've reached the cap
            if MAX_READS is not None and total_reads > MAX_READS:
                print(f"Reached read cap of {MAX_READS:,}. Stopping.", file=sys.stderr)
                break

            # Progress indicator
            if total_reads % 1000000 == 0:
                print(f"Processed {total_reads:,} reads...", file=sys.stderr)

            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Check if read starts in BED region
            if not read_starts_in_bed(read, bed_trees):
                reads_outside_bed += 1
                continue

            reads_in_bed += 1

            # Parse only close indels from CIGAR
            indels = parse_cigar_close_indels(read.cigartuples)
            num_indels = len(indels)

            # Count reads by number of indels
            indel_count_distribution[num_indels] += 1

            # Print examples with >3 close indels
            if num_indels > 3 and many_indels_printed < max_many_indels_examples:
                indel_sum = sum(indels)
                print(f"\n>3 close indels example #{many_indels_printed + 1}:")
                print(f"  Chromosome: {read.reference_name}")
                print(f"  Start: {read.reference_start}")
                print(f"  Read name: {read.query_name}")
                print(f"  CIGAR: {read.cigarstring}")
                print(f"  Close indels: {indels} (count={num_indels}, sum={indel_sum})")
                many_indels_printed += 1

            # Check if indels sum to zero
            if num_indels > 1:
                indel_sum = sum(indels)
                if indel_sum == 0:
                    sum_zero_by_count[num_indels] += 1
                else:
                    # Print first 10 examples of non-zero sum reads
                    if nonzero_examples_printed < max_examples:
                        print(f"\nNon-zero sum example #{nonzero_examples_printed + 1}:")
                        print(f"  Chromosome: {read.reference_name}")
                        print(f"  Start: {read.reference_start}")
                        print(f"  Read name: {read.query_name}")
                        print(f"  CIGAR: {read.cigarstring}")
                        print(f"  Indels: {indels} (sum={indel_sum})")
                        nonzero_examples_printed += 1

    print(f"\nTotal reads processed: {total_reads:,}", file=sys.stderr)
    print(f"Reads in BED regions: {reads_in_bed:,}", file=sys.stderr)
    print(f"Reads outside BED regions: {reads_outside_bed:,}\n", file=sys.stderr)

    # Output text results
    print("=" * 70)
    print("INDEL COUNT DISTRIBUTION WITH SUM=0 ANALYSIS")
    print("(Only reads starting in BED regions)")
    print("=" * 70)
    print(f"{'Num Indels':<15} {'Total Count':<20} {'Sum=0 Count':<20} {'% Sum=0':<15}")
    print("-" * 70)

    max_indels = max(indel_count_distribution.keys()) if indel_count_distribution else 0
    for i in range(max_indels + 1):
        count = indel_count_distribution[i]
        sum_zero = sum_zero_by_count[i]
        pct_sum_zero = (sum_zero / count * 100) if count > 0 else 0
        print(f"{i:<15} {count:<20,} {sum_zero:<20,} {pct_sum_zero:<15.4f}%")

    # Create chart (excluding 0 and 1 indel counts)
    print(f"\nGenerating chart: {output_chart}", file=sys.stderr)

    # Filter data for >1 indels
    x_values = []
    total_counts = []
    sum_zero_counts = []

    for i in range(2, max_indels + 1):
        if indel_count_distribution[i] > 0:
            x_values.append(i)
            total_counts.append(indel_count_distribution[i])
            sum_zero_counts.append(sum_zero_by_count[i])

    if x_values:
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Plot 1: Bar chart with total and sum=0
        x_pos = range(len(x_values))
        width = 0.35

        ax1.bar([p - width/2 for p in x_pos], total_counts, width,
                label='Total Reads', alpha=0.8, color='steelblue')
        ax1.bar([p + width/2 for p in x_pos], sum_zero_counts, width,
                label='Sum = 0 (No Frameshift)', alpha=0.8, color='orange')

        ax1.set_xlabel('Number of Indels per Read', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax1.set_title('Indel Distribution: Total vs Sum=0', fontsize=14, fontweight='bold')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_values)
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)

        # Plot 2: Percentage that sum to zero
        percentages = [(sz / tc * 100) if tc > 0 else 0
                       for sz, tc in zip(sum_zero_counts, total_counts)]

        ax2.bar(x_pos, percentages, color='green', alpha=0.7)
        ax2.set_xlabel('Number of Indels per Read', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Percentage Sum = 0', fontsize=12, fontweight='bold')
        ax2.set_title('Percentage of Reads Where Indels Sum to Zero (No Frameshift)',
                     fontsize=14, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_values)
        ax2.grid(axis='y', alpha=0.3)
        ax2.set_ylim(0, 100)

        # Add percentage labels on bars
        for i, (xp, pct) in enumerate(zip(x_pos, percentages)):
            if pct > 0:
                ax2.text(xp, pct + 1, f'{pct:.1f}%', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        plt.savefig(output_chart, dpi=300, bbox_inches='tight')
        print(f"Chart saved to: {output_chart}", file=sys.stderr)
    else:
        print("No reads with >1 indel found. No chart generated.", file=sys.stderr)

if __name__ == "__main__":
    main()
