#!/usr/bin/env python3
"""
Analyzes base quality around indels for reads in specified BED regions.
For reads with average quality > threshold:
- For insertions: quality of inserted bases
- For insertions/deletions: quality of 5bp before and after the indel
"""

import pysam
import sys
from collections import defaultdict
from intervaltree import IntervalTree
import gzip
import numpy as np

# Global configuration
MAX_READS = 10000000  # Set to an integer to cap the number of reads processed, None for no limit
BED_FILE = "/home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz"  # Filter to reads starting in these regions
MIN_AVG_QUALITY = 30  # Minimum average read quality (Q30)
CONTEXT_SIZE = 5  # Number of bases before/after indel to analyze

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
    return len(bed_trees[read.reference_name][read_start]) > 0

def get_average_quality(qualities):
    """Calculate average quality score."""
    if not qualities or len(qualities) == 0:
        return 0
    return np.mean(qualities)

def analyze_indel_qualities(read):
    """
    Analyze quality scores around indels in a read.
    Returns dict with insertion and deletion quality information.
    """
    if read.cigartuples is None or read.query_qualities is None:
        return None

    seq = read.query_sequence
    quals = read.query_qualities

    results = {
        'insertions': [],  # Each: {'insert_qual': avg, 'before_qual': avg, 'after_qual': avg, 'length': N}
        'deletions': []    # Each: {'before_qual': avg, 'after_qual': avg, 'length': N}
    }

    query_pos = 0  # Position in the read sequence

    for op, length in read.cigartuples:
        if op == 1:  # Insertion
            # Get quality of inserted bases
            insert_start = query_pos
            insert_end = query_pos + length
            insert_quals = quals[insert_start:insert_end]

            # Get quality of CONTEXT_SIZE bases before insertion
            before_start = max(0, query_pos - CONTEXT_SIZE)
            before_quals = quals[before_start:query_pos]

            # Get quality of CONTEXT_SIZE bases after insertion
            after_end = min(len(quals), insert_end + CONTEXT_SIZE)
            after_quals = quals[insert_end:after_end]

            results['insertions'].append({
                'insert_qual': get_average_quality(insert_quals),
                'before_qual': get_average_quality(before_quals),
                'after_qual': get_average_quality(after_quals),
                'length': length
            })

            query_pos += length

        elif op == 2:  # Deletion
            # For deletions, get quality of bases before and after
            # Get quality of CONTEXT_SIZE bases before deletion
            before_start = max(0, query_pos - CONTEXT_SIZE)
            before_quals = quals[before_start:query_pos]

            # Get quality of CONTEXT_SIZE bases after deletion point
            after_end = min(len(quals), query_pos + CONTEXT_SIZE)
            after_quals = quals[query_pos:after_end]

            results['deletions'].append({
                'before_qual': get_average_quality(before_quals),
                'after_qual': get_average_quality(after_quals),
                'length': length
            })
            # Deletions don't consume query sequence

        elif op in [0, 7, 8]:  # Match, =, X (consume query)
            query_pos += length
        # Operations 3 (N), 4 (S), 5 (H), 6 (P) handled implicitly

    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python indel_quality_analysis.py <cram_file> [reference.fasta]")
        sys.exit(1)

    cram_file = sys.argv[1]
    reference = sys.argv[2] if len(sys.argv) > 2 else None

    # Load BED regions into interval trees
    bed_trees = load_bed_intervals(BED_FILE)

    # Quality score counters (count by integer quality values)
    insertion_insert_quals = defaultdict(int)  # Quality of inserted bases
    insertion_before_quals = defaultdict(int)  # Quality before insertions
    insertion_after_quals = defaultdict(int)   # Quality after insertions

    deletion_before_quals = defaultdict(int)   # Quality before deletions
    deletion_after_quals = defaultdict(int)    # Quality after deletions

    # Statistics
    total_reads = 0
    reads_in_bed = 0
    reads_high_quality = 0
    reads_analyzed = 0
    total_insertions = 0
    total_deletions = 0

    print(f"Processing {cram_file}...", file=sys.stderr)
    print(f"Minimum average quality: Q{MIN_AVG_QUALITY}", file=sys.stderr)
    print(f"Context size: {CONTEXT_SIZE}bp before/after indel\n", file=sys.stderr)

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
                continue

            reads_in_bed += 1

            # Check average read quality
            if read.query_qualities is None:
                continue

            avg_qual = get_average_quality(read.query_qualities)
            if avg_qual < MIN_AVG_QUALITY:
                continue

            reads_high_quality += 1

            # Analyze indels
            indel_data = analyze_indel_qualities(read)
            if indel_data is None:
                continue

            if indel_data['insertions'] or indel_data['deletions']:
                reads_analyzed += 1

            # Process insertions
            for ins in indel_data['insertions']:
                total_insertions += 1
                # Round to nearest integer and bin
                insertion_insert_quals[round(ins['insert_qual'])] += 1
                insertion_before_quals[round(ins['before_qual'])] += 1
                insertion_after_quals[round(ins['after_qual'])] += 1

            # Process deletions
            for dels in indel_data['deletions']:
                total_deletions += 1
                deletion_before_quals[round(dels['before_qual'])] += 1
                deletion_after_quals[round(dels['after_qual'])] += 1

    # Print statistics
    print(f"\n{'='*70}", file=sys.stderr)
    print(f"PROCESSING STATISTICS", file=sys.stderr)
    print(f"{'='*70}", file=sys.stderr)
    print(f"Total reads processed: {total_reads:,}", file=sys.stderr)
    print(f"Reads in BED regions: {reads_in_bed:,}", file=sys.stderr)
    print(f"Reads with avg quality >= Q{MIN_AVG_QUALITY}: {reads_high_quality:,}", file=sys.stderr)
    print(f"Reads with indels analyzed: {reads_analyzed:,}", file=sys.stderr)
    print(f"Total insertions found: {total_insertions:,}", file=sys.stderr)
    print(f"Total deletions found: {total_deletions:,}\n", file=sys.stderr)

    # Print results
    print("=" * 90)
    print("INSERTION QUALITY ANALYSIS")
    print("=" * 90)
    print(f"\nQuality of Inserted Bases (n={total_insertions:,}):")
    print(f"{'Quality (Q)':<15} {'Count':<20} {'Percentage':<15}")
    print("-" * 50)
    for q in sorted(insertion_insert_quals.keys()):
        count = insertion_insert_quals[q]
        pct = (count / total_insertions * 100) if total_insertions > 0 else 0
        print(f"Q{q:<14} {count:<20,} {pct:<15.2f}%")

    print(f"\nQuality {CONTEXT_SIZE}bp Before Insertions (n={total_insertions:,}):")
    print(f"{'Quality (Q)':<15} {'Count':<20} {'Percentage':<15}")
    print("-" * 50)
    for q in sorted(insertion_before_quals.keys()):
        count = insertion_before_quals[q]
        pct = (count / total_insertions * 100) if total_insertions > 0 else 0
        print(f"Q{q:<14} {count:<20,} {pct:<15.2f}%")

    print(f"\nQuality {CONTEXT_SIZE}bp After Insertions (n={total_insertions:,}):")
    print(f"{'Quality (Q)':<15} {'Count':<20} {'Percentage':<15}")
    print("-" * 50)
    for q in sorted(insertion_after_quals.keys()):
        count = insertion_after_quals[q]
        pct = (count / total_insertions * 100) if total_insertions > 0 else 0
        print(f"Q{q:<14} {count:<20,} {pct:<15.2f}%")

    print("\n" + "=" * 90)
    print("DELETION QUALITY ANALYSIS")
    print("=" * 90)
    print(f"\nQuality {CONTEXT_SIZE}bp Before Deletions (n={total_deletions:,}):")
    print(f"{'Quality (Q)':<15} {'Count':<20} {'Percentage':<15}")
    print("-" * 50)
    for q in sorted(deletion_before_quals.keys()):
        count = deletion_before_quals[q]
        pct = (count / total_deletions * 100) if total_deletions > 0 else 0
        print(f"Q{q:<14} {count:<20,} {pct:<15.2f}%")

    print(f"\nQuality {CONTEXT_SIZE}bp After Deletions (n={total_deletions:,}):")
    print(f"{'Quality (Q)':<15} {'Count':<20} {'Percentage':<15}")
    print("-" * 50)
    for q in sorted(deletion_after_quals.keys()):
        count = deletion_after_quals[q]
        pct = (count / total_deletions * 100) if total_deletions > 0 else 0
        print(f"Q{q:<14} {count:<20,} {pct:<15.2f}%")

if __name__ == "__main__":
    main()
