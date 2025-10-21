#!/usr/bin/env python3
"""
Explore VCF structure to understand columns, INFO fields, and variant types
"""

import gzip
import sys
from collections import defaultdict

# Configuration
VCF_FILE = "/home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.vcf.gz"
MAX_EXAMPLES = 20  # Number of example rows to show

def main():
    print("=" * 80)
    print("VCF STRUCTURE EXPLORATION")
    print("=" * 80)
    print(f"\nAnalyzing: {VCF_FILE}\n")

    # Storage for analysis
    info_fields = set()
    format_fields = set()
    column_headers = []
    example_rows = []
    info_descriptions = {}
    format_descriptions = {}

    # Track positions with multiple variant types
    position_variants = defaultdict(list)  # {(chrom, pos): [variant_types]}
    mixed_position_lines = defaultdict(list)  # {(chrom, pos): [full_lines]}

    # Counters
    total_lines = 0
    header_lines = 0
    data_lines = 0

    with gzip.open(VCF_FILE, 'rt') as vcf:
        for line in vcf:
            total_lines += 1

            # Process header lines
            if line.startswith('##'):
                header_lines += 1
                # Extract INFO field descriptions
                if line.startswith('##INFO='):
                    # Parse INFO field definition
                    # Example: ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
                    if 'ID=' in line and 'Description=' in line:
                        id_start = line.index('ID=') + 3
                        id_end = line.index(',', id_start)
                        field_id = line[id_start:id_end]

                        desc_start = line.index('Description="') + 13
                        desc_end = line.rindex('"')
                        description = line[desc_start:desc_end]

                        info_descriptions[field_id] = description

                # Extract FORMAT field descriptions
                elif line.startswith('##FORMAT='):
                    if 'ID=' in line and 'Description=' in line:
                        id_start = line.index('ID=') + 3
                        id_end = line.index(',', id_start)
                        field_id = line[id_start:id_end]

                        desc_start = line.index('Description="') + 13
                        desc_end = line.rindex('"')
                        description = line[desc_start:desc_end]

                        format_descriptions[field_id] = description

                continue

            # Process column header line
            if line.startswith('#CHROM'):
                header_lines += 1
                column_headers = line.strip().split('\t')
                continue

            # Process data lines
            data_lines += 1
            fields = line.strip().split('\t')

            # Store first few examples
            if len(example_rows) < MAX_EXAMPLES:
                example_rows.append(fields)

            # Analyze this variant
            if len(fields) >= 8:
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                info = fields[6]
                format_field = fields[8] if len(fields) > 8 else ""

                # Track variant types at this position
                position_key = (chrom, pos)

                # Store the full line for this position
                mixed_position_lines[position_key].append(line.strip())

                # Check each ALT allele
                for alt_allele in alt.split(','):
                    if alt_allele == '.':
                        continue
                    if len(ref) != len(alt_allele):
                        position_variants[position_key].append('INDEL')
                    else:
                        position_variants[position_key].append('SNV')

                # Extract INFO field keys
                if info != '.':
                    for item in info.split(';'):
                        if '=' in item:
                            key = item.split('=')[0]
                            info_fields.add(key)
                        else:
                            info_fields.add(item)  # FLAG fields

                # Extract FORMAT field keys
                if format_field:
                    for key in format_field.split(':'):
                        format_fields.add(key)

            # Stop early for exploration
            if data_lines >= 10000:
                break

    # Print results
    print("FILE STATISTICS:")
    print(f"  Total lines processed: {total_lines:,}")
    print(f"  Header lines: {header_lines:,}")
    print(f"  Data lines analyzed: {data_lines:,}")

    print("\n" + "=" * 80)
    print("COLUMN STRUCTURE:")
    print("=" * 80)
    for i, col in enumerate(column_headers):
        print(f"  Column {i}: {col}")

    print("\n" + "=" * 80)
    print("INFO FIELDS FOUND (with descriptions):")
    print("=" * 80)
    for field in sorted(info_fields):
        desc = info_descriptions.get(field, "No description available")
        print(f"  {field}: {desc}")

    print("\n" + "=" * 80)
    print("FORMAT FIELDS FOUND (with descriptions):")
    print("=" * 80)
    for field in sorted(format_fields):
        desc = format_descriptions.get(field, "No description available")
        print(f"  {field}: {desc}")

    print("\n" + "=" * 80)
    print("EXAMPLE DATA ROWS:")
    print("=" * 80)

    # Show a few example rows
    for i, row in enumerate(example_rows[:5]):
        print(f"\nExample {i+1}:")
        for j, (col_name, value) in enumerate(zip(column_headers, row)):
            # Truncate long values for display
            if len(value) > 50:
                value = value[:47] + "..."
            print(f"  {col_name}: {value}")

        # Analyze variant type
        if len(row) >= 5:
            ref = row[3]
            alt = row[4]
            print(f"  Variant type: ", end="")
            types = []
            for alt_allele in alt.split(','):
                if alt_allele != '.':
                    if len(ref) != len(alt_allele):
                        types.append(f"INDEL ({ref}->{alt_allele})")
                    else:
                        types.append(f"SNV ({ref}->{alt_allele})")
            print(", ".join(types))

    print("\n" + "=" * 80)
    print("POSITIONS WITH MULTIPLE VARIANT TYPES:")
    print("=" * 80)

    # Check for positions with both INDELs and SNVs
    mixed_positions = []
    indel_only_positions = 0
    snv_only_positions = 0

    for pos_key, variants in position_variants.items():
        unique_types = set(variants)
        if 'INDEL' in unique_types and 'SNV' in unique_types:
            mixed_positions.append(pos_key)
        elif 'INDEL' in unique_types:
            indel_only_positions += 1
        elif 'SNV' in unique_types:
            snv_only_positions += 1

    print(f"\nSummary of position types (first 10,000 variants):")
    print(f"  Positions with only INDELs: {indel_only_positions:,}")
    print(f"  Positions with only SNVs: {snv_only_positions:,}")
    print(f"  Positions with both INDELs and SNVs: {len(mixed_positions):,}")

    # Print the actual VCF lines for mixed positions
    if mixed_positions:
        print("\n" + "=" * 80)
        print("VCF LINES FOR POSITIONS WITH BOTH INDELS AND SNVS:")
        print("=" * 80)

        for i, pos_key in enumerate(mixed_positions):
            chrom, pos = pos_key
            print(f"\nPosition {i+1}: {chrom}:{pos}")
            print("-" * 60)

            # Print all VCF lines for this position
            for line in mixed_position_lines[pos_key]:
                fields = line.split('\t')
                if len(fields) >= 5:
                    ref = fields[3]
                    alt = fields[4]
                    filter_val = fields[6] if len(fields) > 6 else "."

                    # Analyze each ALT allele
                    alt_types = []
                    for alt_allele in alt.split(','):
                        if alt_allele != '.':
                            if len(ref) != len(alt_allele):
                                alt_types.append(f"{alt_allele}(INDEL)")
                            else:
                                alt_types.append(f"{alt_allele}(SNV)")

                    print(f"  REF: {ref}, ALT: {', '.join(alt_types)}, FILTER: {filter_val}")

                    # Print full line (truncated if too long)
                    if len(line) > 150:
                        print(f"  Full line: {line[:147]}...")
                    else:
                        print(f"  Full line: {line}")

    # Show multi-allelic examples
    print("\n" + "=" * 80)
    print("MULTI-ALLELIC SITE EXAMPLES:")
    print("=" * 80)

    multi_count = 0
    for row in example_rows:
        if len(row) >= 5:
            alt = row[4]
            if ',' in alt:
                multi_count += 1
                if multi_count <= 3:
                    chrom = row[0]
                    pos = row[1]
                    ref = row[3]
                    print(f"\n  {chrom}:{pos}")
                    print(f"    REF: {ref}")
                    print(f"    ALT: {alt}")
                    alt_list = alt.split(',')
                    for alt_allele in alt_list:
                        if len(ref) != len(alt_allele):
                            print(f"      -> {alt_allele} (INDEL)")
                        else:
                            print(f"      -> {alt_allele} (SNV)")

if __name__ == "__main__":
    main()