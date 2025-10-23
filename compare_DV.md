# Compare Germline Indels Documentation

## Overview
This document describes the process for comparing germline indels between two variant calling methods:
1. Team's custom indel caller (parquet format)
2. DeepVariant (VCF format)

## Key File Types

### Input Files
- **DeepVariant VCF** (`.vcf.gz`): Standard VCF format with all variant types (SNVs and indels)
  - Location: `/home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample*.vcf.gz`
  - Contains: All variants called by DeepVariant
  - Key fields: CHROM, POS, REF, ALT, QUAL, FILTER, GT, VAF, DP

- **Team Indel Caller Parquet** (`.parquet`): Custom format with pre-filtered indels
  - Location: `/home/ubuntu/data/teamb_indel-caller/results/sample*_indels_only.parquet`
  - Contains: Only indel calls with additional annotations
  - Key columns: chrom, pos, allele_seq, allele_type (REF/ALT), status (GERMLINE/SOMATIC), in_notdifficult

- **BED File** (`.bed.gz`): Genomic regions file
  - Location: `/home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz`
  - Contains: Complement of difficult regions (i.e., "not difficult" regions)
  - Format: chrom, start, end

### Output Files
- **Filtered DeepVariant Parquet**: Processed VCF data in parquet format for faster analysis
- **Comparison Results**: Statistics and matched variants

---

## Data Preparation Pipeline

### Step 1: Convert DeepVariant VCF to Parquet with Single-Indel Format and Region Annotation

**Purpose**: Extract indels from DeepVariant VCF, convert to single-indel format for accurate comparison, annotate with genomic region status, and save as parquet for efficient processing.

**Script**: `vcf_to_parquet_indels.py`

#### Command Used
```bash
python vcf_to_parquet_indels.py \
    --input /home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.vcf.gz \
    --output /home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.indels.parquet \
    --bed /home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz \
    --coding-bed /home/ubuntu/data/reference/beds/GRCh38_main3_cds.bed.gz
```

#### Processing Results
- **Input**: 6,277,647 total variants in VCF
- **Output**: 1,557,716 indels extracted (25% of variants)
- **File size**: 30MB parquet (vs 112MB gzipped VCF)
- **Processing time**: ~2 minutes including BED loading

#### Key Processing Steps

1. **BED Loading** (15 seconds)
   - Loaded 4,750,502 not-difficult regions into IntervalTree structure

2. **Indel Extraction**
   - Filtered for variants where `len(REF) != len(ALT)`
   - Created separate rows for multi-allelic sites
   - Result: 976,485 deletions (62.7%), 581,231 insertions (37.3%)

3. **Single-Indel Format Conversion**
   - Applied to 49,883 indels (3.2%) at multi-indel positions
   - 77,327 positions had multiple indels
   - 22,965 positions had mixed insertions/deletions
   - Ensures one allele is always 1bp anchor

4. **Region Annotation**
   - 371,755 indels in not-difficult regions (23.9%)
   - 1,185,961 indels in difficult regions (76.1%)

#### Complete List of Output Columns (21 total)

**Genomic location:**
1. `chrom`: Chromosome (e.g., chr1, chr2, ..., chrX, chrY)
2. `pos`: Position (1-based coordinate)

**Variant representation (after conversion):**
3. `ref`: Reference allele in single-indel format
4. `alt`: Alternate allele in single-indel format

**Original values (before conversion):**
5. `original_ref`: Original REF from VCF (null if no conversion applied)
6. `original_alt`: Original ALT from VCF (null if no conversion applied)
7. `conversion_applied`: Boolean - True if single-indel conversion was applied

**Quality and filter information:**
8. `qual`: Quality score from VCF
9. `filter`: Filter status - **Critical field**
   - **"PASS"** (981,985 = 63%): Real variant calls
   - **"RefCall"** (575,731 = 37%): Homozygous reference (0/0), NOT real variants

**Indel characteristics:**
10. `indel_type`: "insertion" or "deletion"
11. `indel_size`: Size of the indel in base pairs

**Genotype and depth information:**
12. `gt`: Genotype (0/0, 0/1, 1/1, ./.)
13. `gq`: Genotype quality
14. `dp`: Read depth
15. `ad`: Allelic depths (ref,alt counts)
16. `vaf`: Variant allele fraction
17. `pl`: Phred-scaled genotype likelihoods

**Region annotation:**
18. `in_notdifficult`: Boolean - True if in not-difficult region (23.9% are True)

**Multi-allelic flags:**
19. `is_multi_allelic`: Boolean - True if position has multiple alternate alleles
20. `has_multi_indels`: Boolean - True if position has multiple indels
21. `has_mixed_indel_types`: Boolean - True if position has both insertions and deletions

#### Single-Indel Format Verification

After conversion, **100% of indels** follow the correct format:
- **Insertions**: REF = 1bp anchor, ALT > 1bp (e.g., `T` → `TAA`)
- **Deletions**: REF > 1bp, ALT = 1bp anchor (e.g., `TAA` → `T`)
- **No problematic cases**: No SNVs or complex variants

Example conversions at multi-indel positions:
```
chr1:181583  Original: REF=CGGGGG, ALT=CG  →  Converted: REF=CGGGG, ALT=C
chr1:1102614 Original: REF=TA, ALT=TAA     →  Converted: REF=T, ALT=TA
chr1:791554  Original: REF=GGAATGGAATC, ALT=GGAATCGAATGGAATC  →  Converted: REF=G, ALT=GGAATC
```

#### Important Notes for Step 2

1. **Filter for real variants**: Use `filter == 'PASS'` to exclude RefCall entries
2. **Expected comparison counts** (in not-difficult regions):
   - DeepVariant PASS indels: ~233,000 (estimated)
   - Team's germline indels: 189,751
   - Revised ratio: ~1.23x (not 1.96x)
3. **Direct comparison possible**: Both datasets now use same single-indel format
4. **Performance**: No BED loading needed - just pandas boolean filtering

### Step 2: Compare Indels Called Germline By Our Pipeline to Those Called Germline By DV

**Objective**: Compare indels called as germline by the team's pipeline with those called as germline by DeepVariant, focusing on overlap statistics in not-difficult regions.
What do we most care about?:
1. False Positives: How many indels called somatic in Team B that were called Germline by DV
2. False Negatives: How many indels called Germline by Team B that were called somatic by DV (in that DV didn't call them germline)
3. True Negatives: We can also count the number of places where Team B and DV agree on Germline
4. True Positives: We can also count the number of places where Team B and DV agree on somatic (in that DV didn't call them germline)


Algo approach
- 1 & 3 require checking
- 2 & 4 can be counted by positions not in DV plus positions in DV but didn't match.
Make the dfs smaller
- team_b_df = team_b_df after dropping all rows in not difficult == FALSE (need to keep ref and germline and somatic)
- DV_df = DV_df after dropping all rows in not difficult == FALSE and filter != PASS (only keep PASS in not difficult)

Start with position checking
- Makes sense to do one pass
Counters:
- false_positives_count = 0
- true_negative_count = 0    
- false_negatives_pos_in_DV = 0
- false_negatives_pos_not_in_DV = 0
- true_positives_pos_in_DV = 0
- true_positives_pos_not_in_DV = 0

false_positive_register_df initiate empty - chrom, pos, ref, alt
false_negative_register_df initiate empty - chrom, pos, ref, alt

- Positions need checking are:
   - Team B positions where a row exists that is:
      - type is alt, somatic or germline, not in difficult region
   - DV positions with PASS and not in difficult region
   - Check for set of positions that overlap
   - Count not in DV at all:
      - false_negatives_pos_not_in_DV = count of germline alts (count all germline alts not just pos) where pos not in DV
      - true_positives_pos_not_in_DV = count of somatic alts (count all somatic alts not just pos) where pos not in DV
- for positions that overlap:
   - get the rows from Team B df
   - get the ref for the pos from Team B df
   - get the DV rows for the pos
   - for the Team B alt rows for the pos:
      - if the ref - alt pair match one of the DV ref - alt pairs for the position 
         - if alt has status somatic
            false_positives_count ++
            Add the false positive to the register
         - elseif alt has status germline 
            true_negative_count ++
         - else error
      - else
         - if alt has status somatic
            true_positives_pos_in_DV ++
         - elseif alt has status germline
            false_negatives_pos_in_DV ++
            Add the false negative to the register
         - else error
- save the registers
- print out final results for counters

#### Rationale for Denominators

**For False Positives (FP):**
- Denominator: Total Team somatic calls (FP + TP)
- Shows: What % of Team B's somatic calls were incorrectly classified (actually germline per DV)

**For True Negatives (TN):**
- Denominator: Total Team germline calls (TN + FN)
- Shows: What % of Team B's germline calls were correctly identified

**For False Negatives (FN):**
- Denominator: Total Team germline calls
- Shows: What % of Team B's germline calls were not found in DV

**For True Positives (TP):**
- Denominator: Total Team somatic calls
- Shows: What % of Team B's somatic calls were correctly identified

The percentages measure **Team B's accuracy** using Team B's own classifications as the baseline. This helps evaluate:
- **Somatic precision**: How accurate are Team B's somatic calls?
- **Germline precision**: How accurate are Team B's germline calls?

### Step 3: Check if any Indels Called Somatic By Our Pipeline are Called Germline By DV



### Configuration

Edit the global variables at the top of the script:

```python
# Location of called indels parquet files
DATA_DIR = "/home/ubuntu/data/teamb_indel-caller/results"

# Sample parquet file to analyze
SAMPLE_FILE = "sample1_indels_only.parquet"

# DeepVariant VCF file
VCF_FILE = "/home/ubuntu/data/teamb_indel-caller/deep_variant/GTEx-sample1.vcf.gz"

# BED file for not-difficult regions
BED_FILE = "/home/ubuntu/data/reference/beds/GRCh38_difficult_comp.bed.gz"
```

### Running the Script

```bash
python compare_germline_indels.py
```

### Output Interpretation

The script provides three main sections:

1. **Parquet Statistics**:
   - Total germline ALT indels
   - Filtering results for not-difficult regions
   - Breakdown by chromosome
   - Size distribution (single-base vs multi-base)

2. **VCF Statistics**:
   - Total indels in VCF
   - Filtering results for not-difficult regions
   - Multi-allelic site counts

3. **Comparison Summary**:
   - Count comparison
   - Ratio calculation
   - Difference in variant counts

### Current Results

After filtering for not-difficult regions:
- **Team's parquet**: 189,751 germline indels (24% of unfiltered)
- **DeepVariant ALL**: 371,755 total indels (24% of unfiltered)
- **DeepVariant PASS only**: ~233,000 estimated (real variants only)

**Important**: The initial 1.96x ratio includes RefCall entries (homozygous reference). When comparing only real variants (PASS filter), the ratio is approximately 1.23x.

---

## Performance Considerations

### BED File Loading
- Loading 4.7M regions takes significant time (~10-15 seconds)
- Consider pre-filtering VCF to avoid repeated BED loading

### Memory Usage
- Full VCF loading can consume significant memory
- Parquet format is more memory-efficient

### Optimization Strategies
1. Pre-filter and convert VCF to parquet (one-time cost)
2. Use indexed files (tabix for VCF)
3. Process chromosomes separately for large datasets

---

## Next Steps

1. **Create VCF to Parquet converter** (Step 1 above)
2. **Position-level comparison**: Match specific indel positions
3. **Genotype analysis**: Compare germline vs somatic calls
4. **Quality metrics**: Analyze VAF, depth, and quality distributions
5. **Concordance analysis**: Identify variants found by both methods

---

## Notes

- DeepVariant includes all variants (germline + somatic), while parquet pre-filters for germline
- Multi-allelic sites with mixed variant types (SNV + indel) are rare (~0.06% of positions)
- The 1.96x ratio is consistent before and after region filtering, suggesting systematic differences in calling sensitivity