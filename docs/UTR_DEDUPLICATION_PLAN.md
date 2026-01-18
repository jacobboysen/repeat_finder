# UTR TE Hit Deduplication Plan

**STATUS: IMPLEMENTED** (2026-01-15)

## Overview

This document outlines the deduplication approach for 5'UTR and 3'UTR TE hit analyses. Deduplication has been implemented and reveals **significant inflation** in previous raw hit counts.

## Results Summary

| Region | Raw Hits | Unique Hits | Removed | Duplication Rate |
|--------|----------|-------------|---------|------------------|
| **Exon** | 1,722,273 | 1,500,719 | 221,554 | **12.86%** |
| **5'UTR** | 1,903,660 | 1,255,113 | 648,547 | **34.07%** |
| **3'UTR** | 2,573,926 | 1,370,202 | 1,203,724 | **46.77%** |

**CRITICAL**: Nearly half of all 3'UTR hits were duplicates from shared isoform regions!

## Problem Statement

When analyzing UTR TE hits, overlapping transcript isoforms can cause the same genomic region to be BLASTed multiple times. This inflates hit counts when the same TE sequence matches the same genomic location across different transcript isoforms.

**Root cause discovered:** GFF files use comma-separated Parent values (`Parent=FBtr001,FBtr002`) for UTRs shared across transcripts. Many genes have dozens of transcript isoforms sharing the same UTR sequence due to:
- Alternative polyadenylation (3'UTR) - most impactful
- Alternative promoters/TSS (5'UTR)
- Alternative splicing that doesn't affect UTR boundaries

## Deduplication Strategy

### Deduplication Key

Same approach as exon deduplication:
```
(fbgn, chrom, genomic_start, genomic_end, te_id, te_start, te_end)
```

**Key decisions:**
1. **Exact coordinate matching** - no tolerance/fuzzy matching
2. **Include TE coordinates** - same genomic region hitting different TE positions = separate hits
3. **Gene-level deduplication** - deduplicate within genes, not across genes

### Coordinate Conversion

UTR FASTA headers contain location information:
```
>FBtr0070000 type=three_prime_UTR; loc=X:19967461..19968479; parent=FBgn0031081; ...
```

To convert BLAST query coordinates to genomic:
```python
def query_to_genomic_coords(qstart, qend, utr_meta):
    """Convert query-relative coordinates to genomic coordinates."""
    if utr_meta['strand'] == '+':
        genomic_start = utr_meta['start'] + qstart - 1
        genomic_end = utr_meta['start'] + qend - 1
    else:
        # Minus strand: query coords are relative to reverse complement
        genomic_end = utr_meta['end'] - qstart + 1
        genomic_start = utr_meta['end'] - qend + 1
    return genomic_start, genomic_end
```

**Note:** Need to extract strand from GFF since UTR FASTA headers only have coordinates.

## Implementation Plan

### Option A: Single Unified Script

Create `scripts/deduplicate_utr_te_hits.py` that handles both 3'UTR and 5'UTR:

```bash
python scripts/deduplicate_utr_te_hits.py \
    --blast-file results/genome_wide_all_3utrs.tsv \
    --utr-type 3utr \
    --output-dir results/3utr_deduplicated

python scripts/deduplicate_utr_te_hits.py \
    --blast-file results/genome_wide_all_5utrs.tsv \
    --utr-type 5utr \
    --output-dir results/5utr_deduplicated
```

**Pros:**
- Single script to maintain
- Consistent behavior between UTR types

**Cons:**
- More complex argument handling

### Option B: Adapt Exon Script

Generalize `deduplicate_exon_te_hits.py` to handle any sequence type with genomic coordinates:

```bash
python scripts/deduplicate_te_hits.py \
    --blast-file results/genome_wide_all_3utrs.tsv \
    --metadata-source gff \  # or fasta-header
    --output-dir results/3utr_deduplicated
```

**Pros:**
- One script for all deduplication
- Cleaner codebase

**Cons:**
- More refactoring needed

### Option C: Separate Scripts (Recommended)

Create dedicated scripts for each UTR type:
- `scripts/deduplicate_3utr_te_hits.py`
- `scripts/deduplicate_5utr_te_hits.py`

**Pros:**
- Simpler implementation
- Can handle UTR-specific quirks independently
- Follows existing pattern (exon script is separate)

**Cons:**
- Some code duplication (mitigate with shared utility functions)

## Required Inputs

### For 3'UTR Deduplication

1. **BLAST results**: `results/genome_wide_all_3utrs.tsv`
2. **UTR metadata** (one of):
   - Parse from FASTA headers: `data/references/dmel_3utr.fasta`
   - Parse from GFF: `data/references/dmel-all-r6.66.gff`
3. **FBtr → FBgn mapping**: From GFF or existing mapping files

### For 5'UTR Deduplication

1. **BLAST results**: `results/genome_wide_all_5utrs.tsv`
2. **UTR metadata**: Same sources as 3'UTR
3. **FBtr → FBgn mapping**: Same as 3'UTR

## Expected Outputs

For each UTR type:
- `{utr}_te_hits_deduplicated.tsv` - Unique hits with annotations
- `gene_te_summary_deduplicated.tsv` - Gene-level aggregation
- `deduplication_stats.json` - Statistics
- `original_vs_deduplicated.tsv` - Per-gene comparison

## Integration with Existing Pipeline

After deduplication:

1. **Update `analyze_integrated_te_enrichment.py`** to use deduplicated UTR data
2. **Update any downstream analyses** that currently use raw UTR hit counts
3. **Update HTML visualizations** if needed

## Validation

1. **Sanity checks:**
   - Unique hits ≤ raw hits for all genes
   - Genomic coordinates are valid (within chromosome bounds)
   - TE coordinates are valid (within TE length)

2. **Spot-check high-duplication genes:**
   - Verify these are genes with many isoforms/alternative UTRs
   - Confirm duplicates are truly identical (same genomic + TE position)

3. **Compare duplication rates:**
   - 3'UTR vs 5'UTR vs exon
   - Per-gene correlation with isoform count

## Questions for Discussion

1. **Which implementation option to use?** (A, B, or C)

2. **Metadata source preference:**
   - Parse from FASTA headers (simpler but strand info may be missing)
   - Parse from GFF (more complete but requires GFF parsing)

3. **Output location:**
   - Alongside existing results: `results/genome_wide_all_3utrs_deduplicated.tsv`
   - New directory: `results/3utr_deduplicated/`
   - Subdirectory: `results/deduplicated/3utr/`

4. **Should we re-run integrated analysis** after UTR deduplication?

5. **Priority:**
   - Deduplicate both UTR types now?
   - Just 3'UTR first (more data, more impact)?

## Estimated Impact

Based on exon analysis:
- **Exon duplication rate:** 12.86%
- **Expected UTR duplication rate:** Likely lower (5-10%) since UTRs are more isoform-specific

However, certain gene categories may have higher rates:
- Genes with alternative polyadenylation (long vs short 3'UTR)
- Genes with alternative promoters (different 5'UTRs)
- Genes with many transcript isoforms

## Timeline

Implementation should be straightforward given the existing exon deduplication code:
1. Create script(s) - adapt from exon version
2. Run on 3'UTR data
3. Run on 5'UTR data
4. Update integrated analysis
5. Validate results
