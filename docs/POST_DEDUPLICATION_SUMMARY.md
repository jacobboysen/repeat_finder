# Post-Deduplication Analysis Summary

**Date:** 2026-01-15
**Status:** COMPLETED

## Overview

This document summarizes the impact of hit deduplication on TE analysis results. Deduplication removes duplicate hits that arise when overlapping transcript isoforms produce identical BLAST hits to the same genomic and TE positions.

---

## Deduplication Statistics

### Hit-Level Deduplication

| Region | Raw Hits | Unique Hits | Duplicates | Rate |
|--------|----------|-------------|------------|------|
| **3'UTR** | 2,573,926 | 1,370,202 | 1,203,724 | **46.77%** |
| **5'UTR** | 1,903,660 | 1,255,113 | 648,547 | **34.07%** |
| **Exon** | 1,722,273 | 1,500,719 | 221,554 | **12.86%** |

**Key insight:** 3'UTRs had the highest duplication rate (47%) because alternative polyadenylation creates many transcript isoforms sharing the same 3'UTR sequence.

### Shuffled Control Comparison

| Dataset | Raw Hits | Unique Hits | Duplication Rate |
|---------|----------|-------------|------------------|
| Real 3'UTR (10% sample) | 463,094 | 409,751 | **11.52%** |
| Shuffled (mean) | 214,442 | 214,437 | **~0%** |

**Enrichment change:**
- Old (raw): 2.16x
- New (deduplicated): **1.91x**
- Change: -11.5%

---

## Impact on Rankings

### Gene Rankings by TE Density

The **top 100 and bottom 100 gene rankings did NOT change significantly** because:

1. Rankings are based on **density** (hit_bp per kb of UTR), not raw hit count
2. Deduplication removes duplicate *hits*, but preserved hits retain their alignment length
3. Most top-ranked genes have few isoforms and short UTRs (density metric favors these)

### Genes Most Affected by Deduplication

These genes had the highest duplication rates but may not appear in top/bottom rankings:

| Gene (FBgn) | Raw Hits | Unique Hits | Duplicates | Rate |
|-------------|----------|-------------|------------|------|
| FBgn0033159 | 6,472 | 161 | 6,311 | **97.5%** |
| FBgn0003435 (Dscam1) | 9,534 | 905 | 8,629 | **90.5%** |
| FBgn0003345 (Sxl) | 9,219 | 871 | 8,348 | **90.6%** |
| FBgn0264562 | 11,082 | 1,635 | 9,447 | **85.2%** |
| FBgn0264001 | 12,533 | 2,432 | 10,101 | **80.6%** |

These genes have many transcript isoforms, causing massive hit inflation in raw counts.

---

## Impact on Enrichment Analysis

### Functional Gene Set Enrichment

The integrated analysis now uses deduplicated data for all three regions:
- 3'UTR: 13,298 genes with hits
- 5'UTR: 12,314 genes with hits
- Exon: 13,483 genes with hits
- **Total:** 13,979 genes with any TE hits

Top enriched gene sets (using deduplicated data):
1. RNA binding proteins
2. Transcription factors
3. Nervous system development
4. DNA binding proteins

### Shuffled Control Validation

Despite the reduced enrichment ratio (2.16x → 1.91x):
- **High-quality hits remain strongly enriched** (83x for ≥80%/≥50bp)
- **~48% of real hits are genuine TE signal** (not explainable by composition)
- All TE families remain enriched in real sequences

---

## Files Updated

### New Deduplicated Results

| Location | Contents |
|----------|----------|
| `results/3utr_deduplicated/` | Deduplicated 3'UTR analysis |
| `results/5utr_deduplicated/` | Deduplicated 5'UTR analysis |
| `results/exon_analysis/deduplicated/` | Deduplicated exon analysis |
| `results/shuffled_controls/deduplicated/` | Deduplicated shuffled controls |
| `results/integrated_te_analysis/` | Integrated analysis (uses deduplicated data) |

### Key Files

- `*_deduplicated_hits.tsv` - Full deduplicated BLAST output
- `*_deduplication_stats.json` - Per-gene statistics
- `*_original_vs_deduplicated.tsv` - Before/after comparison
- `DEDUPLICATED_COMPARISON.md` - Shuffled control summary

---

## Scripts Updated

| Script | Changes |
|--------|---------|
| `deduplicate_te_hits.py` | Unified script for all region types |
| `analyze_integrated_te_enrichment.py` | Now uses deduplicated UTR data |
| `analyze_shuffled_deduplicated.py` | New script for shuffled control comparison |

---

## Recommendations

1. **Use deduplicated results** for all analyses going forward
2. **Raw results are archived** but should not be used for new analyses
3. **Gene rankings by density** are largely unchanged (density metric is robust)
4. **Total hit counts** should always use deduplicated values
5. **Shuffled control enrichment** is 1.91x (not 2.16x)

---

## Validation Checklist

- [x] Unique hits ≤ raw hits for all genes
- [x] Exon deduplication matches previous results (12.86%)
- [x] Shuffled sequences have ~0% duplication (expected)
- [x] Integrated analysis runs with deduplicated data
- [x] Gene rankings generated from deduplicated data

---

## Conclusion

Deduplication revealed significant hit inflation in UTR analyses:
- **3'UTR was most affected** (47% duplicate hits)
- **Rankings by density are stable** (not significantly changed)
- **Shuffled control enrichment is lower but still strong** (1.91x)
- **High-quality hits remain highly enriched** (83x)

The biological conclusions remain valid, but all quantitative comparisons should use deduplicated data.
