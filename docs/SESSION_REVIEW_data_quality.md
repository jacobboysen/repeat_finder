# Session Review: Key Questions, Answers, Data Locations, and Caveats

## Overview

This document summarizes the TE (transposable element) analysis session, identifying what questions were asked, what answers were found, where data lives, and critically, what data quality issues were discovered and which analyses may need to be re-run.

---

## Key Questions Asked & Answers

### Q1: Do germ plasm 3'UTRs have elevated TE content compared to controls?

**Answer:** Yes, but the pattern is more nuanced than expected.

- **piRNA pathway genes (AGO3, aub, vas, piwi)** rank in the top 6-23% genome-wide for TE density
- **Other germ plasm genes (nos, osk, pgc, gcl)** are near average (45-60th percentile)
- **CycB** is an outlier at only 9.6th percentile (below average)
- Housekeeping genes are depleted (33.8th percentile average)

**Interpretation:** The "germ plasm TE enrichment" is driven primarily by piRNA pathway genes, not a general germ plasm phenomenon.

**Data location:** `results/strand_bias_by_utr.tsv`, gene set analysis in session

---

### Q2: What is the genome-wide strand bias pattern?

**Answer:** 60.2% sense, 39.8% antisense (1.51x sense bias)

**Key finding:** Different insertions of the SAME TE family show opposite strand biases:
- roo{}25, roo{}39: 94-95% sense
- roo{}850, roo{}1631: 72-79% antisense

**Interpretation:** Strand bias reflects **insertion orientation** into the genome, not an intrinsic TE property.

**Data location:**
- `results/strand_bias_by_utr.tsv` - per-UTR strand stats
- `results/strand_bias_by_te.tsv` - per-TE strand stats

---

### Q3: What are the top/bottom genes by TE density?

**Answer:**

**Top genes (genome-wide):**
1. FBgn0040959 (Prt-15a) - chitin-binding, density 325,271
2. FBgn0034403 (CG18190) - microtubule-associated
3. FBgn0067905 (Dso2) - antimicrobial peptide
...

**Bottom genes:** Long UTRs, housekeeping/metabolic functions

**Data location:**
- `results/top_100_te_genes_FIXED.tsv` - GOOD (gene-level, properly aggregated)
- `results/bottom_100_te_genes_FIXED.tsv` - GOOD
- `results/top_100_te_transcripts_CORRECTED.tsv` - GOOD (transcript-level)
- `results/bottom_100_te_transcripts_CORRECTED.tsv` - GOOD

**⚠️ BAD DATA:**
- `results/top_100_te_genes.tsv` - OLD, has gene/transcript mixing bug
- `results/bottom_100_te_genes.tsv` - OLD, has gene/transcript mixing bug

---

### Q4: Do TE hits cluster in specific regions of TEs?

**Answer:** Yes, strong 5' bias in LTR retrotransposons:
- 45.3% of hits at 5' end (0-20% of TE)
- 32.7% in middle
- 22.1% at 3' end

LTRs (regulatory regions) are enriched over CDS (coding regions).

**Data location:** `results/INTERMEDIATE_TE_REGION_ENRICHMENT_ANALYSIS.md`, various `INTERMEDIATE_*_te_regions.tsv` files

**⚠️ NOTE:** This analysis was done on small gene sets (germ plasm, housekeeping, somatic, cleared), NOT genome-wide. Should be re-run on genome-wide data for confirmation.

---

### Q5: What % of TE database is covered by UTR hits?

**Answer:**
- 91.6% of TEs (5,255/5,734) have at least one hit
- 70.0% of TE bases covered
- roo family has highest coverage

**Data location:** Analysis in session (not saved to file - SHOULD BE RE-RUN)

---

### Q6: Do top TE-density genes have TE similarity extending beyond the 3'UTR?

**Answer:** YES - critical finding!

Many "high TE density" genes have hits spanning CDS/UTR boundaries:
- **Dso2**: 81% spanning (antimicrobial peptide - likely domesticated TE)
- **CG33093**: 96% spanning (likely TE-derived gene)
- **dpr21**: 84% spanning
- **Sf3b5**: 88% spanning

These may be **domesticated TEs**, not genes with TE insertions.

**Data location:**
- `results/full_transcript_te/index.html` - HTML visualizations
- `results/full_transcript_te/*_full_transcript.html` - per-gene pages

---

## Data Quality Issues & Caveats

### ⚠️ CRITICAL: Gene vs Transcript Mixing Bug (FIXED)

**Problem:** Original analysis aggregated hits across all transcript isoforms per gene, but stored only the last transcript's UTR length.

**Result:** Impossible cases like genes with 36bp UTR but 515 hits.

**Status:** FIXED - re-ran at transcript level and gene level separately.

**Good files:** `*_FIXED.tsv`, `*_CORRECTED.tsv`
**Bad files:** Original `top_100_te_genes.tsv`, `bottom_100_te_genes.tsv` without suffixes

---

### ⚠️ CAUTION: Ranking Confusion

**Problem:** In one analysis, I reported AGO3 as "#1" when it was only #1 among our 47 gene set genes, not genome-wide.

**Actual ranks:**
- AGO3: #787 genome-wide (94.2nd percentile) - NOT #1
- Prt-15a: #1 genome-wide (density 325,271)

**Status:** Clarified in session, but verify any "rank" claims carefully.

---

### ⚠️ CAUTION: Strand Annotation in Early Analyses

**Problem:** Early in the project (before this session), strand information may not have been consistently annotated. The `sstart > send` convention for antisense was discovered during troubleshooting.

**Status:** Current analyses use strand correctly, but older result files may not have strand columns or may have them mislabeled.

**Verify:** Any file without explicit `strand` or `sense_pct` columns should be treated as potentially incomplete.

---

### ⚠️ CAUTION: "TE-derived genes" vs "genes with TE insertions"

**Problem:** The top TE-density genes may not be what they appear.

**Finding:** Dso2, CG33093, dpr21, etc. have TE similarity spanning their CDS, suggesting they ARE domesticated TEs, not genes that happen to have TE insertions.

**Implication:** Analyses that assume "high TE density = more TE insertions" may be confounded by TE-derived genes.

**Recommendation:** For future analysis, separate:
1. Genes with UTR-confined TE hits (true insertions)
2. Genes with CDS-spanning TE hits (possible domesticated TEs)

---

### ⚠️ INFO: DUST Filtering

**Context:** Original parameter sweep (before this session) found 90% simple repeats. DUST filtering (`dust=yes`) was added to fix this.

**Status:** All current analyses use `dust=yes`. Files in `results/dust_sweep/` are the corrected versions.

**Bad data:** Any BLAST results without DUST filtering may be dominated by AT-repeats.

---

## File Inventory

### GOOD DATA (verified/corrected)

| File | Description | Status |
|------|-------------|--------|
| `results/top_100_te_genes_FIXED.tsv` | Gene-level top 100, properly aggregated | ✓ Good |
| `results/bottom_100_te_genes_FIXED.tsv` | Gene-level bottom 100 | ✓ Good |
| `results/top_100_te_transcripts_CORRECTED.tsv` | Transcript-level top 100 | ✓ Good |
| `results/bottom_100_te_transcripts_CORRECTED.tsv` | Transcript-level bottom 100 | ✓ Good |
| `results/strand_bias_by_utr.tsv` | Per-UTR strand statistics (genome-wide) | ✓ Good |
| `results/strand_bias_by_te.tsv` | Per-TE strand statistics (genome-wide) | ✓ Good |
| `results/genome_wide_all_3utrs.tsv` | Raw BLAST results (2.57M hits) | ✓ Good |
| `results/full_transcript_te/*.html` | Full transcript visualizations | ✓ Good |
| `results/te_annotations/*.html` | Germ plasm gene visualizations | ✓ Good |
| `results/te_annotations_genomewide/*.html` | Top/bottom 10 gene visualizations | ✓ Good |
| `results/GENOME_WIDE_TE_ANALYSIS_CORRECTED.md` | Final genome-wide summary | ✓ Good |

### BAD DATA (labeled with BAD_ prefix)

These files contain **incorrect data** due to bugs or lack of proper filtering.

| File/Directory | Issue |
|----------------|-------|
| `BAD_gene_transcript_mixing_top_100_te_genes.tsv` | Gene/transcript mixing bug |
| `BAD_gene_transcript_mixing_bottom_100_te_genes.tsv` | Gene/transcript mixing bug |
| `BAD_pre_DUST_parameter_sweep/` | Pre-DUST, 90% simple repeats |
| `BAD_pre_DUST_controls/` | Pre-DUST filtering |
| `BAD_pre_DUST_density_analysis/` | Pre-DUST filtering |
| `BAD_pre_DUST_te_family_analysis/` | Pre-DUST filtering |
| `BAD_pre_DUST_antisense_best.tsv` | Pre-DUST filtering |
| `BAD_pre_DUST_shuffled_best.tsv` | Pre-DUST filtering |
| `BAD_wrong_rankings_GENOME_WIDE_COMPARISON.md` | Wrong gene rankings (pre-bug fix) |
| `BAD_references_buggy_data_TOP_BOTTOM_100_ANALYSIS.md` | References buggy top/bottom 100 files |

### INTERMEDIATE DATA (labeled with INTERMEDIATE_ prefix)

These files are from **early exploratory analysis** that was superceded by genome-wide analysis. Not necessarily wrong, but incomplete or superseded.

| File/Directory | Description |
|----------------|-------------|
| `INTERMEDIATE_diverged_controls/` | Small gene set BLAST results (pre-genome-wide) |
| `INTERMEDIATE_dust_sweep/` | Parameter testing data |
| `INTERMEDIATE_DIVERGED_TE_ANALYSIS_SUMMARY.md` | Early DUST parameter analysis summary |
| `INTERMEDIATE_DIVERGED_TE_CONTROL_COMPARISON.md` | Early control comparison |
| `INTERMEDIATE_small_geneset_STRAND_ANALYSIS_SUMMARY.md` | Strand analysis on small gene sets (not genome-wide) |
| `INTERMEDIATE_small_geneset_germ_plasm_strand_summary.tsv` | Small gene set strand data |
| `INTERMEDIATE_TE_REGION_ENRICHMENT_ANALYSIS.md` | TE region analysis (references intermediate data) |
| `INTERMEDIATE_*_te_regions.tsv` | Small gene set TE region mapping |
| `INTERMEDIATE_*_te_regions.blast.tsv` | Small gene set BLAST results |
| `INTERMEDIATE_germ_plasm_dust_on.tsv` | Small gene set DUST-corrected data |
| `INTERMEDIATE_shuffled_dust_optimized.tsv` | Shuffled control data |
| `INTERMEDIATE_Kr_3UTR_TE_annotation*` | Early Krüppel gene visualizations |
| `INTERMEDIATE_te_fossil_candidates.txt` | Early TE fossil list |
| `INTERMEDIATE_te_fossils_diverged.txt` | Early diverged TE list |
| `INTERMEDIATE_top_hits_*.txt` | Early exploratory top hits |
| `INTERMEDIATE_test_data_genome_wide_sample_500.tsv` | Test sample data |

### DOCUMENTATION

| File | Description | Status |
|------|-------------|--------|
| `docs/SESSION_SUMMARY_diverged_TE_analysis.md` | Main session summary | ✓ Current |
| `docs/SESSION_REVIEW_data_quality.md` | This file - data quality audit | ✓ Current |
| `results/GENOME_WIDE_TE_ANALYSIS_CORRECTED.md` | Final genome-wide analysis | ✓ Good |

### SCRIPTS

| Script | Description |
|--------|-------------|
| `scripts/analyze_genome_wide_te.py` | Gene-level aggregation, strand bias |
| `scripts/analyze_full_transcripts.py` | Full transcript (5'UTR+CDS+3'UTR) analysis |
| `scripts/compare_gene_sets.py` | Compare gene sets to genome-wide |
| `scripts/analyze_te_regions.py` | Map hits to TE structural regions |

---

## Recommended Re-Analysis

Before extending analysis, confirm these with fresh runs:

1. **Gene-level density rankings** - Re-run `analyze_genome_wide_te.py` and verify top 100 matches `top_100_te_genes_FIXED.tsv`

2. **TE database coverage stats** - These were calculated in-session but not saved to a file. Re-run and save.

3. **Strand bias patterns** - Verify the 60.2% sense finding with fresh analysis

4. **Full transcript analysis** - The "domesticated TE" finding for Dso2/CG33093 is important - verify independently

---

## Key Takeaways for Future Analysis

1. **Always work at transcript level first**, then aggregate to gene level explicitly

2. **Always include strand information** (check `sstart < send` for sense)

3. **Separate "UTR-only" hits from "spanning" hits** - these represent different biology

4. **DUST filtering is essential** - without it, results are dominated by simple repeats

5. **piRNA pathway genes are outliers** - they drive much of the "germ plasm enrichment" signal

6. **Some top genes are domesticated TEs** - don't assume high TE density = more insertions

---

*Document created: 2026-01-14*
