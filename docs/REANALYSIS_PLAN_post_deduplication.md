# Reanalysis Plan: Post-Deduplication

**Created:** 2026-01-15
**Status:** PLANNING

## Summary of Issue

Deduplication revealed that **47% of 3'UTR hits and 34% of 5'UTR hits were duplicates** from overlapping transcript isoforms. All previous analyses using raw hit counts are inflated and need to be rerun with deduplicated data.

| Region | Raw Hits | Unique Hits | Inflation Factor |
|--------|----------|-------------|------------------|
| **Exon** | 1,722,273 | 1,500,719 | 1.15x |
| **5'UTR** | 1,903,660 | 1,255,113 | 1.52x |
| **3'UTR** | 2,573,926 | 1,370,202 | **1.88x** |

**Impact:** Rankings, enrichment statistics, and hit counts in all previous analyses are incorrect.

---

## Affected Analyses

### HIGH PRIORITY - Core Results

| Analysis | Files Affected | Impact |
|----------|---------------|--------|
| **Gene TE rankings** | `top_100_te_genes_FIXED.tsv`, `bottom_100_te_genes_FIXED.tsv` | Rankings may change significantly |
| **5'UTR analysis** | `5utr_analysis/*.tsv` | All metrics inflated by ~1.5x |
| **Functional enrichment** | `functional_te_enrichment.tsv` | Enrichment ratios affected |
| **FlyFISH enrichment** | `flyfish_te_family_enrichment.tsv`, `flyfish_hit_characteristics.tsv` | Localization correlations affected |
| **Integrated analysis** | `integrated_te_analysis/` | Uses raw UTR data |

### MEDIUM PRIORITY - Validation & Comparison

| Analysis | Files Affected | Impact |
|----------|---------------|--------|
| **Shuffled controls** | `shuffled_controls/` | Need dedup on both real AND shuffled |
| **Strand bias** | `strand_bias_by_*.tsv` | Percentages may shift |
| **RepeatMasker comparison** | `repeatmasker_analysis/` | Hit counts affected |

### LOWER PRIORITY - May Not Need Full Redo

| Analysis | Files Affected | Notes |
|----------|---------------|-------|
| **Conservation analysis** | `CONSERVATION_ANALYSIS_RESULTS.md` | Uses genomic coords, not counts |
| **Synteny analysis** | `SYNTENY_ANALYSIS_RESULTS.md` | Uses genomic coords, not counts |
| **UTR TE loci** | `utr_te_loci/` | Already per-isoform, may be OK |

### REPORTS TO UPDATE

All summary reports contain outdated statistics:
- `COMPLETE_TE_ANALYSIS_REPORT.md`
- `GENOME_WIDE_TE_ANALYSIS_CORRECTED.md`
- `5UTR_TE_ANALYSIS_SUMMARY.md`
- `5UTR_PROJECT_SUMMARY.md`

---

## Reanalysis Approach

### Option A: Modify Existing Scripts (Recommended)

Update each analysis script to load deduplicated data instead of raw BLAST results.

**Pros:** Preserves existing code structure, easy to verify changes
**Cons:** Many scripts to update

### Option B: Create Deduplicated BLAST Files

Generate new BLAST-format TSV files with only unique hits, then rerun all scripts unchanged.

**Pros:** Minimal script changes
**Cons:** Need to output full hit details (not just counts)

### Option C: Hybrid

- For count-based analyses: Load from dedup summary files
- For hit-level analyses: Create deduplicated BLAST files

---

## Implementation Plan

### Phase 1: Generate Deduplicated Hit Files

The current `deduplicate_te_hits.py` only outputs statistics. Need to also output:
1. **Full deduplicated hit TSV** - same format as BLAST output, but deduplicated
2. **Gene-level summary with deduplicated counts**

```bash
# Add --output-hits flag to generate full hit files
python scripts/deduplicate_te_hits.py --seq-type 3utr --output-hits
python scripts/deduplicate_te_hits.py --seq-type 5utr --output-hits
```

### Phase 2: Update Core Analysis Scripts

| Script | Changes Needed |
|--------|---------------|
| `analyze_genome_wide_te.py` | Load deduplicated hits |
| `analyze_integrated_te_enrichment.py` | Already partially done for exon; extend to UTR |
| `analyze_functional_enrichment.py` | Load deduplicated counts |
| `rank_gene_sets.py` | Load deduplicated counts |

### Phase 3: Rerun Shuffled Control Analysis

**Critical:** Must deduplicate shuffled sequences too, not just real sequences.

```bash
# Dedup real data
python scripts/deduplicate_te_hits.py --seq-type 3utr --blast-file results/shuffled_controls/real_blast.tsv

# Dedup each shuffled replicate
for i in {1..10}; do
    python scripts/deduplicate_te_hits.py --seq-type 3utr \
        --blast-file results/shuffled_controls/shuffled_rep${i}_blast.tsv
done
```

### Phase 4: Regenerate Reports

After all analyses rerun:
1. Regenerate all summary MD files
2. Update HTML visualizations if needed
3. Archive old results with clear "PRE_DEDUP" label

---

## Validation Checklist

After reanalysis:

- [ ] Total unique hits match deduplication counts
- [ ] No gene has more unique hits than raw hits
- [ ] Top gene rankings changed (expected - validates the fix)
- [ ] Shuffled control enrichment ratios updated
- [ ] All summary reports reflect new numbers

---

## Expected Impact on Findings

### Rankings Will Change
Genes with many isoforms (Dscam1, Sxl, sd) were artificially inflated. After dedup:
- Dscam1 drops from ~6,500 to ~160 hits (97.5% were duplicates)
- Genes with few isoforms will rise in rankings

### Enrichment Ratios May Change
If certain gene sets (e.g., germline genes) have more isoforms than others, their enrichment was artificially inflated. Need to verify if:
- Germ plasm enrichment holds after dedup
- Tissue-specific patterns persist

### Shuffled Control Comparison
The 2.2x enrichment of real vs shuffled may change because:
- Both real AND shuffled need deduplication
- If duplication rates differ, ratio will shift

---

## Files to Archive

Move to `results/archive/pre_dedup_2026-01-15/`:
- All `*_FIXED.tsv` files
- Current summary MD files
- Current enrichment TSV files

---

## Timeline

1. **Phase 1** (Hit file generation): Modify script, generate files
2. **Phase 2** (Core reanalysis): Update scripts, rerun
3. **Phase 3** (Shuffled controls): Apply dedup to all replicates
4. **Phase 4** (Reports): Regenerate all summaries
5. **Validation**: Verify results make sense

---

## Questions to Resolve

1. Should we keep both raw and deduplicated versions of results?
2. Do we need to regenerate HTML visualizations?
3. Should strand bias be calculated from unique hits or weighted by duplication?
4. For shuffled controls, do we expect similar duplication rates as real data?
