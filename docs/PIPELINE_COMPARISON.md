# Pipeline Comparison: Legacy vs fossil_finder (2026-03-11)

## Purpose

Documents the expected discrepancy between the legacy `scripts/blast_runner.py` results
(Jan 2025) and the new `fossil_finder` package results (Mar 2026). When re-running
via Snakemake, use this to confirm that differences are parameter-driven, not bugs.

## Parameter Differences

| Parameter | Legacy (`scripts/blast_runner.py`) | New (`fossil_finder`) | Impact |
|-----------|------------------------------------|-----------------------|--------|
| word_size | **11** | **7** | Shorter seeds find diverged TE fossils |
| evalue | **1e-5** | **10** | 10x more permissive; captures weak aggregate signal |
| dust | **yes** | **no** | DUST masking discarded real TE signal in simple-repeat-adjacent regions |
| max_target_seqs | 500 | 1000 | More target TEs per query |
| outfmt columns | 14 (no qseq/sseq) | 16 (with qseq/sseq) | New pipeline retains alignment sequences |
| Query IDs | FBtr (transcript) | Coordinate-based | Different gene-mapping approach |
| Dedup strategy | Post-BLAST hit dedup | Pre-BLAST sequence dedup | New: fewer queries, no isoform redundancy |

The parameter revision (word_size=7, dust=no) is documented in `docs/DUST_FILTERING_ANALYSIS.md`
(2026-01-22). The legacy results predate that revision.

## Quantitative Comparison (3'UTR, dm6)

| Metric | Legacy (e=1e-5, ws=11, dust=yes) | New (e=10, ws=7, dust=no) |
|--------|----------------------------------|---------------------------|
| Query regions | ~28K transcripts | 21,661 deduplicated sequences |
| Raw BLAST hits | 2,573,926 | 7,539,059 |
| After dedup | 1,370,202 (hit dedup) | N/A (sequences pre-deduped) |
| Genes with hits | 13,298 | 13,156 |
| TE families | ~5,600 | 5,612 |

## Rank Concordance (Top Genes by TE Density)

All 100 legacy top genes appear in the new results. Rankings shift because the
sensitive parameters reveal previously undetectable diverged fossils.

| Legacy top N | Found in new top ... | Interpretation |
|-------------|----------------------|----------------|
| Top 100 | 27 in top 100 | Rankings reshuffle with more signal |
| Top 100 | 71 in top 500 | 71% remain in top 3.7% |
| Top 100 | 87 in top 1000 | 87% remain in top 7.5% |
| Top 100 | 100 in top 13,375 | None lost entirely |

Legacy #6 (FBgn0037514) → New #1. Legacy #1 (FBgn0040959, 48bp UTR) → New #47.
Short UTRs are disproportionately affected by evalue relaxation (more noise per bp).

## Expected Behavior on Snakemake Re-run

When reproducing with Snakemake using `dmel_r6.66.yaml` blast config
(word_size=7, evalue=0.001 default or evalue=10 for exhaustive):

1. Hit counts should match `results/dmel_3utr_e10/` at same evalue
2. Gene rankings should match within stochastic BLAST variation
3. Legacy results in `results/genome_wide_all_3utrs.tsv` are NOT expected to match
   and should be treated as superseded

## Notes

- `max_target_seqs` saturates for 308/21,661 queries at evalue=10 — BLAST may
  not return the globally best 1000 targets. This is a known BLAST behavior
  (Shah et al. 2019). For reproducibility, consider raising to 5000 or using
  `-subject` mode for small databases.
- The legacy `blast_results.tsv` and new `blast_results.tsv` are not directly
  comparable even at the same evalue threshold, because different word_size
  produces different seed-extension alignments.
