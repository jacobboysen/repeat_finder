# TE Fossil Analysis Architecture

> Comprehensive analysis suite for Drosophila melanogaster TE fossil characterization.
> Designed for the `fossil_finder` pipeline output (evalue=10, word_size=7, dust=no).

**Date**: 2026-03-11
**Pipeline data**: `results/dmel_3utr_e10/`, `results/dmel_5utr_e10/`
**Applies to both UTR types** — same structure, same workstreams, compared in report.

---

## Dependency Graph

```
Pipeline outputs (blast_results.tsv, regions.tsv, gene_stats.tsv, regions.fa)
|
+--> WS1: Match Distributions & Patterns
|    +- PRODUCES: hit_quality_tiers.tsv, alignment_stats.tsv,
|    |            positional_profiles_*.tsv, te_coverage_maps.tsv
|    |
|    +-----------> WS3: TE Family Enrichment (uses quality tiers)
|    +-----------> WS5: RepeatMasker Overlap (uses quality tiers)
|    |               +- PRODUCES: known_hits.tsv, novel_hits.tsv
|    |               |
|    |               +--> WS6: Conservation & Synteny (uses known/novel split)
|    |                     +- PRODUCES: phylop_by_hit.tsv, synteny_results.tsv
|    |                     |
|    |                     +--> WS2: Motif Analysis (uses conservation scores)
|    |
|    +-----------> WS4: Unusual Matches (re-examines WS1 outliers + WS5 novel hits)
|
+--> WS7: GO & Functional Enrichment (independent, reads gene_stats.tsv + annotations)
|
+--> Summary Report (reads all workstream outputs)
```

**Execution order**: WS1 -> (WS3, WS5, WS7 in parallel) -> WS6 -> WS4 -> WS2 -> Report

### Inter-Module Data Flows

| Upstream | Downstream | Data passed | Purpose |
|----------|-----------|-------------|---------|
| WS1 | WS3, WS5, WS6 | `hit_quality_tiers.tsv` | Quality stratification changes every downstream analysis |
| WS5 | WS6, WS4 | `known_hits.tsv`, `novel_hits.tsv` | Novel hits are the interesting ones for conservation |
| WS6 | WS2 | `phylop_by_hit.tsv` | Motif significance depends on whether the region is conserved |
| WS1 | WS4 | `alignment_stats.tsv` | Unusual matches are outliers from the WS1 distribution |
| WS3 | Report | `family_enrichment_by_set.tsv` | Top families per gene set |
| WS1 | Report | All positional/distribution files | Spread-level detail for exploration |

---

## WS1: Match Distributions & Patterns

**Purpose**: Characterize the breadth and shape of TE<>target alignments. Foundation
workstream -- every downstream analysis uses its quality tiers and per-hit statistics.

**Why this is first**: The evalue=10 search produces 7.5M 3'UTR hits spanning a huge
quality range. Without stratification, aggregate statistics are dominated by millions
of marginal hits. WS1 establishes the quality tiers that make every other analysis
interpretable.

### Inputs

| File | Source | Key columns |
|------|--------|-------------|
| `blast_results.tsv` | Pipeline | qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore, qlen, slen, **qseq, sseq** |
| `regions.tsv` | Pipeline | region_id, chrom, start, end, strand, parent_id, length |

### Analyses

#### 1.1 Alignment Quality Distribution
- pident x length 2D histogram (heatmap)
- evalue distribution (log-scale, cumulative)
- Bitscore distribution
- **Quality tier assignment**:
  - `strict`: pident >= 85 AND length >= 100
  - `moderate`: pident >= 75 AND length >= 50
  - `relaxed`: everything else
- Log tier counts and percentages

#### 1.2 Edit Distance & Mismatch Profile
- From qseq/sseq: per-position mismatch rate, gap frequency
- Transition/transversion ratio (diagnostic: neutral decay has ti/tv ~2, selection skews it)
- Gap length distribution
- **Per-hit stats**: mismatch_count, gap_count, ti_count, tv_count, edit_distance, mismatch_rate
- Aggregate: mean/median/std/quartiles per quality tier

#### 1.3 Positional Bias on Target (UTR)
- Normalize hit midpoint by UTR length -> decile bins (0-10%, 10-20%, ..., 90-100%)
- Compare 3'UTR vs 5'UTR positional profiles
- Stratify by quality tier
- Legacy found 2.27x enrichment at 3' end of 3'UTRs -- validate with new params
- **Per-hit stat**: normalized_position (0.0 = UTR start, 1.0 = UTR end)

#### 1.4 Positional Pileup on TE
- For each TE (sseqid), normalize hit positions (sstart/send) by TE length (slen)
- Decile-binned coverage per TE
- Aggregate across all TEs and per top-20 most-hit TEs
- Reveals whether hits cluster at LTR, gag, pol, or other domains
- **Per-hit stat**: te_normalized_start, te_normalized_end

#### 1.5 Per-TE Unique Hit Counts & Gene Breadth
- Per TE: n_unique_queries, n_unique_genes, mean_pident, total_bp
- Gene entropy: Shannon entropy of gene distribution (high = broadly distributed, low = concentrated)
- Rank TEs by: total hits, unique genes, gene entropy
- **Individual-level**: every TE gets a row with all metrics

#### 1.6 Hit Multiplicity
- Hits per query: distribution (how many TEs match each UTR region?)
- Hits per gene: distribution
- Queries per TE: distribution
- Genes per TE: distribution
- Report: mean, median, max, quartiles for each

### Outputs

| File | Granularity | Key columns |
|------|-------------|-------------|
| `hit_quality_tiers.tsv` | Per hit (7.5M rows for 3'UTR) | All BLAST columns + tier, mismatch_rate, gap_rate, edit_distance, normalized_utr_pos, te_normalized_start, te_normalized_end |
| `alignment_stats_summary.tsv` | Per quality tier (3 rows) | tier, n_hits, mean_pident, mean_length, mean_evalue, mean_mismatch_rate, mean_gap_rate, mean_ti_tv, median_*, q25_*, q75_* |
| `positional_profiles_utr.tsv` | Per decile x tier x utr_type (60 rows) | decile, tier, utr_type, hit_count, hit_density, normalized_density |
| `positional_profiles_te.tsv` | Per TE x decile (top 100 TEs, 1000 rows) | te_id, decile, hit_count, coverage_fraction |
| `te_hit_breadth.tsv` | Per TE (5612 rows for 3'UTR) | te_id, n_hits, n_unique_queries, n_unique_genes, gene_entropy, mean_pident, mean_length, total_bp, top_gene_1, top_gene_2, top_gene_3 |
| `hit_multiplicity.tsv` | Summary (4 distributions) | metric, mean, median, q25, q75, max, n |

### Code Approach

- **Reuse**: Legacy `analyze_utr_position_bias.py` for UTR positional normalization logic
- **New**: qseq/sseq alignment parsing (legacy outfmt lacked these columns)
- **Module**: `src/fossil_finder/analysis/match_patterns.py` for reusable functions; analysis script calls these

### Caveats

- Edit distance from qseq/sseq is *aligned* edit distance, not global. Unaligned flanks not captured.
- For short alignments (<30bp), mismatch statistics are noisy -- report N per bin.
- `hit_quality_tiers.tsv` will be ~2GB for 3'UTR. Consider parquet for future iterations.
- Tier thresholds (85/100, 75/50) are heuristic. Document as configurable constants.

---

## WS3: TE Family Enrichment

**Purpose**: Which TE families drive the signal, and do they differ between gene sets?

### Inputs

| File | Source |
|------|--------|
| `blast_results.tsv` | Pipeline |
| `gene_stats.tsv` | Pipeline |
| `hit_quality_tiers.tsv` | WS1 |
| `data/gene_lists/*_fbgn_ids.txt` | Gene sets (5 lists) |

### Analyses

#### 3.1 Global Family Ranking
- All families ranked by: hit_count, total_bp, n_unique_genes, gene_entropy
- Stratified by quality tier (strict/moderate/relaxed)
- **Per-family stats at each tier**: n_hits, total_bp, mean_pident, mean_evalue, n_genes, gene_entropy, sense_hits, antisense_hits

#### 3.2 Family Enrichment by Gene Set
- For each gene set (germ_plasm, housekeeping, somatic, cleared, adult):
  - Compute per-family hit frequency in set vs background
  - Fold enrichment = freq_in_set / freq_in_background
  - Log2 fold enrichment for symmetric comparison
- Uses existing `compute_fold_enrichment()` from `fossil_finder.analysis.families`
- **Per-family x per-set**: fold_enrichment, log2_enrichment, count_in_set, count_in_bg

#### 3.3 TE Class Distribution
- Map TE IDs to class (LTR, LINE, DNA transposon, rolling circle, etc.)
- Class-level hit counts, gene counts, bp totals
- Stratify by quality tier
- **Requires**: TE ID -> family -> class lookup table (see Data Validation)

#### 3.4 Family x Strand Interaction
- Per-family sense/antisense ratio
- Flag families with significant strand bias (binomial test, BH-corrected)
- Uses existing `compute_te_strand_bias()` from `fossil_finder.analysis.strand`
- **Per-family**: sense_hits, antisense_hits, sense_pct, p_value_binomial, fdr

### Outputs

| File | Granularity |
|------|-------------|
| `family_ranking_global.tsv` | Per TE x tier (~16K rows) |
| `family_enrichment_by_set.tsv` | Per TE x gene set (~28K rows) |
| `class_distribution.tsv` | Per class x tier (~30 rows) |
| `family_strand_bias.tsv` | Per TE (~5.6K rows) |
| `te_id_to_family_class.tsv` | Lookup table (all TEs) |

### Code Approach

- **Reuse**: `fossil_finder.analysis.families` module (compute_family_stats, compute_fold_enrichment, compute_class_distribution)
- **Reuse**: Legacy `analyze_te_families.py` for gene-set stratification logic
- **New**: TE ID -> family -> class mapping builder (parse FBti metadata + consensus names)

### Data Validation

- **TE class metadata**: Legacy used `te_annotations.gff` for structural domains. Need to verify coverage of all 5,612 families in new results. FlyBase TE instance IDs (FBti*) need lookup to family name and class. Consensus sequence names in `dmel_te_consensus.fasta` contain family info directly.
- **Action**: Build `te_id_to_family_class.tsv` from FlyBase TE metadata + consensus FASTA headers. Log coverage (how many of 5,612 sseqids can be classified).

### Caveats

- sseqid is a mix of FBti instance IDs and consensus names depending on which DB was used. Need unified mapping.
- Some TE families will have very few hits -- fold enrichment is noisy for small N. Report N alongside enrichment.
- Class distribution depends on mapping quality. Report "unclassified" count.

---

## WS5: RepeatMasker Overlap

**Purpose**: What fraction of BLAST hits overlap RepeatMasker annotations? Novel hits
(not in RM) are the most interesting candidates for uncharacterized TE fossils.

### Inputs

| File | Source |
|------|--------|
| `blast_results.tsv` | Pipeline |
| `regions.tsv` | Pipeline |
| `hit_quality_tiers.tsv` | WS1 |
| `data/references/dm6.fa.out` | UCSC RepeatMasker |

### Analyses

#### 5.1 Overlap Classification
- Convert each BLAST hit to genomic coordinates (region chrom:start + qstart/qend)
- Check overlap with RepeatMasker intervals (same chrom, overlapping coords)
- Uses existing `fossil_finder.analysis.repeatmasker` module
- **Per-hit**: is_known (bool), rm_name, rm_class, rm_family (if overlapping)

#### 5.2 Known vs Novel by Quality Tier
- Cross-tabulate: quality tier x known/novel
- Are strict-quality hits more or less likely to be in RepeatMasker?
- Legacy found 83.2% novel overall -- re-measure with new params
- **Per-tier**: n_known, n_novel, pct_novel

#### 5.3 Novel Hit Characterization
- Novel hits by: pident distribution, length distribution, TE family distribution
- Compare novel vs known on these axes
- **Per-hit in novel set**: all BLAST columns + tier + alignment stats from WS1

### Outputs

| File | Granularity |
|------|-------------|
| `known_hits.tsv` | Per known hit (subset of all hits) |
| `novel_hits.tsv` | Per novel hit (subset of all hits) |
| `rm_overlap_by_tier.tsv` | Per tier (3 rows) |
| `novel_vs_known_stats.tsv` | Summary statistics (pident, length distributions) |
| `rm_stats.json` | Top-level counts |

### Code Approach

- **Reuse**: `fossil_finder.analysis.repeatmasker` (parse_repeatmasker_out, find_overlaps, classify_hits)
- **Adapt**: Genomic coordinate conversion for coordinate-based query IDs (legacy used FBtr with separate coordinate lookup)

### Data Validation

- `dm6.fa.out`: UCSC RepeatMasker output, 137,555 repeat annotations. Verified present. 1-based coordinates matching GFF convention.

### Caveats

- RepeatMasker uses different sensitivity than BLAST (tuned for recent/intact TEs). High novel rate is expected, not indicative of false positives.
- Overlap is positional, not sequence-based. A BLAST hit at the same coordinates as an RM annotation doesn't mean they found the same TE -- could be overlapping but different elements.
- Coordinate conversion from query-relative (qstart/qend) to genomic requires strand awareness.

---

## WS6: Conservation, Synteny & PhyloP (Quality Paradox)

**Purpose**: Reproduce and extend the "quality paradox" -- do more degraded TE fossils
show MORE functional signal? Map conservation across Drosophilids and wider insects.

### Inputs

| File | Source |
|------|--------|
| `blast_results.tsv` | Pipeline |
| `hit_quality_tiers.tsv` | WS1 |
| `known_hits.tsv`, `novel_hits.tsv` | WS5 |
| `data/references/dm6.phyloP27way.bw` | UCSC (396 MB) |
| `data/references/maf/*.maf.gz` | UCSC 27-way MAF (6 chroms) |
| `data/references/tools/bigWigAverageOverBed` | UCSC utility |

### Analyses

#### 6.1 PhyloP Scores by Quality Tier
- Convert each hit to genomic BED coordinates
- Extract phyloP scores using bigWigAverageOverBed
- Mean/median conservation per quality tier
- Compare TE-hit positions vs flanking non-TE positions (same UTR)
- **Per-hit**: mean_phylop, median_phylop, min_phylop, max_phylop

#### 6.2 Synteny by Quality Tier
- For each hit's genomic interval, check synteny in MAF alignments
- Species: droSim1 (D. simulans), droYak2 (D. yakuba), droEre2 (D. erecta), droSec1 (D. sechellia)
- Fraction syntenic per tier
- **Per-hit x per-species**: is_syntenic (bool), aligned_pident (if syntenic)

#### 6.3 Known (RM) vs Novel Conservation
- Compare phyloP distributions: known hits vs novel hits
- Do novel hits (not in RepeatMasker) have different conservation signatures?
- **Summary**: mean/median phyloP for known vs novel, Mann-Whitney test

#### 6.4 Per-Species Conservation Heatmap
- For each Drosophilid + wider insect species in the 27-way alignment:
  - What fraction of TE hits are syntenic?
  - Broken down by quality tier
- Reveals evolutionary depth of TE fossil retention
- **Per-species x per-tier**: n_syntenic, n_total, fraction_syntenic

#### 6.5 Pident vs Conservation Correlation (Core Quality Paradox)
- Scatter: BLAST pident (x) vs mean phyloP (y) for each hit
- Spearman correlation per tier and overall
- The paradox: lower pident (more diverged TE fossil) may correlate with HIGHER conservation (under functional constraint)
- **Per-hit**: pident, mean_phylop (for scatter data)
- **Summary**: rho, p-value, n per tier

### Outputs

| File | Granularity |
|------|-------------|
| `phylop_by_hit.tsv` | Per hit (all hits with phyloP scores) |
| `phylop_by_quality.tsv` | Per tier (3 rows, summary stats) |
| `synteny_by_hit.tsv` | Per hit x species |
| `synteny_by_quality.tsv` | Per tier x species |
| `synteny_by_species.tsv` | Per species (4+ rows) |
| `known_vs_novel_conservation.tsv` | Summary (2 rows + test stats) |
| `pident_vs_conservation.tsv` | Per hit (paired values for plotting) |
| `quality_paradox_stats.json` | Correlation coefficients, p-values |

### Code Approach

- **Reuse**: Legacy `conservation_analysis/convert_all_hits_to_bed.py` (BED conversion), `analyze_all_hits_conservation.py` (phyloP extraction), `fast_synteny_analysis.py` (MAF parsing, O(n+m) merge-sort)
- **Adapt**: Coordinate-based query IDs (legacy assumed FBtr format for coordinate lookup). New pipeline has chrom:start-end in query ID, simplifying genomic coordinate derivation.

### Data Validation

- `dm6.phyloP27way.bw`: UCSC source, 396 MB. Verified present. Covers all 6 major chromosomes.
- `maf/*.maf.gz`: 6 files for 2L, 2R, 3L, 3R, X, 4. Verified present.
- `bigWigAverageOverBed`: UCSC binary in `data/references/tools/`. Verified present.
- **Species in 27-way**: D. simulans, D. yakuba, D. erecta, D. sechellia (close); D. persimilis, D. pseudoobscura, D. willistoni (medium); wider Diptera (distant). The 4 close Drosophilids are most informative for TE fossil dating.

### Caveats

- phyloP 27-way includes non-insect outgroups. Scores reflect deep conservation for universally conserved genes. For Drosophila-specific analysis, per-species synteny from MAF is more appropriate. Report both but interpret carefully.
- PhyloP extraction is I/O-bound (bigWig random access on 7.5M intervals). Consider sampling for initial exploration, full run for final analysis.
- Synteny = positional conservation, not sequence conservation. A region can be syntenic but highly diverged, or non-syntenic due to rearrangement (not loss).

---

## WS7: GO & Functional Enrichment

**Purpose**: What biological functions are enriched among TE-dense genes? Independent
workstream -- does not depend on WS1 quality tiers (operates on gene-level stats).

### Inputs

| File | Source |
|------|--------|
| `gene_stats.tsv` | Pipeline |
| `data/annotations/raw/gene_association.gaf` | FlyBase GO |
| `data/annotations/raw/gene_rpkm_report.tsv` | FlyBase RNA-seq |
| `data/annotations/raw/flyfish_localization.csv` | FlyFISH |
| `data/annotations/gene_annotations.tsv` | Unified annotations |
| `data/references/fbgn_to_symbol.tsv` | Gene symbol mapping |

### Analyses

#### 7.1 GO Term Enrichment
- Define TE-dense genes: top 10% by density (not presence -- 98.4% have hits at evalue=10)
- For each GO term with >= 5 genes in the annotation set:
  - Fisher's exact test (TE-dense vs not)
  - Benjamini-Hochberg FDR correction
- Report top terms by: Biological Process, Molecular Function, Cellular Component
- Also test top 5% and top 25% thresholds for sensitivity analysis
- **Per-GO-term**: go_id, go_term, go_domain, n_in_term, n_te_dense_in_term, odds_ratio, p_value, fdr

#### 7.2 Expression Correlation
- Per-tissue RPKM vs TE density (Spearman correlation)
- Tissues: ovary, testis, larval CNS, adult head, whole body, etc.
- Is TE density associated with maternal/germline expression?
- **Per-tissue**: tissue, rho, p_value, n_genes, mean_rpkm_te_dense, mean_rpkm_background

#### 7.3 Tissue-Specific Enrichment
- For each tissue: define tissue-enriched genes (RPKM > 2x median across tissues)
- Fisher's exact test: tissue-enriched x TE-dense
- **Per-tissue**: tissue, n_enriched, n_te_dense_in_enriched, odds_ratio, p_value, fdr

#### 7.4 FlyFISH Localization Enrichment
- RNA subcellular localization patterns vs TE density
- Patterns: pole cell, posterior, anterior, apical, ubiquitous, etc.
- Fisher's exact test per pattern
- **Per-pattern**: pattern, n_genes, n_te_dense, odds_ratio, p_value, fdr

### Outputs

| File | Granularity |
|------|-------------|
| `go_enrichment.tsv` | Per GO term (~thousands of rows) |
| `go_enrichment_top50.tsv` | Top 50 terms per domain (150 rows, for quick view) |
| `expression_correlation.tsv` | Per tissue (~30 rows) |
| `tissue_enrichment.tsv` | Per tissue (~30 rows) |
| `flyfish_enrichment.tsv` | Per localization pattern (~20 rows) |
| `te_dense_genes.tsv` | Per gene in top 10% by density (~1.3K rows, with symbol/GO/expression) |

### Code Approach

- **Reuse**: Legacy `analyze_functional_enrichment.py` (Fisher/MW logic), `utils/annotation_loaders.py` (GO parser, RPKM loader, FlyFISH loader)
- **Reuse**: `fossil_finder.analysis.enrichment` for statistical tests + FDR correction
- **Adapt**: Legacy filtered on TE presence (binary). Need to adapt for density percentile thresholds since 98.4% of genes have hits at evalue=10.

### Data Validation

| Source | File | Check |
|--------|------|-------|
| GO annotations | `gene_association.gaf` | Verify FlyBase release >= r6.60. Check gene count vs our 13K. |
| RNA-seq RPKM | `gene_rpkm_report.tsv` | Verify present, check tissue count. |
| FlyFISH | `flyfish_localization.csv` | Static dataset (completed experiment). Verify present. |
| Gene symbols | `fbgn_to_symbol.tsv` | 8,580 genes -- may miss some of 13K. Log coverage. |

### Caveats

- **Critical**: At evalue=10, 98.4% of genes are TE-positive. Binary presence/absence enrichment is nearly meaningless. Use density percentile (top 10%) as the threshold, not presence.
- GO enrichment is sensitive to the background set. Use all 13K analyzed genes as background, not the full genome.
- Multiple testing correction across thousands of GO terms requires strict FDR. Report both raw p and FDR q.
- Gene symbol mapping covers ~8.5K of ~13K genes. Report unmapped count.

---

## WS4: Unusual Matches & Patterns

**Purpose**: Flag outliers and interesting match patterns for manual investigation.
This workstream synthesizes WS1 distributions and WS5 novel/known classification
to identify the most biologically interesting hits.

### Inputs

| File | Source |
|------|--------|
| `hit_quality_tiers.tsv` (with alignment_stats) | WS1 |
| `novel_hits.tsv` | WS5 |
| `te_hit_breadth.tsv` | WS1 |

### Analyses

#### 4.1 Outlier Hits
- Unusually long alignments (> 99th percentile length for that TE family)
- Unusually high pident in otherwise diverged families (> 95% pident where family mean < 75%)
- Very low evalue hits at otherwise marginal TEs
- **Per-hit**: outlier_type (long_alignment | high_pident | low_evalue), z_score, context

#### 4.2 Multi-Family Hotspots
- Genes hit by many different TE families
- Threshold: >= 50 distinct TE families (top 1% of gene distribution)
- Could indicate insertion hotspots, chimeric fossils, or simple-repeat-adjacent regions
- **Per-gene**: gene_id, n_te_families, top_families, total_hits, total_bp, density

#### 4.3 Clustered Hits on UTR
- Hits that pile up within 50bp of each other on a single UTR
- Merge overlapping/adjacent hits into clusters
- Report cluster size, span, number of distinct TEs in cluster
- **Per-cluster**: region_id, cluster_start, cluster_end, cluster_span, n_hits, n_distinct_tes, mean_pident

#### 4.4 Antisense-Only Genes
- Genes where ALL TE hits are antisense (strand = "minus")
- At minimum 5 hits to avoid noise
- Potential regulatory insertions (antisense TE in mRNA could form dsRNA)
- **Per-gene**: gene_id, n_antisense_hits, n_total_hits, top_te_families

#### 4.5 Novel High-Quality Candidates
- Strict-tier hits NOT in RepeatMasker
- These are the top exaptation candidates: high-quality TE match that RepeatMasker missed
- **Per-hit**: all BLAST columns + tier + rm_status + phylop (if available from WS6)

### Outputs

| File | Granularity |
|------|-------------|
| `outlier_hits.tsv` | Per outlier hit (variable, ~hundreds) |
| `multifamily_hotspots.tsv` | Per gene (top 1%, ~130 genes) |
| `clustered_hits.tsv` | Per cluster (variable) |
| `antisense_only_genes.tsv` | Per gene (variable, likely ~100-500) |
| `novel_hq_candidates.tsv` | Per hit (strict tier + novel, likely ~thousands) |

### Code Approach

- **Reuse**: Legacy `analyze_full_transcripts.py` for hit clustering logic
- **New**: Outlier detection, multi-family hotspot analysis, antisense-only filtering

### Caveats

- "Unusual" is relative to WS1 distributions. Thresholds are exploratory and should be documented as tunable constants.
- Multi-family hotspots may be enriched for long UTRs (more room for hits). Normalize by UTR length.
- Antisense-only requires sufficient hit count to be meaningful. Genes with 1-2 hits that happen to be antisense are not interesting.

---

## WS2: Motif Analysis

**Purpose**: What sequence motifs are enriched in TE-containing UTR regions?
Last in execution order because it depends on WS1 (quality tiers), WS6 (conservation),
and requires generating shuffled controls (compute-intensive).

### Inputs

| File | Source |
|------|--------|
| `regions.fa` | Pipeline |
| `blast_results.tsv` | Pipeline |
| `hit_quality_tiers.tsv` | WS1 |
| `phylop_by_hit.tsv` | WS6 (optional, enriches motif interpretation) |

### Analyses

#### 2.1 Generate Shuffled Controls
- Dinucleotide shuffle of each query UTR, 10 replicates
- BLAST each replicate against same TE database
- **Compute note**: 10 replicates x 21K queries x evalue=10 = ~16 hours per UTR type. Consider running at evalue=0.001 (sufficient for enrichment comparison) or sampling 5K queries. Flag as compute decision.

#### 2.2 K-mer Enrichment (6-mers)
- Count 6-mer frequency in TE-hit regions vs non-hit regions of same UTRs
- Compare to shuffled control frequencies
- Z-score = (observed - mean_shuffled) / std_shuffled
- FDR correction (Benjamini-Hochberg)
- **Per-kmer**: kmer, count_real, count_shuffled_mean, count_shuffled_std, z_score, p_value, fdr

#### 2.3 Motif Positional Distribution
- For top enriched k-mers: where in the UTR do they occur?
- Decile-binned density, same framework as WS1 positional analysis
- **Per-kmer x decile**: kmer, decile, count, density

#### 2.4 Motif Conservation
- For positions matching top enriched k-mers: average phyloP score
- Compare to: random UTR positions, TE-hit positions without the motif
- Uses WS6 `phylop_by_hit.tsv`
- **Per-kmer**: kmer, mean_phylop_at_motif, mean_phylop_at_random, enrichment_ratio

#### 2.5 Motif x Gene Set Association
- For each top k-mer x gene set: is the motif enriched in that gene set's UTRs?
- Fisher's exact test
- **Per-kmer x gene_set**: kmer, gene_set, odds_ratio, p_value, fdr

### Outputs

| File | Granularity |
|------|-------------|
| `shuffled_blast/` | 10 TSV files (replicate BLAST results) |
| `motif_enrichment_6mer.tsv` | Per k-mer (4,096 rows) |
| `motif_enrichment_top50.tsv` | Top 50 k-mers (quick view) |
| `motif_positions.tsv` | Per k-mer x decile (top 50, 500 rows) |
| `motif_conservation.tsv` | Per k-mer (top 50) |
| `motif_gene_set_association.tsv` | Per k-mer x gene set (top 50 x 5, 250 rows) |

### Code Approach

- **Reuse**: Legacy `analyze_te_motifs.py` (k-mer enrichment logic, Z-score framework)
- **Reuse**: Legacy `shuffle_sequences.py` (dinucleotide shuffle, well-tested)
- **Reuse**: Legacy `analyze_motif_positions.py` (positional distribution)
- **Adapt**: Coordinate-based query IDs; legacy assumed FBtr

### Caveats

- **Compute bottleneck**: Shuffled BLAST is expensive. For initial pass, consider: (a) running shuffles at evalue=0.001 instead of 10, (b) sampling 5K queries, (c) using only 3 replicates. Document the choice.
- K-mer analysis at k=6 gives 4,096 possible k-mers. Many will have zero counts. Report only k-mers with >= 10 occurrences.
- Dinucleotide shuffle preserves GC content and dinucleotide frequency but not higher-order structure. More sophisticated shuffles (e.g., preserving trinucleotide) exist but are computationally expensive.
- Pre-BLAST sequence dedup means some isoform-specific UTR variants are collapsed. If isoform-level motif analysis is needed, re-extract without dedup.

---

## Output Structure

```
results/dmel_3utr_e10/
+-- blast_results.tsv              (pipeline, existing)
+-- gene_stats.tsv                 (pipeline, existing)
+-- regions.fa, regions.tsv        (pipeline, existing)
+-- enrichment.json                (pipeline, existing)
+-- family_stats.tsv               (pipeline, existing)
+-- strand_bias*.tsv               (pipeline, existing)
+-- summary.json                   (pipeline, existing)
|
+-- ws1_match_patterns/
|   +-- hit_quality_tiers.tsv      (LARGE: ~2GB, per-hit)
|   +-- alignment_stats_summary.tsv
|   +-- positional_profiles_utr.tsv
|   +-- positional_profiles_te.tsv
|   +-- te_hit_breadth.tsv
|   +-- hit_multiplicity.tsv
|
+-- ws3_te_families/
|   +-- family_ranking_global.tsv
|   +-- family_enrichment_by_set.tsv
|   +-- class_distribution.tsv
|   +-- family_strand_bias.tsv
|   +-- te_id_to_family_class.tsv
|
+-- ws5_repeatmasker/
|   +-- known_hits.tsv
|   +-- novel_hits.tsv
|   +-- rm_overlap_by_tier.tsv
|   +-- novel_vs_known_stats.tsv
|   +-- rm_stats.json
|
+-- ws6_conservation/
|   +-- phylop_by_hit.tsv          (LARGE: per-hit with scores)
|   +-- phylop_by_quality.tsv
|   +-- synteny_by_hit.tsv
|   +-- synteny_by_quality.tsv
|   +-- synteny_by_species.tsv
|   +-- known_vs_novel_conservation.tsv
|   +-- pident_vs_conservation.tsv
|   +-- quality_paradox_stats.json
|
+-- ws7_functional/
|   +-- go_enrichment.tsv
|   +-- go_enrichment_top50.tsv
|   +-- expression_correlation.tsv
|   +-- tissue_enrichment.tsv
|   +-- flyfish_enrichment.tsv
|   +-- te_dense_genes.tsv
|
+-- ws4_unusual/
|   +-- outlier_hits.tsv
|   +-- multifamily_hotspots.tsv
|   +-- clustered_hits.tsv
|   +-- antisense_only_genes.tsv
|   +-- novel_hq_candidates.tsv
|
+-- ws2_motifs/
|   +-- shuffled_blast/            (10 replicate files)
|   +-- motif_enrichment_6mer.tsv
|   +-- motif_enrichment_top50.tsv
|   +-- motif_positions.tsv
|   +-- motif_conservation.tsv
|   +-- motif_gene_set_association.tsv
|
+-- summary_report.md
```

Same structure under `results/dmel_5utr_e10/`.

---

## Data Validation Checklist

| Source | File | Concern | Action | Status |
|--------|------|---------|--------|--------|
| GO annotations | `gene_association.gaf` | May be outdated | Verify FlyBase release, re-download if pre-r6.60 | TODO |
| RNA-seq RPKM | `gene_rpkm_report.tsv` | Static FlyBase data | Verify present, check gene count | TODO |
| FlyFISH | `flyfish_localization.csv` | Completed experiment, static | Verify present | TODO |
| RepeatMasker | `dm6.fa.out` | UCSC source, stable | Verified (137,555 repeats) | OK |
| phyloP | `dm6.phyloP27way.bw` | UCSC, stable | Verified (396 MB) | OK |
| MAF synteny | `maf/*.maf.gz` | UCSC 27-way | Verified (6 chrom files) | OK |
| TE metadata | `te_annotations.gff` | May not cover all 5,612 families | Check coverage, build lookup | TODO |
| Gene symbols | `fbgn_to_symbol.tsv` | 8,580 genes, may be incomplete | Check coverage vs 13K | TODO |
| TE consensus | `dmel_te_consensus.fasta` | 127 TEs | Verify headers have family/class | TODO |
| TE instances | `dmel_te_flybase.fasta` | 5,734 TEs | Verify FBti ID coverage | TODO |

---

## Global Caveats for Future Iterations

1. **evalue=10 saturates TE-positive gene counts** (98.4% at 3'UTR, 92.2% at 5'UTR). Binary presence/absence enrichment is nearly meaningless. All workstreams use density percentiles or quality tiers instead.

2. **max_target_seqs=1000 clips results for 308/21,661 queries**. Some true hits may be missing for queries that match many TEs. Future: raise to 5000 or use `-subject` mode for small databases.

3. **Coordinate-based query IDs** differ from legacy FBtr IDs. All legacy scripts that assume FBtr format need adaptation. The new format (`chrom:start-end`) actually simplifies genomic coordinate derivation.

4. **Pre-BLAST sequence dedup** collapses identical isoform UTRs. Analyses requiring isoform resolution need re-extraction. This affects WS2 (motif analysis) if isoform-specific patterns are suspected.

5. **hit_quality_tiers.tsv will be ~2GB** for 3'UTR at evalue=10. Consider parquet or chunked processing for future iterations. For this pass, TSV is fine since we're reading with pandas.

6. **phyloP 27-way includes vertebrate outgroups**. Scores are inflated for deeply conserved genes. Per-species Drosophilid synteny is more specific. Report both.

7. **Shuffled controls for motif analysis are compute-intensive** (~16h per UTR type at evalue=10). Initial pass should use a reduced setting (evalue=0.001 or sampled queries). Document the choice and plan full runs for final analysis.

8. **TE family assignment** depends on sseqid format (FBti instance vs consensus name). A unified lookup table is required before WS3 can produce clean class-level results. Build this as a pre-analysis step.

9. **Gene symbol coverage** may be incomplete (~8.5K of ~13K genes). Some outputs will have FBgn IDs without symbols. This is cosmetic but affects readability of reports.

10. **3'UTR vs 5'UTR comparison** should be explicit in the summary report. All workstreams run on both UTR types; differences (sense bias 63% vs 69%, hit density, positional profiles) are biologically meaningful.
