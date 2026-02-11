# NTv3 TE Fossil Regulatory Analysis — Summary

## Approach

We scored 1,900 TE-bearing regulatory regions (promoters and enhancers) in *Drosophila melanogaster* using Nucleotide Transformer v3 (NTv3), a genomic language model trained on multi-species genomic data. For each region, we computed:

- **te_enrichment**: ratio of NTv3 predicted regulatory signal at TE-derived positions vs non-TE positions (>1 means TEs carry more signal)
- **ablation_delta**: change in predicted signal when TE positions are masked with Ns (positive means TEs contribute to regulation)

TE positions were identified by BLAST against the FlyBase *D. melanogaster* TE database (strict tier: e-value ≤ 1e-5), and cross-referenced with RepeatMasker annotations. Conservation was measured using phyloP 27-way Drosophila scores from UCSC.

Regions were ranked by TE base-pair coverage (te_bp_strict) and scored on an A100 GPU via Lightning.ai.

## Key Numbers

| Metric | Value |
|--------|-------|
| Regions scored | 1,900 |
| Functional (enrich > 1 AND ablation_delta > 0) | 651 (40.4%) |
| Mean te_enrichment | 3.43 |
| Mean ablation_delta | +0.0105 |
| True exaptation candidates (functional + conserved) | 145 |
| BLAST-only regions scored (no RepeatMasker overlap) | 18 |

## Finding 1: The Ablation Paradox

**More degraded TE fragments are more load-bearing for regulatory predictions.**

| Measure | rho | p |
|---------|-----|---|
| pident vs ablation_delta (raw) | -0.114 | 9.8e-7 |
| pident vs ablation_delta \| te_fraction | -0.106 | 5.9e-6 |
| pident vs ablation_delta \| te_fraction + n_te_positions | -0.107 | 5.1e-6 |

The relationship survives controlling for TE coverage (te_fraction) and number of TE positions. pident and te_fraction are essentially uncorrelated (rho = +0.04, ns), ruling out the mechanical confound that lower-pident regions simply have more TE to mask.

**The enrichment-ablation dissociation is the mechanistic insight:**

| Correlation | rho | p | Interpretation |
|-------------|-----|---|----------------|
| pident vs te_enrichment | +0.02 | 0.32 (ns) | TEs carry the same signal density regardless of age |
| pident vs ablation_delta | -0.11 | 9.8e-7 | But removing degraded TE signal disrupts predictions more |

TEs arrive with regulatory potential (consistent with Du et al. 2024's "ancestral origin" model for LTR elements). Recent insertions are modular — the surrounding sequence compensates when they're removed. Degraded fragments have become integrated into the regulatory grammar: flanking elements have co-adapted to depend on them. This is what exaptation looks like at the sequence level.

**The effect is strongest in low-TE regions** (Q1 te_fraction quartile: rho = -0.23, p < 0.001), where a small TE fragment sits within a larger regulatory context. At high TE fractions (Q4), where the TE essentially IS the promoter, the effect washes out (rho = -0.05, ns). This is consistent with the co-adaptation model: integration requires surrounding non-TE sequence to co-evolve with the TE fragment.

## Finding 2: Conservation Validates Model Predictions

PhyloP conservation (27-way Drosophila alignment) and NTv3 functional predictions converge on the same regions despite being methodologically independent.

| Metric | Value |
|--------|-------|
| Fisher's exact OR (conserved × functional) | 3.83 |
| p-value | 9.5e-19 |
| Functional rate among conserved TE regions | 68.4% |
| Functional rate among non-conserved | 36.1% |
| te_enrichment vs phyloP_te (Spearman) | rho = +0.32, p = 9.5e-39 |

**Conservation × Function matrix (N = 1,613 with phyloP data):**

| | Functional | Non-functional |
|---|---|---|
| **Conserved (phyloP > 1)** | 145 (9.0%) — true exaptation | 67 (4.2%) — structural constraint? |
| **Not conserved** | 506 (31.4%) — species-specific or young | 895 (55.5%) — neutral passengers |

TE positions are on average less conserved than non-TE positions within the same regions (0.42 vs 0.61, Wilcoxon p = 9.5e-39), consistent with most TE insertions being neutral or deleterious. The 145 conserved-and-functional regions represent the subset where purifying selection has maintained TE-derived regulatory function.

## Finding 3: BLAST-Only Regions

Our BLAST approach detects 2,735 regulatory regions with TE signal invisible to RepeatMasker (18.9% of all TE-bearing regions). These are more diverged (mean pident 69.7% vs 77.1% for RM-confirmed) and dominated by ancient TE families: INE-1 (SINE), roo (LTR), 17.6 (LTR).

18 BLAST-only regions were scored (they first appear at pipeline rank 720). Highlights:

| Region | Gene | TE Family | Pident | Enrichment | Ablation | PhyloP | Status |
|--------|------|-----------|--------|------------|----------|--------|--------|
| TSS_mE1_002238 | CG10465 | 412 (LTR) | 69.4% | 2.11 | +0.019 | 1.17 | CONS + FUNC |
| TSS_mE1_002239 | CG10465 | 412 (LTR) | 72.5% | 2.14 | +0.020 | 0.87 | FUNC |
| TSS_RAMPAGE_010361 | CG10465 | Burdock (LINE) | 67.8% | 2.22 | +0.018 | 0.84 | FUNC |
| TSS_mE1_002255 | Not3 | diver2 | 69.3% | 1.77 | +0.016 | 0.63 | FUNC |
| TSS_mE1_003730 | CG4866 | FB | 61.5% | 0.74 | +0.014 | 2.15 | CONS only |
| FBsf0000872982 | Cyp6t1 | Y | 87.4% | 1.67 | +0.009 | 0.50 | FUNC |

CG10465 encodes a BTB/POZ domain protein (KCTD family); human orthologs KCTD10/KCTD13 are implicated in neurodevelopmental disorders (16p11.2 locus). The gene is maternally deposited. Three independent promoter variants show TE-derived regulatory signal from degraded 412/Burdock elements invisible to RepeatMasker.

## Finding 4: Quality Paradox Does Not Hold for pident

Earlier 3'UTR analysis suggested that stricter BLAST quality filters yield less conserved hits ("Quality Paradox"). In promoters/enhancers, the direction is **opposite**: pident vs phyloP_te is positive (rho = +0.19, p < 1e-14). More recognizable TEs are more conserved. The 3'UTR result was likely driven by the length filter (hit_length vs phyloP: rho = +0.32), not sequence identity.

The ablation paradox (Finding 1) is the non-circular reformulation: the relevant signal isn't "degraded TEs are more conserved" but "degraded TEs are more integrated into regulatory grammar."

## Top Exaptation Candidates

The 25 highest-confidence candidates (conserved phyloP > 1, functional enrichment > 1, positive ablation delta), ranked by phyloP:

| Rank | Gene | Region Type | PhyloP | Enrichment | Ablation Δ | Pident |
|------|------|-------------|--------|------------|------------|--------|
| 1 | eIF3f1 | promoter | 2.44 | 5.31 | +0.026 | 88.6% |
| 2 | eIF3f1 | promoter | 2.41 | 5.84 | +0.023 | 87.6% |
| 3 | scyl | enhancer | 2.31 | 1.68 | +0.056 | 64.5% |
| 4 | Unr | promoter | 2.19 | 1.17 | +0.081 | 66.1% |
| 5 | Ank | promoter | 2.12 | 1.89 | +0.008 | 81.8% |
| 6 | Mo25 | promoter | 2.06 | 1.80 | +0.013 | 86.4% |
| 7 | Not3 | promoter | 1.80 | 4.56 | +0.028 | 73.6% |
| 8 | mol | promoter | 1.79 | 25.33 | +0.003 | 96.5% |
| 9 | eIF3g1 | promoter | 1.77 | 333.3 | +0.013 | 79.5% |
| 10 | Pdzd8 | promoter | 1.76 | 1.43 | +0.014 | 70.5% |

Notable: **Unr** (upstream of N-ras, involved in dosage compensation and mRNA regulation) has the largest ablation effect in the top candidates (+0.081). **eIF3g1** has the highest enrichment (333x — region is nearly 100% TE-derived with extreme signal). **scyl** enhancer combines high conservation with a large ablation effect, suggesting a deeply integrated TE-derived enhancer element.

## Comparison to Du et al. 2024

Du et al. ("Regulatory transposable elements in the encyclopedia of DNA elements," *Nature Communications*) characterized TE-derived cis-regulatory elements in human using RepeatMasker annotations and ENCODE cCRE data. Key parallels and differences:

| Aspect | Du et al. 2024 | This analysis |
|--------|---------------|---------------|
| Species | Human | *D. melanogaster* |
| TE detection | RepeatMasker | BLAST (finds 18.9% more regions) |
| Regulatory assay | ENCODE cCREs, MPRA | NTv3 genomic language model |
| Conservation | Cross-species cCRE overlap | phyloP 27-way Drosophila |
| Key finding | TEs carry TFBSs from birth (ancestral origin) | Confirmed: ~40% functional regardless of pident |
| Novel finding | — | Ablation paradox: degraded TEs more load-bearing |
| Gap they flag | "TEs decayed past recognition... underestimating TE contribution" | We directly address this with BLAST-only regions |

Our ablation paradox (degraded TEs more integrated into regulatory grammar) is complementary to Du et al.'s ancestral origin model: TEs arrive with regulatory potential, but the ones that persist long enough to degrade become wired into their regulatory context. Du et al. describe the input; we describe the integration.

## Methods Note

The candidate classification used in the scoring pipeline (exapted_fossil/young_active/neutral) uses a pident < 80% threshold that creates circularity — the "Quality Paradox" enrichment table it generates is definitional, not biological. The non-circular analyses above (Spearman correlations, partial correlations, Fisher's exact tests) avoid this issue by treating pident, model outputs, and conservation as continuous variables.

## Files

| File | Description |
|------|-------------|
| `results/ntv3_te_scoring/region_scores.tsv` | Per-region NTv3 scores (1,900 regions) |
| `results/ntv3_te_scoring/ablation_results.tsv` | Ablation test results (3,800 rows: mask + shuffle per region) |
| `results/ntv3_te_scoring/integrated_results.tsv` | Merged scores + ablation + pident + classification |
| `results/te_fossil_lm_dataset/` | Full dataset: regions.tsv, hits.tsv, sequences.fasta, masks.npz |
| `data/references/dm6.phyloP27way.bw` | PhyloP 27-way conservation scores (UCSC) |
| `scripts/score_te_fossils_ntv3.py` | NTv3 scoring pipeline |
| `scripts/prepare_te_fossils_for_lm.py` | Dataset preparation pipeline |
