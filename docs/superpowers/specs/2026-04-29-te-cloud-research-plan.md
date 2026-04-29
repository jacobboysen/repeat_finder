# TE Cloud — Phylogenetic Ghost TE Detection

**Date:** 2026-04-29
**Status:** Paused — memorializing for future pickup
**Project:** AlphaGenome / TE Panel

## Vision

Map "last common ancestor" TEs across Drosophilid genomes. TE sequences may have
been excised from some lineages but persist as phylogenetic ghosts in others. By
building cross-species TE databases, we can detect fossils invisible from any single
genome's TE catalog.

The ultimate question: can TE-derived sequences in 3'UTRs be traced to functional
roles (condensate formation, regulatory scaffolding) that explain their conservation
across ~25-40 My of evolution?

## Completed Work

### Comparative pipeline (dmel x dpse, 2026-03-15 through 2026-03-21)

Full pipeline in `workflows/Snakefile` with cross-BLAST support. Results in
`results/comparative/`. Standalone analysis scripts in `scripts/comparative/`.

**Key results:**
- 4,560 conserved (gene, TE) pairs across orthologous dmel/dpse UTRs
- 4,330 (95%) hit the same TE region — insertions predate ~25 Mya split
- 398 "both_conserved" hits (>=80% identity in both species) — strongest
  functional constraint candidates
- 671 dmel UTRs match dpse TEs but NOT dmel TEs — ghost TE fossils
- Symmetric degradation pattern: mean pident dmel 73.1% vs dpse 73.2%,
  query match rate 0.423 (independent drift from shared ancestor)

**Infrastructure built:**
- UTR inference from mRNA-CDS for NCBI Gnomon annotations (`gff.py:infer_utrs`)
- Cross-BLAST Snakemake rule for pairwise genome comparisons
- Simple_repeat/Low_complexity filtering (26x noise reduction)
- Vectorized RM overlap analysis (was O(N*M), now vectorized)
- NCBI ortholog table: 11,515 dmel-dpse pairs (`data/references/orthologs/`)
- Genome configs: `dpse_MV25.yaml`, `dmel_r6.66_rmmode.yaml`

**Spec:** `docs/superpowers/specs/2026-03-14-comparative-te-pipeline-design.md`
**Run log:** `docs/run_logs/2026-03-15-comparative-pipeline-run.md`

### Methodological learnings

Instance-level TE databases (genomic copies) have critical limitations for
cross-species comparison:
- Same TE family across species creates redundant hits, not independent signal
- False positive inflation from searching multiple databases without correction
- RM sensitivity asymmetry (different params/libraries) confounds comparisons
- A UTR matching species B's TE copy but not species A's may be sampling
  artifact, not biology

**Better approaches identified:**
1. **Dfam curated consensus HMMs + nhmmer** — one consensus per TE family,
   HMM sensitivity > BLAST for diverged homology, clean statistics
2. **UCSC multiz27way MAF alignment** — map RM annotations across species via
   syntenic alignment, sidesteps BLAST entirely for ghost TE detection
3. **Pan-Drosophila consensus library** — RepeatModeler2 de novo per genome,
   cluster with CD-HIT-EST at 80%, curate

## Open Research Threads

### 1. Categorical Jacobians for TE identification

Can Jacobian matrices from genomic language models (NTv3, AlphaGenome) identify
TE-like sequences? The Jacobian captures positional sensitivity — TE fossils may
show distinct patterns vs surrounding UTR (structural breaks at TE boundaries,
different selection signatures).

**How to start:** Run NTv3 (already configured in pipeline as
`InstaDeepAI/NTv3_100M_post`) on 5-10 UTRs containing "both_conserved" TE
fossils. Extract per-position Jacobians. Visualize: do TE boundaries show
discontinuities? Do conserved positions within TE fossils have different
sensitivity profiles than drifting positions?

**Config:** `src/fossil_finder/config/genomes/dmel_r6.66.yaml` has NTv3 scoring
spec (model, species_key, window_size).

### 2. Condensate-forming RNAs and TEs

Recent ULE paper on condensate-forming RNAs — could TE-derived sequences provide
scaffolding for RNA phase separation? This would give a mechanistic explanation
for why certain TE fossils are retained in UTRs across 25 My.

**How to start:** Get ULE paper's condensate-forming RNA predictions or sequence
features. Intersect with coordinates of our 398 both_conserved TE fossils. Test
enrichment: are TE fossils more likely to overlap condensate-forming regions than
matched non-TE UTR sequence?

### 3. Kmer spectral matrices for TE classification

(Discussed 2026-03-29) Build matrices of kmer frequencies to find conserved and
distinct spectra across genomic regions. Key questions:
- Do certain regions (UTRs, introns, promoters) absorb different TE kmer signatures?
- Can distinct biological capabilities (e.g., regeneration) be traced to unique
  kmer spectra +/- TE?
- Can spectral analysis distinguish functional TE fossils from neutral insertions?

**How to start:** Build kmer (k=6) frequency matrices for: (a) known TE families,
(b) UTR regions with TE hits, (c) UTR regions without. PCA/UMAP to find clusters.
Compare spectra across dmel/dpse to identify conserved signatures. The pipeline's
motif analysis module (`analysis/motifs.py:count_kmers_from_blast`) already has
kmer counting infrastructure.

## Highest Impact Next Steps (ranked)

1. **Categorical Jacobian pilot** — lowest effort, highest novelty. NTv3
   infrastructure exists. Just need to extract Jacobians and visualize.

2. **Dfam consensus pivot** — replace instance DBs with Dfam 3.9 curated
   consensus (~400 families). Cleaner statistical framework for all downstream
   analysis. Required before expanding to more species.

3. **Kmer spectral analysis** — quick PCA on 398 both_conserved hits. Can be
   done with existing data, no new BLAST runs needed.

4. **Condensate RNA cross-reference** — depends on getting ULE paper data.
   High impact if positive.

5. **Add D. virilis** — most informative species addition (novel TE families:
   PLEs, Helitrons, unusual TE abundance). Pipeline infrastructure ready.

## Species Priority for Expansion

| Priority | Species | Divergence | Unique value |
|----------|---------|-----------|-------------|
| Done | D. pseudoobscura | ~25 Mya | Sweet spot distance, good ortholog coverage |
| 1 | D. virilis | ~40 Mya | Novel TE families (PLEs, Helitrons), high TE content |
| 2 | D. yakuba | ~6 Mya | Close-range validation, HT detection |
| 3 | D. ananassae | ~40 Mya | Different TE landscape (LTR-dominant, unique hAT) |

## Known Issues to Fix

- dmel RM overlap shows `known_hits: 0` — chr prefix mismatch in `find_overlaps`
  (RM uses `chr2L`, extracted regions use `2L`). Same fix pattern as
  `build_te_fasta_from_repeatmasker`.
- dpse UTR deduplication: 25,693 regions but only ~16,029 unique genomic
  coordinates (~9,664 duplicates from multi-isoform inference).
- Cross-BLAST analysis steps not yet run (only self-BLAST analyze completed).

## Key Files

| File | Purpose |
|------|---------|
| `workflows/Snakefile` | Pipeline with cross_blast rule |
| `src/fossil_finder/config/genomes/dpse_MV25.yaml` | dpse genome config |
| `src/fossil_finder/config/genomes/dmel_r6.66_rmmode.yaml` | dmel RM-mode config |
| `src/fossil_finder/io/gff.py` | UTR inference (infer_utrs) |
| `src/fossil_finder/analysis/repeatmasker.py` | Vectorized overlap/classify |
| `src/fossil_finder/repeatmasker/library.py` | TE extraction + simple repeat filter |
| `scripts/comparative/conserved_te_hits.py` | Conserved TE analysis |
| `scripts/comparative/dpse_only_hits.py` | Lineage-specific TE analysis |
| `scripts/comparative/alignment_degradation.py` | Cross-species alignment comparison |
| `data/references/orthologs/dmel_dpse_orthologs.tsv` | 11,515 ortholog pairs |
| `results/comparative/analysis/` | Analysis outputs |
