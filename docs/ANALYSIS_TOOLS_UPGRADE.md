# Analysis Tools Upgrade Reference

> Third-party tools that could enhance each workstream of the TE fossil analysis suite.
> For future iterations -- not required for the current pass.

**Date**: 2026-03-11
**Context**: These tools were evaluated against our 7-workstream architecture
(see `docs/ANALYSIS_ARCHITECTURE.md`). Prioritized by: relevance to TE fossil
analysis, maintenance status, Drosophila support, and ease of integration.

---

## WS1: Match Distributions & Patterns

### Current approach
Custom Python: qseq/sseq parsing, mismatch profiling, positional normalization.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **edlib** | Fast C++ edit distance library with Python bindings | Compute true edit distance (Levenshtein) between qseq/sseq, much faster than pure Python. Supports NW/HW/SHW modes. | `pip install edlib` | Only does edit distance, not full alignment statistics. Use alongside our mismatch parsing. |
| **parasail** | SIMD-accelerated pairwise alignment (SSW, NW, SG) | Re-align qseq/sseq with custom scoring matrices tuned for TE divergence. Could reveal whether BLAST's default scoring misses certain mutation patterns. | `conda install -c bioconda parasail-python` | Overkill for current pass. Useful if we need custom scoring matrices later. |
| **pysamstats** | BAM/alignment statistics | If we ever move to mapping-based TE detection (minimap2), this handles positional stats. Not needed for BLAST output. | `pip install pysamstats` | Only relevant for mapping-based workflows. |
| **MAFFT** | Multiple sequence alignment | Align all BLAST hit sequences (qseq) from a single TE family to visualize mutation accumulation patterns. Outperforms other MSA tools on divergent TE sequences (Hubley et al. 2022). Use L-INS-i mode for <200 sequences. | `conda install -c bioconda mafft` | Compute-intensive for large families. Sample top 100 hits per family. |
| **logomaker** | Sequence logo visualization | Generate position-specific logos from aligned TE hits to visualize which positions are conserved vs diverged within a TE family. | `pip install logomaker` | Visualization only. Useful for figures in the report. |

### Recommended upgrade path
1. Add `edlib` for fast edit distance computation (immediate value, trivial integration)
2. Add `MAFFT` + `logomaker` for per-family alignment visualization (notebook phase)

---

## WS2: Motif Analysis

### Current approach
Custom k-mer counting, dinucleotide shuffle, Z-score enrichment.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **MEME Suite** (MEME, DREME, SEA, AME, CentriMo) | De novo motif discovery, motif enrichment analysis | MEME finds novel motifs in TE-hit regions. DREME finds short (3-8bp) enriched motifs. SEA tests enrichment of known motifs. CentriMo tests positional enrichment. AME for motif enrichment with flexible statistical tests. | `conda install -c bioconda meme` | Large install. Web server at meme-suite.org for one-off analyses. Local install for batch processing. |
| **STREME** | Short motif discovery (MEME Suite successor to DREME) | Faster than DREME, finds enriched short motifs in TE-containing vs non-TE UTR sequences. Directly uses our shuffled controls as the negative set. | Part of MEME Suite | Requires MEME Suite >= 5.4.1. |
| **HOMER** | Motif discovery + enrichment (findMotifsGenome.pl) | Integrated motif discovery with background model. Has Drosophila dm6 genome built in. Can search against known motif databases. | `conda install -c bioconda homer` | Perl-based, large dependency set. Better for ChIP-seq than UTR analysis but works. |
| **Biostrings** (R/Bioconductor) | Efficient biological string operations | Fast k-mer counting, pattern matching, dinucleotide frequency. Alternative to our Python implementation. | R: `BiocManager::install("Biostrings")` | Requires R. Our Python approach is sufficient. |
| **ushuffle** | Efficient k-let preserving shuffle | Faster dinucleotide shuffle than our Python implementation. C library with Python wrapper. Preserves exact k-let frequencies (our shuffle preserves dinucleotide). | `pip install ushuffle` | Minor speedup. Our shuffle works. |

### Recommended upgrade path
1. Add **MEME Suite** (STREME + SEA) for de novo motif discovery alongside our k-mer approach
2. Use STREME output as input to SEA for enrichment testing against known Drosophila motif databases

---

## WS3: TE Family Enrichment

### Current approach
`fossil_finder.analysis.families` module: per-family stats, fold-enrichment, class distribution.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **RepeatMasker** (with RepBase/Dfam) | TE classification | Dfam database provides TE family -> superfamily -> class hierarchy. Could supplement our FBti-based classification. Dfam 3.8 has 273K TE families across species. | `conda install -c bioconda repeatmasker` | Already have dm6.fa.out. Dfam provides the classification hierarchy we need for WS3. Download Dfam HMM profiles separately. |
| **Earl Grey** | Automated TE annotation + analysis pipeline | Full TE annotation including de novo discovery, classification, and landscape analysis. Could discover TEs our database misses. Produces TE landscape divergence plots. | Docker/Singularity | Heavy pipeline. More useful for genome annotation than downstream analysis. Future iteration for re-annotating dm6 TEs. |
| **TEtools** / **TE-Aid** | TE structure visualization and validation | Visualize TE structural domains (LTR, TIR, coding regions). Helps interpret why certain TE regions are overrepresented in hits. | Various | Visualization aid, not analytical. |
| **Dfam** database | TE classification hierarchy | Provides TE family -> class -> order taxonomy. Our FBti IDs can be mapped through Dfam for systematic classification. REST API available. | Web API: dfam.org/api | FlyBase TEs are well-represented in Dfam. Use API to bulk-classify sseqids. |
| **TE Hub** | Community TE database | Curated TE classification for model organisms including Drosophila. Bridges FlyBase and RepBase nomenclature. | Web: tehub.org | Manual lookup. Useful for resolving naming conflicts between FlyBase/RepBase/Dfam. |

### Recommended upgrade path
1. Download **Dfam** classification tables for Drosophila to build the te_id_to_family_class.tsv lookup
2. Consider **Earl Grey** for future genome re-annotation if current TE database gaps are found

---

## WS4: Unusual Matches & Patterns

### Current approach
Custom outlier detection, clustering, multi-family analysis.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **scikit-learn** (IsolationForest, LOF) | Unsupervised outlier detection | IsolationForest on feature vectors (pident, length, evalue, mismatch_rate, gap_rate) to systematically identify unusual hits without manual thresholds. LOF for local outlier detection. | Already in environment.yml | IsolationForest is fast (O(n log n)). LOF is O(n^2) -- sample for 7.5M hits. |
| **BEDTools** | Genomic interval operations | `bedtools cluster` for merging nearby hits into clusters (WS4.3). `bedtools intersect` for overlap detection. Faster than pure Python for interval operations. | `conda install -c bioconda bedtools` | C++ tool, very fast. Would speed up hit clustering significantly. |
| **pybedtools** | Python wrapper for BEDTools | Pythonic interface to BEDTools. Integrates with pandas DataFrames. | `pip install pybedtools` | Requires BEDTools installed. Clean API. |
| **GenomicRanges** (R/Bioconductor) | Genomic interval algebra | findOverlaps, reduce, disjoin for interval operations. Alternative to BEDTools with richer semantics. | R: `BiocManager::install("GenomicRanges")` | Requires R. BEDTools/pybedtools is sufficient in Python. |

### Recommended upgrade path
1. Add **BEDTools** + **pybedtools** for interval clustering (immediate value for WS4.3)
2. Use **IsolationForest** from scikit-learn for systematic outlier detection (already available)

---

## WS5: RepeatMasker Overlap

### Current approach
`fossil_finder.analysis.repeatmasker` module: parse .out file, interval overlap, classify known/novel.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **BEDTools** intersect | Fast genomic interval intersection | Replace our Python overlap detection with `bedtools intersect -a blast_hits.bed -b rm.bed -wo`. Orders of magnitude faster for 7.5M hits x 137K RM annotations. | `conda install -c bioconda bedtools` | Would dramatically speed up WS5. Our Python is correct but slow at scale. |
| **GAT** (Genomic Association Tester) | Permutation-based overlap significance | Tests whether BLAST hit overlap with RepeatMasker is more/less than expected by chance, accounting for genome organization (isochores, mappability). | `pip install gat` | Adds statistical rigor to the overlap comparison. Currently we just count; GAT gives p-values. |
| **UCSC liftOver** | Coordinate conversion between assemblies | If we need to compare RM annotations across genome versions or between species. | `conda install -c bioconda ucsc-liftover` | Only needed if working across assemblies. |
| **RepeatMasker** (re-run) | Custom RM analysis | Re-run RepeatMasker with sensitive settings on our UTR sequences specifically. Our dm6.fa.out is genome-wide with default sensitivity. | `conda install -c bioconda repeatmasker` | Compute-intensive. The existing dm6.fa.out is standard and sufficient. |

### Recommended upgrade path
1. Add **BEDTools** for fast interval intersection (immediate value)
2. Add **GAT** for permutation-based significance testing (enriches the novel-hit narrative)

---

## WS6: Conservation, Synteny & PhyloP

### Current approach
bigWigAverageOverBed for phyloP, custom MAF parsing for synteny.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **PHAST** (phastCons, phyloP, phyloFit) | Phylogenetic conservation scoring | We already use phyloP scores from UCSC. PHAST lets us compute custom conservation scores with our own species tree, or test specific evolutionary models (acceleration, conservation). `phyloP` can test per-element acceleration (faster evolution than neutral). | `conda install -c bioconda phast` | Could test whether TE fossils show accelerated evolution (positive selection after exaptation). Powerful but requires careful model specification. |
| **HAL tools** (halLiftover, halStats) | Hierarchical alignment format tools | Convert between MAF and HAL for more efficient random access to synteny data. HAL supports progressive alignment queries. | `conda install -c bioconda hal` | Only useful if we switch from MAF to HAL format. Current MAF approach works. |
| **GERP++** | Constrained element detection | Identifies genomic elements under evolutionary constraint. Unlike phyloP (per-base), GERP identifies constrained *elements* (contiguous runs of conservation). | Standalone binary | Would identify which TE fossils fall within constrained elements rather than just per-base scores. Complements phyloP. |
| **Cactus** / **Progressive Cactus** | Whole-genome alignment | State-of-the-art multiple genome aligner. Could produce higher-quality alignments than UCSC multiz for close Drosophilid comparisons. | Docker/Singularity | Very heavy pipeline. Only worth it if UCSC alignments are insufficient. |
| **pyBigWig** | Python interface to bigWig files | Faster and more flexible than bigWigAverageOverBed for per-hit phyloP extraction. Can extract per-base scores for detailed conservation profiles along TE hits. | `pip install pyBigWig` | Direct Python integration. Could replace the shell-out to bigWigAverageOverBed. |
| **UCSC 124-insect alignment** | Expanded phyloP for insects | UCSC now provides 124-insect multiz alignments and phyloP scores for dm6. Much broader species tree than the 27-way we currently use. | Download from UCSC | Larger file, more species. Better power for detecting insect-specific conservation vs vertebrate-contaminated signal in 27-way. |

### Recommended upgrade path
1. Add **pyBigWig** for direct Python phyloP extraction (replaces shell-out, faster iteration)
2. Download **UCSC 124-insect phyloP** for insect-specific conservation (more relevant than 27-way)
3. Consider **PHAST phyloP acceleration test** for testing whether TE fossils evolve faster/slower than neutral

---

## WS7: GO & Functional Enrichment

### Current approach
Custom Fisher's exact test + BH correction, legacy annotation loaders.

### Tools to consider

| Tool | What it does | How it helps us | Install | Caveat |
|------|-------------|-----------------|---------|--------|
| **goatools** | Python GO enrichment analysis | Full GO DAG traversal, Fisher's exact + multiple testing, GO slim mapping. More robust than our custom implementation -- handles GO hierarchy (propagation), which our flat approach misses. | `pip install goatools` | Pure Python, well-maintained (>1K GitHub stars). Handles GO term hierarchy properly. |
| **g:Profiler** (gprofiler2) | Multi-source enrichment (GO, KEGG, Reactome, TF, miRNA) | Web API + Python client. Tests enrichment against multiple databases simultaneously. Has Drosophila melanogaster support. | `pip install gprofiler-official` | Web API dependency. Local mode available. Broader than GO-only analysis. |
| **clusterProfiler** (R/Bioconductor) | GO/KEGG enrichment + visualization | Gold standard for GO enrichment in R. Rich visualization (dot plots, gene-concept networks, enrichment maps). Has Drosophila annotation package (org.Dm.eg.db). | R: `BiocManager::install("clusterProfiler")` | Requires R. Best visualization of any GO tool. Consider for notebook figures. |
| **DAVID** | Functional annotation clustering | Groups enriched terms by semantic similarity. Reduces redundancy in GO results. | Web: david.ncifcrf.gov | Web-only. Gene list size limits. Outdated relative to goatools/g:Profiler. |
| **FlyEnrichr** | Drosophila-specific enrichment | Enrichr variant trained on Drosophila gene sets and pathways. | Web: maayanlab.cloud/FlyEnrichr | Web-only. Useful for quick lookups. |
| **rrvgo** (R) | Reduce GO term redundancy | Semantic similarity-based clustering of enriched GO terms. Treemap visualization. Eliminates the "wall of GO terms" problem. | R: `BiocManager::install("rrvgo")` | Requires R. Post-processing step after enrichment. |

### Recommended upgrade path
1. Replace custom Fisher test with **goatools** for proper GO hierarchy handling (immediate improvement)
2. Add **g:Profiler** for multi-database enrichment (GO + KEGG + TF targets in one call)
3. Use **clusterProfiler** / **rrvgo** in notebooks for publication-quality GO visualization

---

## Cross-Cutting Tools

These tools benefit multiple workstreams.

| Tool | Workstreams | What it adds | Install |
|------|-------------|-------------|---------|
| **BEDTools** + **pybedtools** | WS4, WS5, WS6 | Fast genomic interval operations. Replaces slow Python interval overlap code. | `conda install -c bioconda bedtools && pip install pybedtools` |
| **pyBigWig** | WS6, WS2 | Direct Python access to bigWig files (phyloP). Replaces shell-out. | `pip install pyBigWig` |
| **polars** | All | DataFrame library 10-50x faster than pandas for large files. Our 7.5M-hit TSVs would benefit. | `pip install polars` |
| **pyarrow** / **parquet** | WS1, WS5, WS6 | Columnar storage for large per-hit files. 3-5x smaller than TSV, much faster I/O. | `pip install pyarrow` |
| **plotnine** / **seaborn** | All | Grammar-of-graphics plotting for Python. Publication-quality figures for notebooks. | `pip install plotnine` (or seaborn already installed) |
| **Dfam API** | WS3, WS4 | TE classification hierarchy. REST API for bulk TE ID -> family -> class lookup. | Web API: dfam.org/api/families |

---

## TE-Specific Research Tools

These are specialized tools for TE exaptation analysis, published in the literature.

| Tool / Method | Citation | Relevance |
|--------------|----------|-----------|
| **TEtranscripts** | Jin et al. 2015, Bioinformatics | Quantifies TE expression from RNA-seq. Could correlate TE fossil density with TE transcription. |
| **SQuIRE** | Yang et al. 2019, NAR | Locus-specific TE expression quantification. More precise than TEtranscripts for individual TE copies. |
| **TERA** (TE Regulatory Activity) | Sundaram et al. 2018, MBE | Tests whether TE copies in regulatory regions show stronger conservation than expected. Directly relevant to our quality paradox (WS6). Key finding: "exaptation is rare, influenced by evolutionary age, and subject to pleiotropic constraints." |
| **Kimura distance landscape** | Standard TE analysis | Plot TE copy divergence from consensus as a proxy for insertion age. Our BLAST pident approximates this. |
| **Earl Grey** | Baril et al. 2024, MBE | Full TE annotation pipeline with landscape analysis, clustering, and divergence estimation. |

---

## Drosophila-Specific Databases

| Resource | URL | Contents | Use case |
|----------|-----|----------|----------|
| **FlyBase** | flybase.org | Definitive Drosophila annotation. Gene models, TE annotations, expression, GO. | Primary reference (already used). |
| **Dfam** | dfam.org | TE HMM profiles + classification hierarchy. 273K families. Drosophila well-represented. | TE family -> class mapping (WS3). |
| **FlyTF** | flytf.org | Drosophila transcription factor database. | TF binding in TE regions (WS2 extension). |
| **REDfly** | redfly.ccr.buffalo.edu | Regulatory Element Database for Fly. Experimentally validated CRMs. | Overlap with TE fossils (WS4 extension). |
| **modENCODE** | modencode.org | Drosophila chromatin states, histone marks, TF binding. | Chromatin context of TE fossils (future). |
| **UCSC 124-insect alignment** | hgdownload.soe.ucsc.edu | Expanded multiz alignment + phyloP for 124 insects. | Better conservation analysis (WS6). |

---

## Priority Installation Order

For the next iteration, install these tools in this order (highest ROI first):

```bash
# Tier 1: Immediate analytical improvements
pip install edlib pyBigWig goatools pybedtools
conda install -c bioconda bedtools meme

# Tier 2: Performance improvements
pip install polars pyarrow

# Tier 3: Extended analysis capabilities
pip install gprofiler-official logomaker
conda install -c bioconda mafft phast

# Tier 4: Future research directions (evaluate first)
# Earl Grey (Docker), TERA method (custom), 124-insect phyloP (download)
```

---

## Notes

- All tools listed are free for academic use. MEME Suite and RepeatMasker have
  non-commercial licenses. Others are MIT/BSD/Apache.
- Our current Python implementation is correct and sufficient for the first pass.
  These tools add speed, statistical rigor, or capabilities we don't have.
- The biggest single upgrade would be **BEDTools + pybedtools** -- it speeds up
  interval operations across WS4, WS5, and WS6 by orders of magnitude.
- The biggest analytical upgrade would be **goatools** -- proper GO hierarchy
  handling is important for interpretable enrichment results.
