# Data Files and Visualizations Reference

This document provides a complete audit of all data files and visualization outputs in the `repeat_finder` project.

---

## Data Directory Structure

```
data/
├── references/          # Downloaded reference files (FlyBase, Dfam)
├── blastdb/             # BLAST nucleotide databases
├── gene_lists/          # Curated gene lists for analysis
├── queries/             # Extracted query sequences for BLAST
├── te_clusters/         # Detected TE signal clusters
└── cache/               # Cached API responses
```

---

## Reference Files (`data/references/`)

These are the core reference files downloaded from FlyBase (release r6.66) and Dfam.

| File | Size | Sequences | Description |
|------|------|-----------|-------------|
| `dmel_genome.fasta` | 139 MB | 1,870 | D. melanogaster genome (chromosomes + scaffolds) |
| `dmel_annotation.gff` | 6.3 GB | — | Full genome annotation (genes, transcripts, etc.) |
| `dmel_3utr.fasta` | 22 MB | **30,324** | 3'UTR sequences (primary query sequences) |
| `dmel_5utr.fasta` | 14 MB | ~28,000 | 5'UTR sequences |
| `dmel_cds.fasta` | 71 MB | ~31,000 | Coding sequences |
| `dmel_intron.fasta` | 129 MB | ~60,000 | Intron sequences |
| `dmel_te_flybase.fasta` | 9.4 MB | **5,734** | FlyBase TE insertions (primary TE database) |
| `dmel_te_consensus.fasta` | 628 KB | **127** | Bergman Lab TE consensus sequences |
| `te_annotations.gff` | 188 KB | — | TE structural annotations (LTR, CDS regions) |
| `gene_groups.tsv` | 980 KB | — | FlyBase gene group memberships |
| `gene_info.tsv` | 1.2 MB | — | FBgn/FBtr/FBpp ID mappings |

### Key Statistics

- **30,324** 3'UTR sequences available for analysis
- **5,734** individual TE insertions in FlyBase database
- **127** TE consensus sequences (family representatives)

---

## BLAST Databases (`data/blastdb/`)

Pre-built BLAST nucleotide databases for searching.

| Database | Source | Sequences | Description |
|----------|--------|-----------|-------------|
| `dmel_genome` | dmel_genome.fasta | 1,870 | Full genome search |
| `dmel_te_flybase` | dmel_te_flybase.fasta | 5,734 | **Primary TE database** |
| `dmel_te_consensus` | dmel_te_consensus.fasta | 127 | TE family representatives |
| `dmel_te_combined` | Combined FASTA | ~5,861 | FlyBase + consensus merged |

### Database Files

Each database consists of multiple files:
- `.nhr` - Header information
- `.nin` - Index
- `.nsq` - Sequence data
- `.ndb`, `.nog`, `.nos`, `.not`, `.nto`, `.ntf`, `.njs` - Auxiliary files

---

## Gene Lists (`data/gene_lists/`)

Curated gene lists for comparative analysis.

### Germ Plasm Genes

| File | Description |
|------|-------------|
| `germ_plasm_genes_consolidated.tsv` | Full list with FBgn IDs, symbols, functions |
| `germ_plasm_genes_tier1.txt` | Core genes only: nos, osk, pgc, gcl, vas, aub, CycB |
| `germ_plasm_fbgn_ids.txt` | FBgn IDs only (for scripting) |
| `gene_list_status.json` | Build metadata |

**Tier 1 genes:** nos, osk, pgc, gcl, vas, aub, CycB (7 genes)
**Tier 2 genes:** tudor, piwi, AGO3, dhd, Hsp83 (5 genes)
**Total:** 12 genes

### Housekeeping Genes (Control)

| File | Description |
|------|-------------|
| `housekeeping_genes_consolidated.tsv` | Full list with metadata |
| `housekeeping_fbgn_ids.txt` | FBgn IDs only |
| `housekeeping_status.json` | Build metadata |

**Genes:** Act5C, RpL32, Gapdh1, alphaTub84B, Ef1alpha48D, RpS17, Tbp, eIF4A, SdhA, Rpl4 (10 genes)

### Control Gene Sets

| File | Description | Genes |
|------|-------------|-------|
| `somatic_genes_consolidated.tsv` | Somatic-specific genes | 8 genes (en, hh, dpp, etc.) |
| `somatic_fbgn_ids.txt` | FBgn IDs only | |
| `cleared_genes_consolidated.tsv` | Maternally-cleared mRNAs | 10 genes (Hsp83, stg, etc.) |
| `cleared_fbgn_ids.txt` | FBgn IDs only | |
| `adult_genes_consolidated.tsv` | Adult-specific genes | 9 genes (Adh, Lsp1, etc.) |
| `adult_fbgn_ids.txt` | FBgn IDs only | |
| `control_lists_status.json` | Build metadata | |

---

## Query Sequences (`data/queries/`)

Extracted sequences ready for BLAST searching.

### Per-Gene-Set Structure

Each gene set has its own directory with consistent file structure:

```
data/queries/{gene_set}/
├── 3UTR_sense.fasta           # Sense strand sequences
├── 3UTR_antisense.fasta       # Reverse complement (control)
├── 3UTR_sense_tier1.fasta     # Tier 1 genes only
├── 3UTR_antisense_tier1.fasta # Tier 1 reverse complement
└── isoform_map.json           # Transcript → Gene mapping
```

### Directory Sizes

| Gene Set | Size | Description |
|----------|------|-------------|
| `germ_plasm/` | 64 KB | 12 genes, includes shuffled controls |
| `housekeeping/` | 28 KB | 10 genes |
| `somatic/` | 120 KB | 8 genes |
| `cleared/` | 36 KB | 10 genes |
| `adult/` | 20 KB | 9 genes |

### Special Files

| File | Description |
|------|-------------|
| `germ_plasm/3UTR_shuffled.fasta` | Dinucleotide-shuffled controls |
| `germ_plasm/3UTR_shuffled_tier1.fasta` | Tier 1 shuffled controls |
| `genome_wide_sample_500.fasta` | Random sample for testing |

---

## TE Clusters (`data/te_clusters/`)

Output from `detect_te_clusters.py` - identified regions of high TE signal.

| File | Description |
|------|-------------|
| `clusters_summary.tsv` | Summary table of all detected clusters |
| `clusters_detail.json` | Detailed cluster information (coordinates, scores) |
| `cluster_sequences.fasta` | Extracted sequences from cluster regions |

---

## Cache Files (`data/cache/`)

Cached API responses (not critical, can be regenerated).

| File | Description |
|------|-------------|
| `datasets_api_eukaryota.parquet` | NCBI Datasets API cache |
| `ensembl_vertebrates_species.json` | Ensembl species list cache |

---

# Visualization Files

## HTML Visualizations (`results/`)

Interactive HTML files for exploring TE annotations.

### Germ Plasm Gene Annotations (`results/te_annotations/`)

**Index:** `results/te_annotations/index.html`

| File | Gene | UTR Length | Description |
|------|------|------------|-------------|
| `nos_3UTR_TE_annotation.html` | nanos | 880 bp | Posterior patterning |
| `osk_3UTR_TE_annotation.html` | oskar | 1,028 bp | Pole plasm organizer |
| `piwi_3UTR_TE_annotation.html` | piwi | 453 bp | piRNA pathway |
| `tud_3UTR_TE_annotation.html` | tudor | 787 bp | Polar granule component |
| `vas_3UTR_TE_annotation.html` | vasa | 104 bp | RNA helicase |
| `gcl_3UTR_TE_annotation.html` | germ cell-less | 524 bp | Germ cell formation |
| `aub_3UTR_TE_annotation.html` | aubergine | 113 bp | piRNA pathway |
| `pgc_3UTR_TE_annotation.html` | polar granule component | 401 bp | Germ cell marker |
| `CycB_3UTR_TE_annotation.html` | Cyclin B | varies | Cell cycle |
| `Kr_3UTR_TE_annotation.html` | Krüppel | ~800 bp | Gap gene (control) |

**TE-centric views:** `results/te_annotations/tes/`
- Shows which UTRs match each TE
- 12 TE-specific HTML files

### Genome-Wide Top/Bottom Genes (`results/te_annotations_genomewide/`)

**Index:** `results/te_annotations_genomewide/index.html`

#### Top 10 Highest TE Density Genes

| File | Rank | Gene | Symbol | Density |
|------|------|------|--------|---------|
| `top_01_FBgn0040959_te_annotation.html` | #1 | FBgn0040959 | Prt-15a | 325,271 |
| `top_02_FBgn0034403_te_annotation.html` | #2 | FBgn0034403 | CG18190 | 254,580 |
| `top_03_FBgn0067905_te_annotation.html` | #3 | FBgn0067905 | Dso2 | 253,347 |
| `top_04_FBgn0033948_te_annotation.html` | #4 | FBgn0033948 | CG12863 | 212,571 |
| `top_05_FBgn0053093_te_annotation.html` | #5 | FBgn0053093 | CG33093 | 210,638 |
| `top_06_FBgn0037514_te_annotation.html` | #6 | FBgn0037514 | CG10919 | 194,651 |
| `top_07_FBgn0260995_te_annotation.html` | #7 | FBgn0260995 | dpr21 | 171,804 |
| `top_08_FBgn0040534_te_annotation.html` | #8 | FBgn0040534 | Sf3b5 | 168,743 |
| `top_09_FBgn0053458_te_annotation.html` | #9 | FBgn0053458 | CG33458 | 165,500 |
| `top_10_FBgn0054003_te_annotation.html` | #10 | FBgn0054003 | NimB3 | 165,478 |

#### Bottom 10 Lowest TE Density Genes

| File | Rank | Gene |
|------|------|------|
| `bottom_01_FBgn0039833_te_annotation.html` | #13,289 | FBgn0039833 |
| `bottom_02_FBgn0038051_te_annotation.html` | #13,290 | FBgn0038051 |
| ... | ... | ... |
| `bottom_10_FBgn0034543_te_annotation.html` | #13,298 | FBgn0034543 |

### Full Transcript Analysis (`results/full_transcript_te/`)

**Index:** `results/full_transcript_te/index.html`

Shows TE hits across complete transcripts (5'UTR + CDS + 3'UTR) for top 10 genes.

| File | Gene | Key Finding |
|------|------|-------------|
| `Prt-15a_FBgn0040959_full_transcript.html` | Prt-15a | Chitin-binding protein |
| `CG18190_FBgn0034403_full_transcript.html` | CG18190 | Microtubule-associated |
| `Dso2_FBgn0067905_full_transcript.html` | Dso2 | **81% spanning** - likely domesticated TE |
| `CG12863_FBgn0033948_full_transcript.html` | CG12863 | Zinc finger protein |
| `CG33093_FBgn0053093_full_transcript.html` | CG33093 | **96% spanning** - likely TE-derived |
| `CG10919_FBgn0037514_full_transcript.html` | CG10919 | Testis-enriched |
| `dpr21_FBgn0260995_full_transcript.html` | dpr21 | **84% spanning** - cell adhesion |
| `Sf3b5_FBgn0040534_full_transcript.html` | Sf3b5 | **88% spanning** - splicing factor |
| `CG33458_FBgn0053458_full_transcript.html` | CG33458 | Serine protease |
| `NimB3_FBgn0054003_full_transcript.html` | NimB3 | Innate immunity |

**Key insight:** Many top genes have TE similarity spanning the CDS, suggesting they may be domesticated TEs rather than genes with TE insertions.

---

## Static Figures (`figures/`)

Publication-quality figures generated by visualization scripts.

### TE Signal Plots (`figures/te_signal/`)

Multi-panel figures showing TE signal density for germ plasm gene isoforms.

**Format:** `{gene}_{transcript}_te_signal.{png,pdf}`

| Gene | Isoforms | Files |
|------|----------|-------|
| nos | FBtr0083732, FBtr0335019 | 4 files |
| osk | FBtr0081954, FBtr0081956 | 4 files |
| vas | FBtr0445185, FBtr0445186, FBtr0445187 | 6 files |
| aub | FBtr0080165, FBtr0112793 | 4 files |
| gcl | FBtr0088710, FBtr0339337 | 4 files |
| pgc | FBtr0112519, FBtr0112520, FBtr0342987 | 6 files |
| CycB | FBtr0071911, FBtr0071913, FBtr0071914, FBtr0309858 | 8 files |

**Total:** 36 files (18 PNG + 18 PDF)

**Figure panels:**
1. Signal density line plot (Gaussian smoothed)
2. Nucleotide-resolution hit map (colored by TE family)
3. Windowed density heatmap (50bp bins)
4. GC content track

### Parameter Sweep Figures (`figures/parameter_sweep/`)

Visualizations of BLAST parameter optimization.

| File | Description |
|------|-------------|
| `parameter_heatmap.png` | Hit counts across parameter combinations |
| `task_comparison.png` | blastn vs blastn-short comparison |
| `reward_penalty_comparison.png` | Match/mismatch scoring effects |
| `top_combinations.png` | Best performing parameter sets |
| `evalue_distribution.png` | E-value distribution violin plots |

### Gene Comparison Figures (`figures/gene_comparison/`)

Cross-gene comparison visualizations.

| File | Description |
|------|-------------|
| `gene_summary.tsv` | Tabular summary of all genes |
| `gene_ranking.png` | Bar chart of TE density by gene |
| `te_family_breakdown.png` | Stacked bar chart of TE families |
| `length_vs_signal.png` | Scatter plot: UTR length vs TE signal |

### RepeatMasker Comparison Figures (`figures/repeatmasker_comparison/`)

Conservation, synteny, and quality analysis figures.

#### Conservation & Synteny Validation

| File | Description |
|------|-------------|
| `15_all_vs_hq_conservation.png` | Conservation comparison: all hits vs high-quality hits |
| `16_conservation_by_quality_tier.png` | Conservation stratified by identity and length tiers |
| `17_synteny_analysis.png` | Synteny validation across Drosophila species |

#### Ancient TE Candidates

**⚠️ ARCHIVED - Data Quality Issues Identified**

The following figures were based on v1 ancient candidate data with systematic false positives
(CDS overlap, multi-family hits, AT-rich spurious matches). Archived to:
`results/archive/ancient_te_candidates_v1_FLAWED/figures/`

| File | Description | Status |
|------|-------------|--------|
| `19_ancient_candidates_overview.png` | Overview: TE families and synteny distribution | **ARCHIVED** |
| `20_ancient_candidates_details.png` | Conservation and length distributions | **ARCHIVED** |

See `results/archive/ancient_te_candidates_v1_FLAWED/DATA_QUALITY_ISSUES.md` for full details.

New v2 analysis available at: `results/repeatmasker_analysis_v2/`

#### Strand Bias Analysis

| File | Description |
|------|-------------|
| `21_strand_vs_conservation.png` | Four-panel analysis showing no strand-conservation correlation |

#### Quality Paradox (Publication Figures)

| File | Description |
|------|-------------|
| `22_quality_vs_conservation_synteny.png` | Scatter plot with marginal histograms - two populations |
| `23_quality_conservation_density.png` | Side-by-side density: syntenic vs non-syntenic |
| `24_quality_metrics_comparison.png` | Multiple quality metrics compared |
| `25_identical_bp_vs_conservation.png` | **Publication figure** - clean two-population visualization |

**Key insight from Figure 25**: Two distinct populations visible when plotting BLAST quality (identical bases) against conservation (phyloP), colored by synteny:
- **Ancient fossils**: Low quality (~41 bp), high conservation (phyloP 1.43), syntenic
- **Recent insertions**: High quality (~99 bp), low conservation (phyloP 0.47), non-syntenic

#### Shuffled Control Analysis

Validation using dinucleotide-shuffled sequences (10 replicates, 10% sample).

| File | Description |
|------|-------------|
| `26_shuffled_control_comparison.png` | Bar charts: real vs shuffled (2.2x total, 92x HQ enrichment) |
| `27_shuffled_replicate_distributions.png` | Distribution of hits across 10 shuffled replicates |
| `28_te_family_enrichment.png` | TE family enrichment analysis (all families enriched in real) |
| `29_real_vs_shuffled_distributions.png` | Identity, length, E-value distribution comparisons |

**Key findings**:
- **2.2x enrichment** for total hits (Z=104.5)
- **92x enrichment** for high-quality hits (≥80%, ≥50bp)
- All TE families enriched in real sequences (1.7-14.8x)
- E-value < 10⁻⁵ essentially exclusive to real sequences

---

## Intermediate/Bad Visualizations

These files are marked as intermediate or contain incorrect data:

| File | Status | Issue |
|------|--------|-------|
| `INTERMEDIATE_Kr_3UTR_TE_annotation.html` | Intermediate | Early exploratory version |
| `INTERMEDIATE_Kr_3UTR_TE_annotation_v2.html` | Intermediate | Early exploratory version |

---

## How to View HTML Visualizations

### Option 1: Open directly in browser

```bash
open results/te_annotations/index.html
open results/te_annotations_genomewide/index.html
open results/full_transcript_te/index.html
```

### Option 2: Start local server

```bash
cd results
python -m http.server 8000
# Then open http://localhost:8000 in browser
```

### Option 3: VS Code Live Server

If using VS Code, right-click on any HTML file and select "Open with Live Server".

---

## File Size Summary

| Category | Total Size | Files |
|----------|------------|-------|
| References | ~6.6 GB | 11 files |
| BLAST databases | ~50 MB | 40+ files |
| Gene lists | ~60 KB | 14 files |
| Query sequences | ~270 KB | 25+ files |
| HTML visualizations | ~2.5 MB | 60+ files |
| Static figures | ~15 MB | 45+ files |

**Note:** The largest file is `dmel_annotation.gff` at 6.3 GB. This can be deleted if not needed for other analyses.

---

## Data Provenance

| Data Type | Source | Version/Release |
|-----------|--------|-----------------|
| Genome & annotations | FlyBase | r6.66 |
| 3'UTR sequences | FlyBase | r6.66 |
| TE insertions | FlyBase | r6.66 |
| TE consensus | Bergman Lab | drosophila-transposons |
| Gene groups | FlyBase | r6.66 |

---

*Document created: 2026-01-14*
