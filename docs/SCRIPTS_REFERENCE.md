# Scripts Reference Guide

This document provides a complete audit of all Python scripts in the `repeat_finder` project, organized by their role in the analysis pipeline.

> **Last updated**: 2026-01-14
> **Note**: Scripts renamed for clarity on 2026-01-14. Old names in parentheses.

## Quick Reference

| Category | Scripts | Purpose |
|----------|---------|---------|
| **Setup** | 4 scripts | Download data, build databases |
| **Gene Lists** | 4 scripts | Build gene lists for analysis |
| **Query Prep** | 3 scripts | Extract and prepare sequences |
| **BLAST** | 3 scripts | Run BLAST searches |
| **Analysis** | 7 scripts | Analyze BLAST results |
| **Visualization** | 4 scripts | Generate figures and HTML |
| **Utilities** | 4 scripts | Shared modules, list info, summarize |

## Shared Utilities (`scripts/utils/`)

New shared modules added to reduce code duplication:

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `paths.py` | Centralized path resolution | `get_project_root()`, `get_results_dir()`, `get_gene_lists_dir()` |
| `blast_io.py` | BLAST parsing | `BLAST_COLUMNS`, `load_blast_results()`, `classify_strand()` |
| `data_loaders.py` | Data loading | `load_gene_list()`, `parse_fasta()`, `load_te_database()` |

**Usage:**
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root
from utils.blast_io import load_blast_results
```

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           SETUP PHASE                                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   download_references.py  ──►  build_databases.py                           │
│         │                            │                                       │
│         ▼                            ▼                                       │
│   data/references/            data/blastdb/                                 │
│   ├── dmel_3utr.fasta        ├── dmel_te_flybase.*                         │
│   ├── dmel_5utr.fasta        └── dmel_te_consensus.*                       │
│   ├── dmel_cds.fasta                                                        │
│   └── dmel_te_flybase.fasta                                                 │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         GENE LIST PHASE                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   build_germ_plasm_genelist.py  ──►  data/gene_lists/germ_plasm_*.tsv      │
│   build_housekeeping_genelist.py ──► data/gene_lists/housekeeping_*.tsv    │
│   build_control_genelists.py    ──►  data/gene_lists/{somatic,cleared,     │
│                                                        adult}_*.tsv         │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        QUERY PREPARATION                                     │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   extract_germ_plasm_3utrs.py  ──►  data/queries/germ_plasm/               │
│         │                           ├── 3UTR_sense.fasta                    │
│         │                           ├── 3UTR_antisense.fasta                │
│         │                           └── 3UTR_sense_tier1.fasta              │
│         │                                                                    │
│         └──► shuffle_sequences.py  ──►  3UTR_shuffled.fasta                │
│                                                                              │
│   filter_queries.py  (optional: filter by FlyBase gene groups)              │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          BLAST PHASE                                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   te_parameter_sweep.py  ──►  results/parameter_sweep/                      │
│         │                     ├── sweep_summary.tsv                         │
│         │                     └── combo_*/blast_results.tsv                 │
│         │                                                                    │
│         │   (or use blast_runner.py for individual runs)                    │
│         │                                                                    │
│         └──► Best parameters: word_size=7, gapopen=2, gapextend=1,         │
│                               dust=yes, evalue=10                           │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         ANALYSIS PHASE                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   analyze_genome_wide_te.py  ──►  results/                                  │
│         │                         ├── top_100_te_genes_FIXED.tsv            │
│         │                         ├── strand_bias_by_utr.tsv                │
│         │                         └── te_annotations_genomewide/            │
│         │                                                                    │
│   analyze_full_transcripts.py ──► results/full_transcript_te/               │
│         │                                                                    │
│   rank_gene_sets.py  ──►  (console output: gene set comparisons)         │
│                                                                              │
│   calculate_te_signal_density.py  ──►  density pickle + peaks               │
│         │                                                                    │
│         └──► detect_te_clusters.py  ──►  cluster summaries                  │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       VISUALIZATION PHASE                                    │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   plot_te_signal.py          ──►  figures/te_signal/*.png                   │
│   plot_parameter_sweep.py    ──►  figures/parameter_sweep/*.png             │
│   visualize_gene_comparison.py           ──►  figures/*.png (rankings, breakdowns)      │
│   generate_te_fossil_report.py ──► reports/te_fossil_analysis.html          │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Script Documentation

### SETUP SCRIPTS

#### `download_references.py`
Downloads reference files from FlyBase and Dfam.

```bash
# Download all reference files
python scripts/download_references.py --output-dir data/references

# List available files
python scripts/download_references.py --list

# Download specific files
python scripts/download_references.py --files genome 3utr te_flybase
```

**Outputs:**
- `dmel_3utr.fasta` - 3'UTR sequences (30,324 sequences)
- `dmel_5utr.fasta` - 5'UTR sequences
- `dmel_cds.fasta` - Coding sequences
- `dmel_te_flybase.fasta` - FlyBase TE sequences (5,734 TEs)
- `gene_groups.tsv.gz` - FlyBase gene group memberships

---

#### `build_databases.py`
Builds BLAST nucleotide databases from reference FASTA files.

```bash
# Build all databases
python scripts/build_databases.py --reference-dir data/references --output-dir data/blastdb

# Build specific database
python scripts/build_databases.py --databases te_flybase

# Force rebuild
python scripts/build_databases.py --force
```

**Outputs:**
- `data/blastdb/dmel_te_flybase.{nhr,nin,nsq}` - TE BLAST database

---

### GENE LIST SCRIPTS

#### `build_germ_plasm_genelist.py`
Builds the canonical germ plasm gene list.

```bash
python scripts/build_germ_plasm_genelist.py --output-dir data/gene_lists
```

**Genes included:**
- **Tier 1** (core germ plasm): nos, osk, pgc, gcl, vas, aub, CycB
- **Tier 2** (associated): tudor, piwi, AGO3, dhd, Hsp83

**Outputs:**
- `germ_plasm_genes_consolidated.tsv` - Full gene list with FBgn IDs
- `germ_plasm_genes_tier1.txt` - Tier 1 genes only
- `germ_plasm_fbgn_ids.txt` - FBgn IDs only

---

#### `build_housekeeping_genelist.py`
Builds housekeeping gene list for controls.

```bash
python scripts/build_housekeeping_genelist.py --output-dir data/gene_lists
```

**Genes:** Act5C, RpL32, Gapdh1, alphaTub84B, Ef1alpha48D, RpS17, Tbp, eIF4A, SdhA, Rpl4

---

#### `build_control_genelists.py`
Builds control gene lists (somatic, cleared, adult).

```bash
python scripts/build_control_genelists.py --output-dir data/gene_lists

# List available gene sets
python scripts/build_control_genelists.py --list
```

---

### QUERY PREPARATION SCRIPTS

#### `extract_germ_plasm_3utrs.py`
Extracts 3'UTRs for genes in a gene list.

```bash
python scripts/extract_germ_plasm_3utrs.py \
    --reference data/references/dmel_3utr.fasta \
    --gene-list data/gene_lists/germ_plasm_genes_consolidated.tsv \
    --output-dir data/queries/germ_plasm
```

**Outputs:**
- `3UTR_sense.fasta` - Sense strand sequences
- `3UTR_antisense.fasta` - Reverse complement (for antisense control)
- `3UTR_sense_tier1.fasta` - Tier 1 genes only
- `isoform_mapping.json` - Transcript-to-gene mapping

---

#### `shuffle_sequences.py`
Creates dinucleotide-shuffled control sequences (Altschul-Erickson algorithm).

```bash
python scripts/shuffle_sequences.py \
    --input data/queries/germ_plasm/3UTR_sense.fasta \
    --output data/queries/germ_plasm/3UTR_shuffled.fasta \
    --seed 42
```

**Key feature:** Preserves dinucleotide frequencies, providing biologically realistic null sequences.

---

#### `filter_queries.py`
Filters FASTA by FlyBase gene group membership.

```bash
# List available gene groups
python scripts/filter_queries.py --list-groups

# Filter by specific group
python scripts/filter_queries.py \
    --gene-group "GERM_CELL" \
    --segment-types 3utr \
    --output-dir data/queries/germ_cell
```

---

### BLAST SCRIPTS

#### `te_parameter_sweep.py` ⭐ **PRIMARY BLAST SCRIPT**
Systematic parameter optimization for TE detection.

```bash
# Run full sweep (324 combinations)
python scripts/te_parameter_sweep.py \
    --query data/queries/germ_plasm/3UTR_sense.fasta \
    --database data/blastdb/dmel_te_flybase \
    --output-dir results/parameter_sweep \
    --threads 8

# Run on tier1 genes only
python scripts/te_parameter_sweep.py --tier tier1

# List parameter combinations
python scripts/te_parameter_sweep.py --list-params
```

**Parameter grid:**
- `word_size`: [7, 9, 11]
- `gapopen`: [2, 5, 10]
- `gapextend`: [1, 2, 4]
- `reward/penalty`: [(1,-1), (1,-2), (1,-3), (2,-3)]
- `task`: ['blastn', 'blastn-short']
- **Always:** `dust=yes`, `evalue=10`

**Best parameters found:** `word_size=7, gapopen=2, gapextend=1, penalty=-1, reward=1, dust=yes`

---

#### `blast_runner.py`
Flexible BLAST runner with YAML configuration.

```bash
python scripts/blast_runner.py \
    --config configs/blast_config.yaml \
    --database dmel_te_flybase \
    --output results/my_blast.tsv
```

**Use case:** When you need custom BLAST parameters or want config file management.

---

#### `analyze_te_regions.py`
Maps BLAST hits to TE structural regions (LTR, CDS, etc.).

```bash
python scripts/analyze_te_regions.py \
    --query data/queries/germ_plasm/3UTR_sense.fasta \
    --db data/blastdb/dmel_te_consensus \
    --gff data/references/te_annotations.gff \
    --output results/te_regions.tsv
```

---

### ANALYSIS SCRIPTS

#### `analyze_genome_wide_te.py` ⭐ **PRIMARY ANALYSIS SCRIPT**
Genome-wide TE content analysis with proper gene-level aggregation.

```bash
python scripts/analyze_genome_wide_te.py \
    --blast results/genome_wide_all_3utrs.tsv \
    --utr-fasta data/references/dmel_3utr.fasta \
    --te-fasta data/references/dmel_te_flybase.fasta \
    --output results/ \
    --html-top 10 \
    --html-bottom 10
```

**Outputs:**
- `top_100_te_genes_FIXED.tsv` - Top 100 genes by TE density
- `bottom_100_te_genes_FIXED.tsv` - Bottom 100 genes
- `strand_bias_by_utr.tsv` - Per-UTR strand statistics
- `strand_bias_by_te.tsv` - Per-TE strand statistics
- `te_annotations_genomewide/*.html` - HTML visualizations

**Key features:**
- Proper gene-level aggregation (sums UTR lengths across isoforms)
- Strand bias analysis (sense vs antisense)
- HTML visualization generation

---

#### `analyze_full_transcripts.py`
Analyzes TE content across complete transcripts (5'UTR + CDS + 3'UTR).

```bash
python scripts/analyze_full_transcripts.py
```

**Purpose:** Identifies whether high TE density extends beyond 3'UTR into CDS (possible domesticated TEs).

**Outputs:** HTML visualizations showing TE hits across full transcript structure.

---

#### `rank_gene_sets.py` *(renamed from compare_gene_sets.py)*
Ranks predefined gene sets by genome-wide TE metrics.

```bash
python scripts/rank_gene_sets.py
```

**Outputs:** Console tables showing:
- Per-gene TE density rankings
- Percentile positions
- Strand bias by gene
- Cross-group comparisons

---

#### `calculate_te_signal_density.py`
Calculates position-wise TE signal density with Gaussian smoothing.

```bash
python scripts/calculate_te_signal_density.py \
    results/blast_results.tsv \
    --output-dir results/density \
    --sigma 25 \
    --query-fasta data/queries/germ_plasm/3UTR_sense.fasta
```

**Outputs:**
- `density_data.pkl` - Pickle file with density arrays
- `density_summary.tsv` - Summary statistics
- `peaks.json` - Peak coordinates and annotations

---

#### `detect_te_clusters.py`
Identifies clusters of TE signal from density data.

```bash
python scripts/detect_te_clusters.py \
    results/density \
    --query-fasta data/queries/germ_plasm/3UTR_sense.fasta \
    --output-dir results/clusters \
    --threshold 2.0 \
    --merge-distance 50
```

---

#### `analyze_te_families.py`
Analyzes TE family enrichment across gene sets.

```bash
python scripts/analyze_te_families.py \
    --germ-plasm results/germ_plasm.tsv \
    --housekeeping results/housekeeping.tsv \
    --output-dir results/te_families
```

---

### CONSERVATION & VALIDATION SCRIPTS (`scripts/conservation_analysis/` and `scripts/`)

Scripts for validating TE hits using cross-species conservation and synteny data.

#### `scripts/conservation_analysis/analyze_conservation.py`
Extracts phyloP conservation scores for TE hit regions.

```bash
python scripts/conservation_analysis/analyze_conservation.py
```

**Inputs:**
- UCSC dm6.phyloP27way.bw (download required)
- TE hit BED files

**Outputs:**
- `te_hits_all_conservation.tab` - Conservation scores for all hits

---

#### `scripts/conservation_analysis/fast_synteny_analysis.py`
Fast synteny analysis using MAF alignments with binary search optimization.

```bash
python scripts/conservation_analysis/fast_synteny_analysis.py
```

**Inputs:**
- UCSC dm6.multiz27way.maf files
- TE hit coordinates

**Outputs:**
- `te_hits_all_synteny_sampled.tsv` - Synteny scores across species
- `te_hits_hq_synteny.tsv` - High-quality subset

---

#### `scripts/filter_ancient_te_candidates.py` (DEPRECATED - use v2)

**⚠️ DEPRECATED**: This script has systematic false positives. Use `filter_ancient_te_candidates_v2.py` instead.

Issues identified in v1 output:
- Top candidates overlap CDS in alternative isoforms (coding-level phyloP)
- Multi-family hits to same region (AT-rich spurious matches)
- HeT-A (telomeric) hits are false positives

---

#### `scripts/filter_ancient_te_candidates_v2.py` (SUPERSEDED)

**⚠️ SUPERSEDED**: Filtering approach was too aggressive. Use `build_utr_te_loci.py` instead.

This script filtered out CDS overlaps and multi-family hits, but this was overly aggressive:
- CDS overlap doesn't invalidate UTR function (dual-function sequences exist)
- Multi-family hits aren't necessarily spurious (related TEs share sequences)

See `results/utr_te_loci/METHODOLOGY_NOTES.md` for detailed rationale.

---

#### `scripts/build_utr_te_loci.py` (CURRENT)

Builds TE loci per UTR isoform with annotations (no filtering).

```bash
python scripts/build_utr_te_loci.py
```

**Approach**: Annotate, don't filter
- One UTR model per transcript isoform
- CDS overlap **annotated** (`cds_overlap: yes/no`) not filtered
- Multi-family hits **annotated** (`n_te_families`) not filtered
- User decides what to filter based on their biological question

**Outputs** (to `results/utr_te_loci/`):
- `ancient_te_loci_by_isoform.tsv` - 25,240 loci across 14,423 UTR isoforms
- `METHODOLOGY_NOTES.md` - Rationale for annotate-not-filter approach

**Naming convention**: `{transcript}|{chrom}|{start}-{end}`

---

#### `scripts/analyze_strand_conservation.py`
Analyzes correlation between strand bias and conservation.

```bash
python scripts/analyze_strand_conservation.py
```

**Key finding:** No meaningful correlation (r = 0.016-0.029)

---

#### `scripts/plot_strand_conservation.py`
Visualizes strand bias vs conservation relationships.

**Output:** `figures/repeatmasker_comparison/21_strand_vs_conservation.png`

---

#### `scripts/plot_quality_vs_conservation.py`
Creates the "killer figure" showing two distinct populations.

```bash
python scripts/plot_quality_vs_conservation.py
```

**Key output:** `figures/repeatmasker_comparison/25_identical_bp_vs_conservation.png`

Shows:
- **Ancient fossils**: Low quality (~41 bp), high conservation (phyloP 1.43), syntenic
- **Recent insertions**: High quality (~99 bp), low conservation (phyloP 0.47), non-syntenic

---

#### `scripts/plot_quality_vs_conservation_v2.py`
Extended version with multiple quality metrics compared.

**Outputs:**
- `24_quality_metrics_comparison.png` - Four-panel comparison
- `25_identical_bp_vs_conservation.png` - Clean publication figure

---

#### `scripts/visualize_ancient_candidates.py`
Visualizes the ancient TE candidate dataset.

**Outputs:**
- `19_ancient_candidates_overview.png`
- `20_ancient_candidates_details.png`

---

### SHUFFLED CONTROL SCRIPTS

Scripts for validating TE detection against dinucleotide-shuffled control sequences.

#### `scripts/run_shuffled_control_analysis.py`
Main shuffled control analysis pipeline.

```bash
# Quick analysis (10% sample, 10 replicates)
python scripts/run_shuffled_control_analysis.py --sample-frac 0.1 --n-shuffles 10

# Full analysis (all UTRs, 10 replicates)
python scripts/run_shuffled_control_analysis.py --full --n-shuffles 10
```

**Key findings:**
- 2.2x enrichment for total hits
- 92x enrichment for high-quality hits (≥80%, ≥50bp)

**Outputs** (to `results/shuffled_controls/`):
- `control_comparison.tsv` - Summary statistics
- `shuffled_replicates.tsv` - Per-replicate data

---

#### `scripts/analyze_shuffled_te_families.py`
Analyzes TE family enrichment between real and shuffled sequences.

**Outputs:**
- `te_family_enrichment.tsv` - Family-level enrichment data
- `28_te_family_enrichment.png` - TE family enrichment visualization
- `29_real_vs_shuffled_distributions.png` - Distribution comparisons

---

#### `scripts/plot_shuffled_control_results.py`
Visualizes shuffled control results.

**Outputs:**
- `26_shuffled_control_comparison.png`
- `27_shuffled_replicate_distributions.png`

---

### VISUALIZATION SCRIPTS

#### `plot_te_signal.py`
Publication-quality TE signal plots.

```bash
# Plot all sequences
python scripts/plot_te_signal.py \
    results/density \
    --query-fasta data/queries/germ_plasm/3UTR_sense.fasta \
    --output-dir figures/te_signal \
    --all

# Plot specific sequence
python scripts/plot_te_signal.py ... --sequence FBtr0073419
```

**Outputs:** Multi-panel figures with:
- Signal density line plot
- Nucleotide-resolution hit map
- Windowed density heatmap
- GC content track

---

#### `plot_parameter_sweep.py`
Visualizes parameter sweep results.

```bash
python scripts/plot_parameter_sweep.py \
    results/parameter_sweep \
    --output-dir figures/parameter_sweep
```

**Outputs:** Heatmaps showing hit counts across parameter combinations.

---

#### `visualize_gene_comparison.py` *(renamed from compare_genes.py)*
Publication-quality cross-gene comparison visualizations.

```bash
python scripts/visualize_gene_comparison.py \
    results/density \
    --output-dir figures
```

**Outputs:**
- Ranking bar charts
- TE family breakdown stacked bars
- Length vs signal scatter plots

---

#### `analyze_utr_variants.py` *(renamed from compare_5utr_3utr.py)*
Compares TE content between 5'UTR and 3'UTR regions.

```bash
python scripts/analyze_utr_variants.py
```

**Outputs:** Console tables comparing TE density, strand bias, and hit patterns between UTR types.

---

#### `summarize_all_gene_sets.py` *(renamed from compare_all_gene_sets.py)*
Comprehensive summary of 5'UTR and 3'UTR TE content across all gene sets.

```bash
python scripts/summarize_all_gene_sets.py
```

**Outputs:** Detailed tables showing TE hit counts, strand bias, and UTR lengths for each gene.

---

#### `generate_te_fossil_report.py`
Comprehensive HTML/Markdown report.

```bash
python scripts/generate_te_fossil_report.py \
    --density-dir results/density \
    --clusters-file results/clusters/cluster_summary.tsv \
    --gene-list data/gene_lists/germ_plasm_genes_consolidated.tsv \
    --output-dir reports
```

---

### UTILITY SCRIPTS

#### `list_databases.py`
Lists available BLAST databases.

```bash
python scripts/list_databases.py --db-dir data/blastdb --verbose
```

---

#### `list_gene_groups.py`
Lists FlyBase gene groups.

```bash
python scripts/list_gene_groups.py --search "germ" --min-genes 5
```

---

#### `summarize_results.py`
Quick summary of BLAST results.

```bash
python scripts/summarize_results.py results/blast.tsv --plot --top 20
```

---

#### `extract_top_hits.py`
Display top BLAST alignments with filters.

```bash
python scripts/extract_top_hits.py results/blast.tsv \
    --te-fasta data/references/dmel_te_flybase.fasta \
    --top 30 \
    --min-length 100 \
    --gene nos
```

---

## Typical Workflows

### Full Pipeline (New Analysis)

```bash
# 1. Setup
python scripts/download_references.py --output-dir data/references
python scripts/build_databases.py --reference-dir data/references --output-dir data/blastdb

# 2. Gene lists
python scripts/build_germ_plasm_genelist.py --output-dir data/gene_lists
python scripts/build_housekeeping_genelist.py --output-dir data/gene_lists
python scripts/build_control_genelists.py --output-dir data/gene_lists

# 3. Extract queries
python scripts/extract_germ_plasm_3utrs.py \
    --reference data/references/dmel_3utr.fasta \
    --gene-list data/gene_lists/germ_plasm_genes_consolidated.tsv \
    --output-dir data/queries/germ_plasm

# 4. Create shuffled controls
python scripts/shuffle_sequences.py \
    --input data/queries/germ_plasm/3UTR_sense.fasta \
    --output data/queries/germ_plasm/3UTR_shuffled.fasta

# 5. Run BLAST (genome-wide)
# Use blast_runner.py with optimized parameters on all 3'UTRs

# 6. Analyze
python scripts/analyze_genome_wide_te.py \
    --blast results/genome_wide_all_3utrs.tsv \
    --utr-fasta data/references/dmel_3utr.fasta \
    --te-fasta data/references/dmel_te_flybase.fasta \
    --output results/

# 7. Compare gene sets
python scripts/rank_gene_sets.py
```

### Quick Analysis (Existing Data)

```bash
# Summarize existing BLAST results
python scripts/summarize_results.py results/blast.tsv --plot

# View top hits
python scripts/extract_top_hits.py results/blast.tsv --top 20 --min-length 100

# Generate report
python scripts/generate_te_fossil_report.py --density-dir results/density --output-dir reports
```

---

## Script Categories Summary

### Core Pipeline Scripts (run in order)
1. `download_references.py`
2. `build_databases.py`
3. `build_*_genelist.py`
4. `extract_germ_plasm_3utrs.py`
5. `te_parameter_sweep.py` or `blast_runner.py`
6. `analyze_genome_wide_te.py`

### Conservation & Validation Scripts (run after BLAST)
- `conservation_analysis/analyze_conservation.py` - Extract phyloP scores
- `conservation_analysis/fast_synteny_analysis.py` - Cross-species synteny
- `build_utr_te_loci.py` - Build annotated TE loci per UTR isoform **(CURRENT)**
- `analyze_strand_conservation.py` - Strand bias analysis
- `plot_quality_vs_conservation.py` - **Publication figure**
- `plot_quality_vs_conservation_v2.py` - Extended metrics comparison

### One-Off Analysis Scripts
- `analyze_full_transcripts.py` - Check if TE extends into CDS
- `analyze_te_regions.py` - Map to TE structural regions
- `analyze_te_families.py` - TE family enrichment
- `rank_gene_sets.py` - Compare gene sets

### Shuffled Control Scripts (validation)
- `run_shuffled_control_analysis.py` - Main shuffled control pipeline
- `analyze_shuffled_te_families.py` - TE family enrichment analysis
- `plot_shuffled_control_results.py` - Control comparison visualization

### Visualization Scripts (run after analysis)
- `plot_te_signal.py`
- `plot_parameter_sweep.py`
- `visualize_gene_comparison.py`
- `generate_te_fossil_report.py`
- `plot_strand_conservation.py` - Strand bias visualizations
- `visualize_ancient_candidates.py` - Ancient candidate figures

### Utility Scripts (run anytime)
- `list_databases.py`
- `list_gene_groups.py`
- `summarize_results.py`
- `extract_top_hits.py`
- `filter_queries.py`

---

*Document created: 2026-01-14*
