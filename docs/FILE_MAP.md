# Repository File Map

> **Last updated**: 2026-01-13
> **Purpose**: Track what each file contains and how the repo has evolved
> **Maintenance**: Update this file when adding new scripts, data, or results

---

## Quick Reference

| What you want | Where to find it |
|---------------|------------------|
| Latest germ plasm TE hits | `results/diverged_controls/germ_plasm_sense.tsv` |
| Analysis summary | `docs/SESSION_SUMMARY_diverged_TE_analysis.md` |
| BLAST parameters | `results/DIVERGED_TE_ANALYSIS_SUMMARY.md` |
| Gene lists | `data/gene_lists/` |
| 3'UTR sequences | `data/queries/{group}/3UTR_sense.fasta` |
| Run BLAST | `scripts/blast_runner.py` |

---

## Directory Structure

```
repeat_finder/
├── data/                    # Input data (tracked, some downloaded)
│   ├── gene_lists/          # Gene sets for different tissue groups
│   ├── queries/             # Extracted 3'UTR sequences for BLAST
│   ├── references/          # Reference genomes, TEs, annotations
│   ├── blastdb/             # BLAST databases (gitignored)
│   └── cache/               # API response caches
├── scripts/                 # Analysis scripts (tracked)
├── results/                 # BLAST outputs, analyses (gitignored)
├── reports/                 # HTML/MD reports (gitignored)
├── figures/                 # Generated plots (gitignored)
├── docs/                    # Documentation (tracked)
├── src/                     # Legacy/template code
└── workflows/               # Snakemake workflows
```

---

## Scripts

### Core Pipeline Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `download_references.py` | Download FlyBase/Dfam reference data | URLs | `data/references/` |
| `build_databases.py` | Build BLAST databases | FASTAs | `data/blastdb/` |
| `blast_runner.py` | Run BLAST searches | Query FASTA, DB | TSV results |

### Gene List Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `build_germ_plasm_genelist.py` | Build canonical germ plasm gene list | Hardcoded + FlyBase | `data/gene_lists/germ_plasm_*` |
| `build_housekeeping_genelist.py` | Build housekeeping control list | FlyBase API | `data/gene_lists/housekeeping_*` |
| `build_control_genelists.py` | Build somatic/cleared/adult lists | FlyBase API | `data/gene_lists/{group}_*` |
| `list_gene_groups.py` | Query available FlyBase gene groups | FlyBase API | stdout |

### 3'UTR Extraction Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `extract_germ_plasm_3utrs.py` | Extract 3'UTRs for gene lists | Gene list + reference FASTA | `data/queries/{group}/` |
| `shuffle_sequences.py` | Dinucleotide shuffle for controls | FASTA | `*_shuffled.fasta` |
| `filter_queries.py` | Filter queries by criteria | FASTA | Filtered FASTA |

### Analysis Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `te_parameter_sweep.py` | Test BLAST parameter combinations | Queries + DB | `results/parameter_sweep/` |
| `extract_top_hits.py` | **NEW** - Visualize top BLAST alignments | BLAST TSV | Formatted text |
| `analyze_te_families.py` | Analyze TE family distribution | BLAST TSV | Summary stats |
| `calculate_te_signal_density.py` | Position-wise TE signal | BLAST TSV | Density arrays |
| `detect_te_clusters.py` | Find TE signal hotspots | Density data | Cluster coords |
| `compare_genes.py` | Cross-gene comparison | Multiple TSVs | Summary table |
| `summarize_results.py` | Generate result summaries | BLAST TSV | Stats |

### Visualization Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `plot_te_signal.py` | Per-gene signal density plots | `figures/te_signal/` |
| `plot_parameter_sweep.py` | Parameter sweep heatmaps | `figures/parameter_sweep/` |

### Report Generation

| Script | Purpose | Output |
|--------|---------|--------|
| `generate_te_fossil_report.py` | Comprehensive HTML/MD report | `reports/` |

---

## Data Files

### Gene Lists (`data/gene_lists/`)

| File Pattern | Contents | Genes |
|--------------|----------|-------|
| `germ_plasm_genes_consolidated.tsv` | Canonical germ plasm genes | nos, osk, pgc, gcl, vas, aub, tud, piwi, AGO3, dhd, CycB, Hsp83 |
| `germ_plasm_genes_tier1.txt` | Priority subset | nos, osk, pgc, gcl, vas, aub, CycB |
| `housekeeping_genes_consolidated.tsv` | Negative control genes | RpL32, Act5C, Gapdh1, etc. |
| `somatic_genes_consolidated.tsv` | Somatically-localized genes | bcd, nos, osk, stau, etc. |
| `cleared_genes_consolidated.tsv` | Posteriorly-cleared genes | hb, Kr, kni, etc. |
| `adult_genes_consolidated.tsv` | Adult-expressed genes | Various |
| `*_fbgn_ids.txt` | FlyBase gene IDs only | FBgn IDs |
| `*_status.json` | Build metadata | Timestamps, counts |

### Query Sequences (`data/queries/{group}/`)

Each group (germ_plasm, housekeeping, somatic, cleared, adult) contains:

| File | Contents |
|------|----------|
| `3UTR_sense.fasta` | All 3'UTRs, sense strand |
| `3UTR_antisense.fasta` | Reverse complement (control) |
| `3UTR_sense_tier1.fasta` | Priority genes only |
| `3UTR_antisense_tier1.fasta` | Priority genes, antisense |
| `3UTR_shuffled.fasta` | Dinucleotide-shuffled (germ_plasm only) |
| `isoform_map.json` | Gene → transcript mapping |

### Reference Data (`data/references/`)

| File | Source | Contents |
|------|--------|----------|
| `dmel_3utr.fasta` | FlyBase | 30,324 3'UTR sequences |
| `dmel_5utr.fasta` | FlyBase | 5'UTR sequences |
| `dmel_cds.fasta` | FlyBase | Coding sequences |
| `dmel_te_flybase.fasta` | FlyBase | 5,734 TE instances |
| `dmel_genome.fasta` | FlyBase | Genome assembly |
| `dmel_annotation.gff` | FlyBase | Gene annotations |
| `gene_info.tsv` | FlyBase | FBgn/FBtr/FBpp mapping |
| `gene_groups.tsv` | FlyBase | Gene group definitions |

---

## Results (gitignored, local only)

### Current Best Results

| File | What it contains | Parameters |
|------|------------------|------------|
| `results/diverged_controls/germ_plasm_sense.tsv` | **CURRENT** - Germ plasm TE hits with DUST | word_size=7, gapopen=2, gapextend=1, dust=yes |
| `results/diverged_controls/germ_plasm_antisense.tsv` | Antisense control | Same |
| `results/diverged_controls/housekeeping_sense.tsv` | Negative control | Same |
| `results/diverged_controls/somatic_sense.tsv` | Somatic comparison | Same |
| `results/diverged_controls/cleared_sense.tsv` | Cleared comparison | Same |
| `results/diverged_controls/adult_sense.tsv` | Adult comparison | Same |
| `results/diverged_controls/shuffled.tsv` | Statistical control | Same |

### Superseded Results (historical)

| File/Directory | Status | Why superseded |
|----------------|--------|----------------|
| `results/parameter_sweep/` | OLD | Used dust=no, captured simple repeats |
| `results/dust_sweep/` | INTERMEDIATE | Parameter sweep with dust=yes |
| `results/controls/` | OLD | Pre-DUST control runs |
| `results/germ_plasm_dust_on.tsv` | OLD | Single test run |
| `results/antisense_best.tsv` | OLD | Pre-optimized parameters |
| `results/shuffled_best.tsv` | OLD | Pre-optimized parameters |

### Analysis Outputs

| File | Contents |
|------|----------|
| `results/DIVERGED_TE_ANALYSIS_SUMMARY.md` | Parameter comparison, top candidates |
| `results/DIVERGED_TE_CONTROL_COMPARISON.md` | Cross-group enrichment analysis |
| `results/te_fossil_candidates.txt` | Curated candidate list |
| `results/te_fossils_diverged.txt` | Top 50 diverged hits with alignments |
| `results/top_hits_*.txt` | Per-gene hit summaries |

---

## Reports (`reports/`)

Multiple timestamped versions exist. Latest is most current:

| Pattern | Contents |
|---------|----------|
| `te_fossil_analysis_{timestamp}.html` | Full HTML report |
| `te_fossil_analysis_{timestamp}.md` | Markdown version |

---

## Evolution / Version History

### Version 1: Simple Repeat Problem (superseded)
- **Parameters**: dust=no
- **Problem**: 90% of hits were AT-repeats, poly-A
- **Files**: `results/parameter_sweep/`, `results/controls/`

### Version 2: DUST-Filtered (current)
- **Parameters**: dust=yes, word_size=7, gapopen=2, gapextend=1
- **Improvement**: 99.9% complex sequences, 60% with gaps
- **Files**: `results/diverged_controls/`

---

## BLAST Output Format

All TSV files use 17-column format:

```
qseqid, sseqid, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, evalue, bitscore,
qlen, slen, qseq, sseq, strand
```

The `strand` column indicates TE orientation:
- `plus`: sstart < send (TE in sense orientation)
- `minus`: sstart > send (TE in antisense orientation)

---

## Maintenance Notes

When adding new files:
1. Update the relevant section in this document
2. If superseding old results, move old files to "Superseded" section
3. Document the parameters/methods used
4. Update the "Quick Reference" table if adding key outputs
