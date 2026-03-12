# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **TE (Transposable Element) Fossil Mining Pipeline** for detecting ancient TE insertions in Drosophila 3'UTRs, particularly in germ plasm-localized mRNAs.

**Key documentation:**
- `docs/FILE_MAP.md` - Comprehensive map of all files and their purposes
- `docs/SCRIPTS_REFERENCE.md` - Detailed script documentation
- `docs/SESSION_SUMMARY_diverged_TE_analysis.md` - Current analysis status and roadmap

## File Map Maintenance

**IMPORTANT**: When adding new scripts, data files, or results:
1. Update `docs/FILE_MAP.md` with the new file's purpose
2. If superseding old results, move them to `results/archive/`
3. Document parameters/methods used to generate new outputs

## Current Analysis State

**Best BLAST parameters** (revised 2026-01-22, see `docs/DUST_FILTERING_ANALYSIS.md`):
```
word_size=7, gapopen=2, gapextend=1, penalty=-1, reward=1, dust=no
```

**E-value strategy:**
- `dust=no` is correct — DUST filtering discards real TE signal along with simple repeats.
  The earlier finding that "DUST is critical" was superseded by the full analysis in
  `docs/DUST_FILTERING_ANALYSIS.md` which showed `dust=no` captures 52% more high-quality
  hits with better real-vs-shuffled enrichment ratios.
- **High stringency / quick compute**: `evalue=0.001` — good default for most analyses.
- **Exhaustive search**: `evalue` up to `10` — justified by the hypothesis that TE fossils
  are highly diverged and may recur at multiple loci, so weaker hits in aggregate can be
  informative. Use higher e-values when compute budget allows and when analyzing the full
  "fossil distribution" space.
- Choose e-value based on analysis goals and available compute, not as a fixed parameter.

**Current results location**: `results/` (archived results in `results/archive/`)

## Refactoring Status

**Active refactor**: Migrating from flat `scripts/` to installable `fossil_finder` package.
- Plan: `docs/superpowers/plans/2026-03-10-fossil-finder-phase1-foundation.md`
- New package: `src/fossil_finder/` (config-driven, multi-genome capable)
- Old scripts in `scripts/` remain functional during migration
- Python upgraded from 3.8 → 3.11

## Environment Setup

### New package (recommended)
```bash
uv pip install -e ".[dev]"
```

### With conda environment
```bash
conda env create -f environment.yml
conda activate fossil-finder
uv pip install -e ".[dev]"
```

### Legacy conda env
```bash
# Old env still works for scripts/ during migration
conda activate bioinformatics-program
```

### Using uv (standalone)
```bash
uv pip install -r requirements.txt
```

## Common Commands

### Running analysis scripts
```bash
# Run genome-wide analysis
python scripts/analyze_genome_wide_te.py --blast results/genome_wide_all_3utrs.tsv

# Run gene set comparison
python scripts/rank_gene_sets.py

# Compare 5'UTR vs 3'UTR
python scripts/analyze_utr_variants.py
```

### Running the pipeline
```bash
bash scripts/run_pipeline.sh
```

## Architecture

See `docs/FILE_MAP.md` for the complete file inventory.

### Key Directories

| Directory | Contents | Git Status |
|-----------|----------|------------|
| `scripts/` | Analysis scripts | Tracked |
| `scripts/utils/` | Shared utility modules | Tracked |
| `data/gene_lists/` | Gene sets by tissue/function | Tracked |
| `data/queries/` | Extracted 3'UTR sequences | Tracked |
| `data/references/` | FlyBase reference data | Tracked |
| `results/` | BLAST outputs, analyses | **Gitignored** |
| `results/archive/` | Superseded/intermediate results | **Gitignored** |
| `reports/` | Generated HTML/MD reports | **Gitignored** |
| `figures/` | Generated plots | **Gitignored** |
| `docs/` | Documentation | Tracked |

### Key Scripts

| Task | Script |
|------|--------|
| Run BLAST | `scripts/blast_runner.py` |
| Parameter sweep | `scripts/te_parameter_sweep.py` |
| Build gene lists | `scripts/build_*_genelist.py` |
| Extract 3'UTRs | `scripts/extract_germ_plasm_3utrs.py` |
| Genome-wide analysis | `scripts/analyze_genome_wide_te.py` |
| Rank gene sets | `scripts/rank_gene_sets.py` |
| Compare UTR types | `scripts/analyze_utr_variants.py` |
| Generate reports | `scripts/generate_te_fossil_report.py` |

### Shared Utilities (`scripts/utils/`)

| Module | Purpose |
|--------|---------|
| `paths.py` | Centralized path resolution (use instead of hardcoding) |
| `blast_io.py` | BLAST result parsing (BLAST_COLUMNS, load_blast_results) |
| `data_loaders.py` | Gene list and FASTA loading utilities |

**Usage example:**
```python
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import load_blast_results, BLAST_COLUMNS
from utils.data_loaders import load_gene_list_with_symbols
```

### BLAST Output Format

All TSV results use 17 columns (defined in `utils/blast_io.BLAST_COLUMNS`):
```
qseqid, sseqid, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, evalue, bitscore,
qlen, slen, qseq, sseq, strand
```

### Development Tools
- **notebooks/exploratory.ipynb**: Jupyter notebook for exploratory data analysis
- Dependencies: numpy, pandas, biopython, scikit-learn, matplotlib, seaborn
