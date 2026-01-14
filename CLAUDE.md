# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **TE (Transposable Element) Fossil Mining Pipeline** for detecting ancient TE insertions in Drosophila 3'UTRs, particularly in germ plasm-localized mRNAs.

**Key documentation:**
- `docs/FILE_MAP.md` - Comprehensive map of all files and their purposes
- `docs/SESSION_SUMMARY_diverged_TE_analysis.md` - Current analysis status and roadmap

## File Map Maintenance

**IMPORTANT**: When adding new scripts, data files, or results:
1. Update `docs/FILE_MAP.md` with the new file's purpose
2. If superseding old results, note this in the "Superseded Results" section
3. Document parameters/methods used to generate new outputs

## Current Analysis State

**Best BLAST parameters** (as of 2026-01-13):
```
word_size=7, gapopen=2, gapextend=1, penalty=-1, reward=1, dust=yes
```

**Current results location**: `results/diverged_controls/`

**Key finding**: DUST filtering is critical - without it, 90% of hits are simple repeats.

## Environment Setup

### Using Conda (recommended)
```bash
conda env create -f environment.yml
conda activate bioinformatics-program
```

### Using pip
```bash
pip install -r requirements.txt
```

## Common Commands

### Running the program
```bash
python src/bioinfo/main.py
```

### Running tests
```bash
# Run all tests
pytest src/tests/

# Run specific test file
pytest src/tests/test_core.py

# Run with verbose output
pytest -v src/tests/
```

### Running the pipeline
```bash
bash scripts/run_pipeline.sh
```

### Running Snakemake workflows
```bash
snakemake -s workflows/Snakefile --cores all
```

## Architecture

See `docs/FILE_MAP.md` for the complete file inventory.

### Key Directories

| Directory | Contents | Git Status |
|-----------|----------|------------|
| `scripts/` | Analysis scripts | Tracked |
| `data/gene_lists/` | Gene sets by tissue/function | Tracked |
| `data/queries/` | Extracted 3'UTR sequences | Tracked |
| `data/references/` | FlyBase reference data | Tracked |
| `results/` | BLAST outputs, analyses | **Gitignored** |
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
| Analyze results | `scripts/analyze_te_families.py`, `scripts/extract_top_hits.py` |
| Generate reports | `scripts/generate_te_fossil_report.py` |

### BLAST Output Format

All TSV results use 17 columns:
```
qseqid, sseqid, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, evalue, bitscore,
qlen, slen, qseq, sseq, strand
```

### Development Tools
- **notebooks/exploratory.ipynb**: Jupyter notebook for exploratory data analysis
- Dependencies: numpy, pandas, biopython, scikit-learn, matplotlib, seaborn
