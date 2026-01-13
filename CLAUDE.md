# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

This is a modular bioinformatics program with the following key components:

### Core Modules (src/bioinfo/)
- **main.py**: Entry point that orchestrates program execution
- **core.py**: Contains core algorithms and data processing functions
- **io.py**: Handles all file I/O operations with support for CSV and JSON formats

### Data Organization
- **data/raw/**: Unprocessed biological data files
- **data/processed/**: Cleaned and transformed data ready for analysis

### Workflow System
The program uses Snakemake for defining data processing pipelines. Workflows are defined in `workflows/Snakefile` and coordinated via the `scripts/run_pipeline.sh` script.

### Testing
Tests are located in `src/tests/` and use unittest. The modular architecture allows each component to be tested independently.

### Development Tools
- **notebooks/exploratory.ipynb**: Jupyter notebook for exploratory data analysis and visualization
- Dependencies include numpy, pandas, biopython, scikit-learn, matplotlib, and seaborn
