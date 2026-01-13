# Drosophila BLAST Environment

## Overview

This project provides a comprehensive local BLAST environment for analyzing *Drosophila melanogaster* genic segments (CDS, 5'UTR, 3'UTR, introns) against genome and transposable element (TE) databases. It supports filtering queries by FlyBase Gene Groups and provides full parameter control via config file with CLI overrides.

## Features

- **Automated Setup**: One-command setup downloads references and builds BLAST databases
- **Gene Group Filtering**: Filter sequences by FlyBase Gene Groups for targeted analysis
- **Flexible Configuration**: YAML config file with complete CLI override support
- **Multiple Databases**: Query against genome, FlyBase TEs, Dfam TEs, or combined TE database
- **Comprehensive Analysis**: Built-in tools for result summarization and visualization
- **Organized Output**: Results automatically organized by segment type, database, and timestamp

## Quick Start

### 1. Environment Setup

Create the conda environment with all dependencies:

```bash
conda env create -f environment.yml
conda activate bioinformatics-program
```

### 2. Download Data and Build Databases

Run the automated setup script:

```bash
bash setup.sh
```

This will:
- Download genome, annotations, and genic segments from FlyBase
- Download TE sequences from FlyBase and Dfam
- Build BLAST databases for genome and TEs
- Create combined, deduplicated TE database

### 3. Run Your First BLAST Search

```bash
# Query CDS sequences against genome
python scripts/blast_runner.py

# Query specific gene group against TE database
python scripts/blast_runner.py --segment-type CDS --gene-group "GERM PLASM" --database dmel_te_dfam

# Use custom parameters
python scripts/blast_runner.py --evalue 1e-10 --word-size 7 --num-threads 8
```

## Detailed Usage

### Listing Available Resources

**List BLAST databases:**
```bash
python scripts/list_databases.py
```

**List gene groups:**
```bash
python scripts/list_gene_groups.py

# Search for specific groups
python scripts/list_gene_groups.py --search "germ"

# Show only large groups
python scripts/list_gene_groups.py --min-genes 50
```

### Filtering Queries by Gene Group

Before running BLAST on a specific gene group, filter the sequences:

```bash
# Filter all segment types for a gene group
python scripts/filter_queries.py --gene-group "GERM PLASM"

# Filter only CDS
python scripts/filter_queries.py --gene-group "GERM PLASM" --segment-types CDS

# List available groups
python scripts/filter_queries.py --list-groups
```

Filtered sequences are saved to `data/queries/{gene_group}/{segment_type}.fasta`

### Running BLAST

The BLAST runner supports both config file and CLI configuration:

**Using config file (config.yaml):**
```bash
python scripts/blast_runner.py
```

**Override with CLI arguments:**
```bash
python scripts/blast_runner.py \
  --segment-type CDS \
  --gene-group "GERM PLASM" \
  --database dmel_te_combined \
  --evalue 1e-10 \
  --word-size 7 \
  --num-threads 8
```

**Use custom config file:**
```bash
python scripts/blast_runner.py --config sensitive.yaml
```

### Analyzing Results

**Summarize BLAST results:**
```bash
python scripts/summarize_results.py results/CDS/dmel_genome/20240115_120000_GERM_PLASM/blast_results.tsv

# With distribution plots
python scripts/summarize_results.py results/.../blast_results.tsv --plot

# Show top 20 hits
python scripts/summarize_results.py results/.../blast_results.tsv --top 20
```

## Configuration

The `config.yaml` file contains default BLAST parameters:

```yaml
blast:
  program: blastn
  evalue: 1e-5
  word_size: 11
  dust: "yes"
  perc_identity: 0
  max_target_seqs: 500
  max_hsps: 10
  outfmt: "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
  num_threads: 4

query:
  segment_type: CDS
  gene_group: null

subject:
  database: dmel_genome
```

All parameters can be overridden via command line.

## Directory Structure

```
├── config.yaml              # BLAST configuration
├── environment.yml          # Conda environment
├── setup.sh                # Automated setup script
├── scripts/
│   ├── download_references.py   # Download FlyBase/Dfam data
│   ├── build_databases.py       # Build BLAST databases
│   ├── filter_queries.py        # Filter by gene groups
│   ├── blast_runner.py          # Main BLAST runner
│   ├── list_databases.py        # List available databases
│   ├── list_gene_groups.py      # List gene groups
│   └── summarize_results.py     # Analyze BLAST results
├── data/
│   ├── references/          # Downloaded FASTA and annotations
│   ├── blastdb/            # BLAST databases
│   └── queries/            # Filtered query sequences
└── results/                # BLAST results organized by date
```

## Output Structure

Results are automatically organized:

```
results/{segment_type}/{database}/{timestamp}_{gene_group}/
├── blast_results.tsv    # BLAST hits (tabular format)
├── config_used.yaml     # Configuration for this run
├── query_stats.txt      # Query sequence statistics
└── run.log             # Execution log
```

## Scripts Reference

### download_references.py
Downloads reference files from FlyBase and Dfam.

```bash
# Download all files
python scripts/download_references.py

# Download specific files
python scripts/download_references.py --files genome cds transposons

# List available files
python scripts/download_references.py --list
```

### build_databases.py
Builds BLAST databases from reference FASTA files.

```bash
# Build all databases
python scripts/build_databases.py

# Build specific databases
python scripts/build_databases.py --databases genome te_combined

# Rebuild existing databases
python scripts/build_databases.py --force
```

### filter_queries.py
Filters genic segment sequences by FlyBase Gene Group membership.

```bash
# Filter by gene group
python scripts/filter_queries.py --gene-group "GERM PLASM"

# Filter specific segment types
python scripts/filter_queries.py --gene-group "GERM PLASM" --segment-types CDS 5UTR

# List all available groups
python scripts/filter_queries.py --list-groups
```

### blast_runner.py
Runs BLAST with flexible configuration management.

```bash
# Use default config
python scripts/blast_runner.py

# Override parameters
python scripts/blast_runner.py --evalue 1e-10 --num-threads 8

# Use different config file
python scripts/blast_runner.py --config sensitive.yaml

# Run on gene group
python scripts/blast_runner.py --segment-type CDS --gene-group "GERM PLASM"
```

### list_databases.py
Lists available BLAST databases with statistics.

```bash
# List all databases
python scripts/list_databases.py

# Show detailed information
python scripts/list_databases.py --verbose
```

### list_gene_groups.py
Lists FlyBase gene groups with gene counts.

```bash
# List all groups
python scripts/list_gene_groups.py

# Search for specific groups
python scripts/list_gene_groups.py --search "germ"

# Sort by name
python scripts/list_gene_groups.py --sort-by name

# Filter by minimum size
python scripts/list_gene_groups.py --min-genes 100
```

### summarize_results.py
Analyzes BLAST results with statistics and visualizations.

```bash
# Basic summary
python scripts/summarize_results.py results/.../blast_results.tsv

# Generate distribution plots
python scripts/summarize_results.py results/.../blast_results.tsv --plot

# Show top 20 queries/subjects
python scripts/summarize_results.py results/.../blast_results.tsv --top 20
```

## Requirements

- Python 3.8+
- BLAST+ (via conda)
- Biopython
- PyYAML
- pandas
- numpy
- matplotlib (for plotting)

All dependencies are specified in `environment.yml`.

## Data Sources

- **FlyBase** (r6.55): Genome, annotations, genic segments, gene groups
- **Dfam** (3.7): Transposable element consensus sequences

## Troubleshooting

**BLAST not found:**
```bash
conda install -c bioconda blast
```

**Missing reference files:**
```bash
python scripts/download_references.py
```

**Database not found:**
```bash
python scripts/build_databases.py
```

**Gene group not found:**
```bash
# List available groups
python scripts/list_gene_groups.py
```

## Contributing

Contributions are welcome! Please submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License.