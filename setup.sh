#!/bin/bash
#
# Setup script for Drosophila BLAST environment
#
# This script:
# 1. Creates necessary directories
# 2. Downloads reference files from FlyBase and Dfam
# 3. Builds BLAST databases
#
# Usage:
#   bash setup.sh [options]
#
# Options:
#   --skip-download    Skip downloading reference files
#   --skip-databases   Skip building BLAST databases
#   --help            Show this help message

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
SKIP_DOWNLOAD=false
SKIP_DATABASES=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-databases)
            SKIP_DATABASES=true
            shift
            ;;
        --help)
            head -n 15 "$0" | tail -n +2 | sed 's/^# //'
            exit 0
            ;;
        *)
            echo -e "${RED}Error: Unknown option: $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Print banner
echo "=========================================================================="
echo "  Drosophila BLAST Environment Setup"
echo "=========================================================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo -e "${YELLOW}Warning: conda not found${NC}"
    echo "You may need to install conda or activate your conda environment"
    echo ""
fi

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo -e "${RED}Error: Python not found${NC}"
    echo "Please install Python or activate your conda environment"
    exit 1
fi

# Step 1: Create directories
echo -e "${BLUE}[1/3] Creating directory structure...${NC}"
mkdir -p data/references
mkdir -p data/blastdb
mkdir -p data/queries
mkdir -p results
echo -e "${GREEN}✓ Directories created${NC}"
echo ""

# Step 2: Download reference files
if [ "$SKIP_DOWNLOAD" = false ]; then
    echo -e "${BLUE}[2/3] Downloading reference files from FlyBase and Dfam...${NC}"
    echo "This may take a while depending on your internet connection."
    echo ""

    if python scripts/download_references.py; then
        echo -e "${GREEN}✓ Downloads completed successfully${NC}"
    else
        echo -e "${RED}✗ Download failed${NC}"
        echo "You can retry later by running:"
        echo "  python scripts/download_references.py"
        exit 1
    fi
    echo ""
else
    echo -e "${YELLOW}[2/3] Skipping downloads (--skip-download)${NC}"
    echo ""
fi

# Step 3: Build BLAST databases
if [ "$SKIP_DATABASES" = false ]; then
    echo -e "${BLUE}[3/3] Building BLAST databases...${NC}"
    echo ""

    # Check if makeblastdb is available
    if ! command -v makeblastdb &> /dev/null; then
        echo -e "${RED}Error: makeblastdb not found${NC}"
        echo "Please install BLAST+ with:"
        echo "  conda install -c bioconda blast"
        exit 1
    fi

    if python scripts/build_databases.py; then
        echo -e "${GREEN}✓ Databases built successfully${NC}"
    else
        echo -e "${RED}✗ Database build failed${NC}"
        echo "You can retry later by running:"
        echo "  python scripts/build_databases.py"
        exit 1
    fi
    echo ""
else
    echo -e "${YELLOW}[3/3] Skipping database build (--skip-databases)${NC}"
    echo ""
fi

# Summary
echo "=========================================================================="
echo -e "${GREEN}Setup Complete!${NC}"
echo "=========================================================================="
echo ""
echo "Directory structure:"
echo "  data/references/  - Reference FASTA and annotation files"
echo "  data/blastdb/     - BLAST databases"
echo "  data/queries/     - Filtered query sequences (gene groups)"
echo "  results/          - BLAST results"
echo ""
echo "Next steps:"
echo ""
echo "1. List available databases:"
echo "   python scripts/list_databases.py"
echo ""
echo "2. List available gene groups:"
echo "   python scripts/list_gene_groups.py"
echo ""
echo "3. Filter queries by gene group:"
echo "   python scripts/filter_queries.py --gene-group \"GERM PLASM\""
echo ""
echo "4. Run BLAST:"
echo "   python scripts/blast_runner.py --segment-type CDS --gene-group \"GERM PLASM\""
echo ""
echo "5. Summarize results:"
echo "   python scripts/summarize_results.py results/.../blast_results.tsv --plot"
echo ""
echo "For help with any script, use --help:"
echo "   python scripts/blast_runner.py --help"
echo ""
echo "=========================================================================="
