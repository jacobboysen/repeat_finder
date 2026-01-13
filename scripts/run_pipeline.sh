#!/bin/bash

# This script is used to run the entire bioinformatics pipeline.

# Activate the conda environment
# conda activate your_environment_name

# Run the main bioinformatics program
python src/bioinfo/main.py

# Execute the Snakemake workflow
# snakemake -s workflows/Snakefile --cores all

# Add any additional commands needed to run the pipeline here.