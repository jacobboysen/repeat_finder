"""
Shared utilities for TE fossil mining pipeline.

This module provides common functionality used across analysis scripts:
- paths: Centralized path resolution
- blast_io: BLAST result parsing
- data_loaders: Gene list and FASTA loading
"""

from .paths import get_project_root, get_data_dir, get_results_dir
from .blast_io import BLAST_COLUMNS, load_blast_results, classify_strand
from .data_loaders import (
    load_gene_list,
    load_gene_list_with_symbols,
    parse_fasta,
    parse_fasta_by_parent,
    load_te_database,
)
