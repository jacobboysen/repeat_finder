"""
Centralized path resolution for the TE fossil mining pipeline.

All scripts should use these functions instead of hardcoding paths.
This ensures portability across different machines and setups.
"""

from pathlib import Path


def get_project_root() -> Path:
    """
    Get the project root directory.

    Works by finding the directory containing this utils package,
    then going up to the project root.

    Returns:
        Path to project root (repeat_finder/)
    """
    # This file is at: scripts/utils/paths.py
    # Project root is two levels up
    return Path(__file__).parent.parent.parent.resolve()


def get_data_dir() -> Path:
    """Get the data directory."""
    return get_project_root() / 'data'


def get_results_dir() -> Path:
    """Get the results directory."""
    return get_project_root() / 'results'


def get_gene_lists_dir() -> Path:
    """Get the gene lists directory."""
    return get_data_dir() / 'gene_lists'


def get_queries_dir() -> Path:
    """Get the queries directory."""
    return get_data_dir() / 'queries'


def get_references_dir() -> Path:
    """Get the references directory."""
    return get_data_dir() / 'references'


def get_blastdb_dir() -> Path:
    """Get the BLAST database directory."""
    return get_data_dir() / 'blastdb'


def get_figures_dir() -> Path:
    """Get the figures directory."""
    return get_project_root() / 'figures'


def get_reports_dir() -> Path:
    """Get the reports directory."""
    return get_project_root() / 'reports'


# Common reference file paths
def get_te_database_path() -> Path:
    """Get path to FlyBase TE database FASTA."""
    return get_references_dir() / 'dmel_te_flybase.fasta'


def get_te_consensus_path() -> Path:
    """Get path to TE consensus sequences FASTA."""
    return get_references_dir() / 'dmel_te_consensus.fasta'


def get_3utr_fasta_path() -> Path:
    """Get path to 3'UTR reference FASTA."""
    return get_references_dir() / 'dmel_3utr.fasta'


def get_5utr_fasta_path() -> Path:
    """Get path to 5'UTR reference FASTA."""
    return get_references_dir() / 'dmel_5utr.fasta'


def get_gene_annotation_path() -> Path:
    """Get path to gene annotation file."""
    return get_references_dir() / 'fbgn_annotation_ID.tsv'


# Gene list file paths
def get_gene_list_path(group: str, file_type: str = 'consolidated') -> Path:
    """
    Get path to a gene list file.

    Args:
        group: Gene set name (germ_plasm, housekeeping, somatic, cleared, adult)
        file_type: File type (consolidated, fbgn_ids, status)

    Returns:
        Path to the gene list file
    """
    if file_type == 'consolidated':
        return get_gene_lists_dir() / f'{group}_genes_consolidated.tsv'
    elif file_type == 'fbgn_ids':
        return get_gene_lists_dir() / f'{group}_fbgn_ids.txt'
    elif file_type == 'status':
        return get_gene_lists_dir() / f'{group}_status.json'
    else:
        raise ValueError(f"Unknown file_type: {file_type}")


# Query sequence paths
def get_query_fasta_path(group: str, strand: str = 'sense', tier: str = None) -> Path:
    """
    Get path to a query FASTA file.

    Args:
        group: Gene set name
        strand: 'sense', 'antisense', or 'shuffled'
        tier: Optional tier (e.g., 'tier1')

    Returns:
        Path to the query FASTA file
    """
    queries_dir = get_queries_dir() / group

    if strand == 'shuffled':
        return queries_dir / '3UTR_shuffled.fasta'

    filename = f'3UTR_{strand}'
    if tier:
        filename += f'_{tier}'
    filename += '.fasta'

    return queries_dir / filename
