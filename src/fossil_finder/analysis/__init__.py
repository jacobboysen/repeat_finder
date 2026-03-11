"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density
from .families import compute_class_distribution, compute_family_stats, compute_fold_enrichment
from .strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)

__all__ = [
    "aggregate_by_gene",
    "classify_bias",
    "compute_class_distribution",
    "compute_density",
    "compute_family_stats",
    "compute_fold_enrichment",
    "compute_gene_strand_bias",
    "compute_genome_strand_bias",
    "compute_te_strand_bias",
]
