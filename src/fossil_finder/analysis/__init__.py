"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density
from .enrichment import (
    correct_multiple_testing,
    fisher_exact_enrichment,
    mannwhitney_enrichment,
    test_gene_set_enrichment,
)
from .families import compute_class_distribution, compute_family_stats, compute_fold_enrichment
from .repeatmasker import classify_hits, find_overlaps, parse_repeatmasker_out
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
    "classify_hits",
    "correct_multiple_testing",
    "find_overlaps",
    "fisher_exact_enrichment",
    "mannwhitney_enrichment",
    "parse_repeatmasker_out",
    "test_gene_set_enrichment",
]
