"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density
from .conservation import (
    compute_pident_conservation_correlation,
    hits_to_genomic_bed,
    score_with_bigwig,
    summarize_conservation_by_group,
)
from .enrichment import (
    correct_multiple_testing,
    fisher_exact_enrichment,
    mannwhitney_enrichment,
    test_gene_set_enrichment,
)
from .families import compute_class_distribution, compute_family_stats, compute_fold_enrichment
from .multiplicity import compute_hit_multiplicity, compute_query_hit_counts, compute_te_breadth
from .positional import compute_end_bias, compute_positional_profile, compute_te_position, compute_utr_position
from .quality_tiers import assign_quality_tiers, compute_edit_stats, compute_tier_edit_summary, summarize_tiers
from .repeatmasker import classify_hits, find_overlaps, parse_repeatmasker_out
from .strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)

__all__ = [
    "aggregate_by_gene",
    "assign_quality_tiers",
    "classify_bias",
    "classify_hits",
    "compute_class_distribution",
    "compute_density",
    "compute_edit_stats",
    "compute_end_bias",
    "compute_family_stats",
    "compute_fold_enrichment",
    "compute_gene_strand_bias",
    "compute_genome_strand_bias",
    "compute_hit_multiplicity",
    "compute_pident_conservation_correlation",
    "compute_positional_profile",
    "compute_query_hit_counts",
    "compute_te_breadth",
    "compute_te_position",
    "compute_te_strand_bias",
    "compute_tier_edit_summary",
    "compute_utr_position",
    "correct_multiple_testing",
    "find_overlaps",
    "fisher_exact_enrichment",
    "hits_to_genomic_bed",
    "mannwhitney_enrichment",
    "parse_repeatmasker_out",
    "score_with_bigwig",
    "summarize_conservation_by_group",
    "summarize_tiers",
    "test_gene_set_enrichment",
]
