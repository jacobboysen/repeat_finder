"""Statistical enrichment testing for TE fossil analysis.

Provides Fisher's exact test (TE presence enrichment), Mann-Whitney U
(density distribution comparison), and multiple testing correction.
"""

import math

import numpy as np
from scipy import stats


def fisher_exact_enrichment(
    set_positive: int,
    set_negative: int,
    bg_positive: int,
    bg_negative: int,
) -> dict:
    """Fisher's exact test for TE presence enrichment.

    Tests whether a gene set has more TE-positive genes than expected
    by the background rate.

    Args:
        set_positive: Genes in set WITH TE hits.
        set_negative: Genes in set WITHOUT TE hits.
        bg_positive: Background genes WITH TE hits (excluding set).
        bg_negative: Background genes WITHOUT TE hits (excluding set).

    Returns:
        Dict with odds_ratio, p_value.
    """
    table = [[set_positive, set_negative], [bg_positive, bg_negative]]
    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
    return {"odds_ratio": odds_ratio, "p_value": p_value}


def mannwhitney_enrichment(
    group_a: list[float],
    group_b: list[float],
) -> dict:
    """Mann-Whitney U test for density distribution comparison.

    Args:
        group_a: Density values for gene set.
        group_b: Density values for background.

    Returns:
        Dict with u_statistic, p_value.
    """
    if len(group_a) < 2 or len(group_b) < 2:
        return {"u_statistic": float("nan"), "p_value": float("nan")}

    u_stat, p_value = stats.mannwhitneyu(
        group_a, group_b, alternative="two-sided",
    )
    return {"u_statistic": u_stat, "p_value": p_value}


def correct_multiple_testing(
    p_values: list[float],
    method: str = "fdr",
) -> list[float]:
    """Apply multiple testing correction.

    Args:
        p_values: List of raw p-values.
        method: "fdr" (Benjamini-Hochberg) or "bonferroni".

    Returns:
        List of corrected p-values (q-values).
    """
    if not p_values:
        return []

    n = len(p_values)

    if method == "bonferroni":
        return [min(p * n, 1.0) for p in p_values]

    # Benjamini-Hochberg FDR
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [0.0] * n

    for rank, (orig_idx, p) in enumerate(indexed, 1):
        q = min(p * n / rank, 1.0)
        q_values[orig_idx] = q

    # Enforce monotonicity (q-values must be non-decreasing in sorted order)
    sorted_indices = [i for i, _ in indexed]
    for k in range(n - 2, -1, -1):
        idx = sorted_indices[k]
        next_idx = sorted_indices[k + 1]
        if q_values[idx] > q_values[next_idx]:
            q_values[idx] = q_values[next_idx]

    return q_values


def test_gene_set_enrichment(
    gene_set: set[str],
    te_positive_genes: set[str],
    gene_densities: "pd.Series",
) -> dict:
    """Run full enrichment analysis for a gene set.

    Combines Fisher's exact (TE presence) and Mann-Whitney U (density)
    into a single result.

    Args:
        gene_set: Gene IDs in the set of interest.
        te_positive_genes: All gene IDs with at least one TE hit.
        gene_densities: Series mapping gene_id -> density value.

    Returns:
        Dict with 'fisher' and 'mannwhitney' sub-dicts.
    """
    all_genes = set(gene_densities.index)

    # Fisher's exact: TE presence
    set_pos = len(gene_set & te_positive_genes)
    set_neg = len(gene_set - te_positive_genes)
    bg_genes = all_genes - gene_set
    bg_pos = len(bg_genes & te_positive_genes)
    bg_neg = len(bg_genes - te_positive_genes)

    fisher_result = fisher_exact_enrichment(set_pos, set_neg, bg_pos, bg_neg)

    # Mann-Whitney U: density distribution
    set_densities = [
        gene_densities[g] for g in gene_set if g in gene_densities.index
    ]
    bg_densities = [
        gene_densities[g] for g in bg_genes if g in gene_densities.index
    ]

    mw_result = mannwhitney_enrichment(set_densities, bg_densities)

    return {"fisher": fisher_result, "mannwhitney": mw_result}


# Prevent pytest from collecting this function as a test
test_gene_set_enrichment.__test__ = False
