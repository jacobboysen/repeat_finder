"""Per-gene BLAST hit aggregation.

Maps BLAST hits (keyed by query/transcript) to genes and computes
per-gene summary metrics: hit_count, hit_bp, density, strand counts,
and TE family diversity.
"""

import pandas as pd


def aggregate_by_gene(
    df: pd.DataFrame,
    query_to_gene: dict[str, str],
) -> pd.DataFrame:
    """Aggregate BLAST hits by gene.

    Args:
        df: BLAST results DataFrame with at least columns:
            qseqid, length, strand, qlen, sseqid.
        query_to_gene: Mapping from query ID (qseqid) to gene ID.

    Returns:
        DataFrame indexed by gene_id with columns:
        hit_count, hit_bp, sense_hits, antisense_hits,
        query_len, n_te_families.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["hit_count", "hit_bp", "sense_hits", "antisense_hits",
                      "query_len", "n_te_families"],
        )

    # Map queries to genes, drop unmapped
    work = df.copy()
    work["gene_id"] = work["qseqid"].map(query_to_gene)
    work = work.dropna(subset=["gene_id"])

    if work.empty:
        return pd.DataFrame(
            columns=["hit_count", "hit_bp", "sense_hits", "antisense_hits",
                      "query_len", "n_te_families"],
        )

    grouped = work.groupby("gene_id")

    result = pd.DataFrame({
        "hit_count": grouped.size(),
        "hit_bp": grouped["length"].sum(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
        "query_len": grouped["qlen"].max(),
        "n_te_families": grouped["sseqid"].nunique(),
    })

    return result


def compute_density(gene_stats: pd.DataFrame) -> pd.DataFrame:
    """Add TE density column (hit_bp per kb of query sequence).

    Args:
        gene_stats: DataFrame with hit_bp and query_len columns.

    Returns:
        Same DataFrame with added 'density' column.
    """
    result = gene_stats.copy()
    mask = result["query_len"] > 0
    result["density"] = 0.0
    result.loc[mask, "density"] = (
        result.loc[mask, "hit_bp"] / result.loc[mask, "query_len"] * 1000
    )
    return result
