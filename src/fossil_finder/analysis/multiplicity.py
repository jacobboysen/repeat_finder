"""Hit multiplicity analysis for TE fossil mining.

Computes distribution statistics for how many BLAST hits each entity
(query, gene, TE family) accumulates, and per-TE breadth summaries
showing how broadly each TE hits across queries and genes.
"""

import pandas as pd


def _summarize_series(series: pd.Series) -> dict:
    """Compute standard summary statistics for a count series.

    Returns:
        Dict with median (float), mean (float), max (int), n (int).
    """
    return {
        "median": float(series.median()),
        "mean": float(series.mean()),
        "max": int(series.max()),
        "n": int(len(series)),
    }


def compute_hit_multiplicity(
    df: pd.DataFrame,
    query_to_gene: dict[str, str] | None = None,
) -> dict:
    """Compute distribution statistics for hit multiplicity.

    Measures how many hits each query/gene accumulates and how many
    unique queries/genes each TE family hits.

    Args:
        df: BLAST results with at least columns: qseqid, sseqid.
        query_to_gene: Optional mapping from qseqid to gene ID.
            When provided, gene-level metrics are included.

    Returns:
        Dict with keys:
        - "hits_per_query": summary stats (median, mean, max, n)
        - "queries_per_te": summary stats (median, mean, max, n)
        - "genes_per_te": summary stats (only if query_to_gene provided)
        - "hits_per_gene": summary stats (only if query_to_gene provided)
    """
    if df.empty:
        empty = {"median": 0.0, "mean": 0.0, "max": 0, "n": 0}
        result = {
            "hits_per_query": empty.copy(),
            "queries_per_te": empty.copy(),
        }
        if query_to_gene is not None:
            result["genes_per_te"] = empty.copy()
            result["hits_per_gene"] = empty.copy()
        return result

    hits_per_query = df.groupby("qseqid")["sseqid"].count()
    queries_per_te = df.groupby("sseqid")["qseqid"].nunique()

    result: dict = {
        "hits_per_query": _summarize_series(hits_per_query),
        "queries_per_te": _summarize_series(queries_per_te),
    }

    if query_to_gene is not None:
        work = df.copy()
        work["gene_id"] = work["qseqid"].map(query_to_gene)
        mapped = work.dropna(subset=["gene_id"])

        if mapped.empty:
            empty = {"median": 0.0, "mean": 0.0, "max": 0, "n": 0}
            result["genes_per_te"] = empty.copy()
            result["hits_per_gene"] = empty.copy()
        else:
            genes_per_te = mapped.groupby("sseqid")["gene_id"].nunique()
            hits_per_gene = mapped.groupby("gene_id")["sseqid"].count()
            result["genes_per_te"] = _summarize_series(genes_per_te)
            result["hits_per_gene"] = _summarize_series(hits_per_gene)

    return result


def compute_te_breadth(
    df: pd.DataFrame,
    query_to_gene: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Per-TE summary of how broadly each TE hits across queries/genes.

    Args:
        df: BLAST results with at least columns: qseqid, sseqid, pident.
        query_to_gene: Optional mapping from qseqid to gene ID.
            When provided, an n_genes column is included.

    Returns:
        DataFrame indexed by sseqid with columns:
        n_hits, n_queries, n_genes (if query_to_gene provided),
        median_pident, mean_pident.
    """
    if df.empty:
        cols = ["n_hits", "n_queries", "median_pident", "mean_pident"]
        if query_to_gene is not None:
            cols.insert(2, "n_genes")
        return pd.DataFrame(columns=cols)

    grouped = df.groupby("sseqid")

    result = pd.DataFrame({
        "n_hits": grouped["qseqid"].count(),
        "n_queries": grouped["qseqid"].nunique(),
        "median_pident": grouped["pident"].median(),
        "mean_pident": grouped["pident"].mean(),
    })

    if query_to_gene is not None:
        work = df.copy()
        work["gene_id"] = work["qseqid"].map(query_to_gene)
        mapped = work.dropna(subset=["gene_id"])
        if mapped.empty:
            result.insert(2, "n_genes", 0)
        else:
            n_genes = mapped.groupby("sseqid")["gene_id"].nunique()
            result.insert(2, "n_genes", n_genes)
            # TEs with no mapped genes get 0
            result["n_genes"] = result["n_genes"].fillna(0).astype(int)

    return result


def compute_query_hit_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Count BLAST hits per query sequence.

    Args:
        df: BLAST results with at least a qseqid column.

    Returns:
        DataFrame indexed by qseqid with a single hit_count column.
    """
    if df.empty:
        return pd.DataFrame(columns=["hit_count"])

    counts = df.groupby("qseqid").size().rename("hit_count")
    return counts.to_frame()
