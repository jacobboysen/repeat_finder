"""Strand bias analysis for TE fossil hits.

Computes sense/antisense bias at gene, TE-family, and genome-wide levels.
Strand convention: "plus" = sense (same orientation as query), "minus" = antisense.
"""

import pandas as pd


def classify_bias(sense_fraction: float) -> str:
    """Classify strand bias from sense fraction.

    Args:
        sense_fraction: Fraction of sense hits [0.0, 1.0].

    Returns:
        One of: "strong_sense", "sense", "balanced",
        "antisense", "strong_antisense".
    """
    if sense_fraction >= 0.70:
        return "strong_sense"
    elif sense_fraction >= 0.55:
        return "sense"
    elif sense_fraction >= 0.45:
        return "balanced"
    elif sense_fraction >= 0.30:
        return "antisense"
    else:
        return "strong_antisense"


def compute_gene_strand_bias(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-gene strand bias.

    Args:
        df: BLAST results with gene_id and strand columns.

    Returns:
        DataFrame indexed by gene_id with columns:
        total_hits, sense_hits, antisense_hits, sense_pct, bias.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["total_hits", "sense_hits", "antisense_hits",
                      "sense_pct", "bias"],
        )

    grouped = df.groupby("gene_id")

    result = pd.DataFrame({
        "total_hits": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
    })

    result["sense_pct"] = result.apply(
        lambda r: (r["sense_hits"] / r["total_hits"] * 100)
        if r["total_hits"] > 0 else 0.0,
        axis=1,
    )
    result["bias"] = result["sense_pct"].apply(
        lambda pct: classify_bias(pct / 100)
    )

    return result


def compute_te_strand_bias(
    df: pd.DataFrame, min_hits: int = 1,
) -> pd.DataFrame:
    """Compute per-TE-family strand bias.

    Args:
        df: BLAST results with sseqid and strand columns.
        min_hits: Minimum total hits to include a TE family.

    Returns:
        DataFrame indexed by TE ID with strand bias columns.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["total_hits", "sense_hits", "antisense_hits",
                      "sense_pct", "bias"],
        )

    grouped = df.groupby("sseqid")

    result = pd.DataFrame({
        "total_hits": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
    })

    result = result[result["total_hits"] >= min_hits]

    result["sense_pct"] = result.apply(
        lambda r: (r["sense_hits"] / r["total_hits"] * 100)
        if r["total_hits"] > 0 else 0.0,
        axis=1,
    )
    result["bias"] = result["sense_pct"].apply(
        lambda pct: classify_bias(pct / 100)
    )

    return result


def compute_genome_strand_bias(df: pd.DataFrame) -> dict:
    """Compute genome-wide strand bias summary.

    Args:
        df: BLAST results with strand column.

    Returns:
        Dict with total_hits, sense_hits, antisense_hits, sense_pct, bias.
    """
    if df.empty:
        return {
            "total_hits": 0, "sense_hits": 0, "antisense_hits": 0,
            "sense_pct": 0.0, "bias": "balanced",
        }

    total = len(df)
    sense = (df["strand"] == "plus").sum()
    antisense = (df["strand"] == "minus").sum()
    pct = sense / total * 100 if total > 0 else 0.0

    return {
        "total_hits": total,
        "sense_hits": int(sense),
        "antisense_hits": int(antisense),
        "sense_pct": pct,
        "bias": classify_bias(pct / 100),
    }
