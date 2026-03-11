"""TE family distribution and fold-enrichment analysis.

Computes per-TE-family hit statistics and comparative enrichment
between gene sets (e.g., germ plasm vs housekeeping).
"""

import math

import pandas as pd


def compute_family_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-TE-family summary statistics.

    Args:
        df: BLAST results with sseqid, strand, pident, evalue, length columns.

    Returns:
        DataFrame indexed by TE ID with columns:
        hit_count, sense_hits, antisense_hits, mean_pident,
        mean_evalue, total_bp, frequency.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["hit_count", "sense_hits", "antisense_hits",
                      "mean_pident", "mean_evalue", "total_bp", "frequency"],
        )

    grouped = df.groupby("sseqid")
    total = len(df)

    result = pd.DataFrame({
        "hit_count": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
        "mean_pident": grouped["pident"].mean(),
        "mean_evalue": grouped["evalue"].mean(),
        "total_bp": grouped["length"].sum(),
        "frequency": grouped.size() / total,
    })

    return result


def compute_class_distribution(
    df: pd.DataFrame,
    te_metadata: dict[str, dict],
) -> dict[str, int]:
    """Count BLAST hits by TE class.

    Args:
        df: BLAST results with sseqid column.
        te_metadata: Dict mapping TE ID to metadata dict
            containing at least a 'te_class' key.

    Returns:
        Dict mapping TE class name to hit count.
    """
    counts: dict[str, int] = {}
    for te_id in df["sseqid"]:
        meta = te_metadata.get(te_id, {})
        te_class = meta.get("te_class", "Unknown")
        counts[te_class] = counts.get(te_class, 0) + 1
    return counts


def compute_fold_enrichment(
    set_a: pd.DataFrame,
    set_b: pd.DataFrame,
) -> pd.DataFrame:
    """Compute fold-enrichment of TE families between two hit sets.

    For each TE family, fold_enrichment = freq_in_a / freq_in_b.
    A value > 1 means enriched in set A; < 1 means enriched in set B.

    Args:
        set_a: BLAST results for gene set A.
        set_b: BLAST results for gene set B.

    Returns:
        DataFrame indexed by TE ID with columns:
        count_a, count_b, freq_a, freq_b, fold_enrichment, log2_enrichment.
    """
    stats_a = compute_family_stats(set_a) if not set_a.empty else pd.DataFrame()
    stats_b = compute_family_stats(set_b) if not set_b.empty else pd.DataFrame()

    all_families = set()
    if not stats_a.empty:
        all_families.update(stats_a.index)
    if not stats_b.empty:
        all_families.update(stats_b.index)

    rows = []
    for fam in sorted(all_families):
        count_a = int(stats_a.loc[fam, "hit_count"]) if fam in stats_a.index else 0
        count_b = int(stats_b.loc[fam, "hit_count"]) if fam in stats_b.index else 0
        freq_a = float(stats_a.loc[fam, "frequency"]) if fam in stats_a.index else 0.0
        freq_b = float(stats_b.loc[fam, "frequency"]) if fam in stats_b.index else 0.0

        if freq_b > 0:
            fold = freq_a / freq_b
        elif freq_a > 0:
            fold = float("inf")
        else:
            fold = 0.0

        if fold == float("inf"):
            log2 = float("inf")
        elif fold > 0:
            log2 = math.log2(fold)
        else:
            log2 = float("-inf")

        rows.append({
            "te_id": fam,
            "count_a": count_a,
            "count_b": count_b,
            "freq_a": freq_a,
            "freq_b": freq_b,
            "fold_enrichment": fold,
            "log2_enrichment": log2,
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.set_index("te_id")
    return result
