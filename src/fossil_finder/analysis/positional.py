"""Positional profiling of TE fossil hits within UTRs and TE elements.

Computes normalized hit positions along UTR queries and TE subjects,
bins them into positional profiles, and measures 3'/5' end bias.
All functions are pure, stateless, and operate on pandas DataFrames.
"""

import numpy as np
import pandas as pd


def compute_utr_position(df: pd.DataFrame) -> pd.DataFrame:
    """Add normalized UTR position column based on query start.

    Computes ``normalized_utr_pos = qstart / qlen``, giving the position
    of each hit within the UTR query on a 0-1 scale.

    Args:
        df: BLAST results with ``qstart`` and ``qlen`` columns.

    Returns:
        Copy of *df* with ``normalized_utr_pos`` column added.
        Rows where ``qlen == 0`` receive ``NaN``.
    """
    result = df.copy()
    if result.empty:
        result["normalized_utr_pos"] = pd.Series(dtype=float)
        return result

    qlen = result["qlen"].astype(float)
    result["normalized_utr_pos"] = np.where(
        qlen == 0,
        np.nan,
        result["qstart"].astype(float) / qlen,
    )
    return result


def compute_te_position(df: pd.DataFrame) -> pd.DataFrame:
    """Add normalized TE position column based on subject start and strand.

    For plus-strand hits: ``sstart / slen``.
    For minus-strand hits: ``(slen - sstart) / slen``.

    Args:
        df: BLAST results with ``sstart``, ``slen``, and ``strand`` columns.

    Returns:
        Copy of *df* with ``normalized_te_pos`` column added.
        Rows where ``slen == 0`` receive ``NaN``.
    """
    result = df.copy()
    if result.empty:
        result["normalized_te_pos"] = pd.Series(dtype=float)
        return result

    slen = result["slen"].astype(float)
    sstart = result["sstart"].astype(float)
    is_minus = result["strand"] == "minus"

    raw_pos = np.where(is_minus, slen - sstart, sstart)
    result["normalized_te_pos"] = np.where(
        slen == 0,
        np.nan,
        raw_pos / slen,
    )
    return result


def compute_positional_profile(
    df: pd.DataFrame,
    position_col: str,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Bin a normalized position column into equal-width bins.

    Args:
        df: DataFrame containing *position_col* with values in [0, 1].
        position_col: Name of the column to bin.
        n_bins: Number of equal-width bins (default 10).

    Returns:
        DataFrame with columns: ``bin_start``, ``bin_end``, ``n_hits``,
        ``pct_hits``.  One row per bin, ordered from 0 to 1.
    """
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)

    # Drop NaN values for binning
    valid = df[position_col].dropna()

    counts, _ = np.histogram(valid, bins=bin_edges)
    total = counts.sum()

    rows = []
    for i in range(n_bins):
        rows.append({
            "bin_start": round(bin_edges[i], 6),
            "bin_end": round(bin_edges[i + 1], 6),
            "n_hits": int(counts[i]),
            "pct_hits": (counts[i] / total * 100) if total > 0 else 0.0,
        })

    return pd.DataFrame(rows)


def compute_end_bias(df: pd.DataFrame) -> dict:
    """Compute 3'/5' end bias from normalized UTR positions.

    Measures the fraction of hits in the first 20% (5' end) vs the last
    20% (3' end) of the UTR.

    Args:
        df: DataFrame with ``normalized_utr_pos`` column (0-1 scale).

    Returns:
        Dict with keys:
        - ``five_prime_pct``: Percentage of hits in the first 20%.
        - ``three_prime_pct``: Percentage of hits in the last 20%.
        - ``end_ratio``: ``three_prime_pct / five_prime_pct``
          (``float('inf')`` if five_prime_pct is 0 and three_prime_pct > 0,
          ``0.0`` if both are 0).
        - ``n_total``: Number of non-NaN positions considered.
    """
    valid = df["normalized_utr_pos"].dropna()
    n_total = len(valid)

    if n_total == 0:
        return {
            "five_prime_pct": 0.0,
            "three_prime_pct": 0.0,
            "end_ratio": 0.0,
            "n_total": 0,
        }

    five_prime_n = (valid < 0.2).sum()
    three_prime_n = (valid >= 0.8).sum()

    five_prime_pct = five_prime_n / n_total * 100
    three_prime_pct = three_prime_n / n_total * 100

    if five_prime_pct == 0:
        end_ratio = float("inf") if three_prime_pct > 0 else 0.0
    else:
        end_ratio = three_prime_pct / five_prime_pct

    return {
        "five_prime_pct": five_prime_pct,
        "three_prime_pct": three_prime_pct,
        "end_ratio": end_ratio,
        "n_total": n_total,
    }
