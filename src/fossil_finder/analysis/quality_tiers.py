"""Quality tier assignment and edit statistics for BLAST hits.

Ports the quality tier and edit stats logic from scripts/run_full_analysis.py
into the fossil_finder package.  All functions are pure, stateless, and operate
on pandas DataFrames using vectorized operations.

Quality tiers classify BLAST hits into three categories based on percent
identity and alignment length:

- **strict**: high-confidence matches (default: pident >= 85, length >= 100)
- **moderate**: mid-confidence matches (default: pident >= 75, length >= 50)
- **relaxed**: everything else
"""

import pandas as pd


def assign_quality_tiers(
    df: pd.DataFrame,
    strict_pident: float = 85,
    strict_length: int = 100,
    moderate_pident: float = 75,
    moderate_length: int = 50,
) -> pd.DataFrame:
    """Add a ``tier`` column classifying each hit as strict/moderate/relaxed.

    Assignment rules (applied in order of decreasing stringency):

    - **strict**: ``pident >= strict_pident`` AND ``length >= strict_length``
    - **moderate**: ``pident >= moderate_pident`` AND ``length >= moderate_length``
      (and not already strict)
    - **relaxed**: everything else

    Args:
        df: BLAST results with at least ``pident`` and ``length`` columns.
        strict_pident: Minimum percent identity for strict tier.
        strict_length: Minimum alignment length for strict tier.
        moderate_pident: Minimum percent identity for moderate tier.
        moderate_length: Minimum alignment length for moderate tier.

    Returns:
        Copy of *df* with an added ``tier`` column.
    """
    if df.empty:
        result = df.copy()
        result["tier"] = pd.Series(dtype=str)
        return result

    result = df.copy()
    strict_mask = (result["pident"] >= strict_pident) & (result["length"] >= strict_length)
    moderate_mask = (result["pident"] >= moderate_pident) & (result["length"] >= moderate_length)

    result["tier"] = "relaxed"
    result.loc[moderate_mask, "tier"] = "moderate"
    result.loc[strict_mask, "tier"] = "strict"
    return result


def compute_edit_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-hit edit metrics from BLAST mismatch and gapopen columns.

    Adds three columns:

    - ``mismatch_rate``: ``mismatch / length`` (alignment-length-normalized)
    - ``gap_rate``: ``gapopen / length``
    - ``edit_distance``: ``mismatch + gapopen``

    Division by zero is avoided by clipping ``length`` to a minimum of 1.

    Args:
        df: BLAST results with ``mismatch``, ``gapopen``, and ``length`` columns.

    Returns:
        Copy of *df* with the three new columns added.
    """
    if df.empty:
        result = df.copy()
        result["mismatch_rate"] = pd.Series(dtype=float)
        result["gap_rate"] = pd.Series(dtype=float)
        result["edit_distance"] = pd.Series(dtype=int)
        return result

    result = df.copy()
    safe_length = result["length"].clip(lower=1)
    result["mismatch_rate"] = result["mismatch"] / safe_length
    result["gap_rate"] = result["gapopen"] / safe_length
    result["edit_distance"] = result["mismatch"] + result["gapopen"]
    return result


def summarize_tiers(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize hit statistics grouped by quality tier.

    Args:
        df: BLAST results with a ``tier`` column and at least
            ``pident``, ``length``, ``evalue`` columns.

    Returns:
        DataFrame indexed by tier with columns:
        ``n_hits``, ``pct``, ``mean_pident``, ``mean_length``, ``mean_evalue``.
    """
    columns = ["n_hits", "pct", "mean_pident", "mean_length", "mean_evalue"]

    if df.empty:
        return pd.DataFrame(columns=columns)

    total = len(df)
    grouped = df.groupby("tier")

    summary = pd.DataFrame({
        "n_hits": grouped.size(),
        "pct": grouped.size() / total * 100,
        "mean_pident": grouped["pident"].mean(),
        "mean_length": grouped["length"].mean(),
        "mean_evalue": grouped["evalue"].mean(),
    })

    return summary


def compute_tier_edit_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize edit rates grouped by quality tier.

    Requires that :func:`compute_edit_stats` has already been called on *df*
    (i.e., ``mismatch_rate`` and ``gap_rate`` columns must be present).

    Args:
        df: BLAST results with ``tier``, ``mismatch_rate``, and ``gap_rate``
            columns.

    Returns:
        DataFrame indexed by tier with columns:
        ``mean_mismatch_rate``, ``mean_gap_rate``.
    """
    columns = ["mean_mismatch_rate", "mean_gap_rate"]

    if df.empty:
        return pd.DataFrame(columns=columns)

    grouped = df.groupby("tier")

    summary = pd.DataFrame({
        "mean_mismatch_rate": grouped["mismatch_rate"].mean(),
        "mean_gap_rate": grouped["gap_rate"].mean(),
    })

    return summary
