"""Composable BLAST hit quality filtering.

Pure functions that accept a DataFrame and return a filtered copy.
Replaces ad-hoc filtering scattered across legacy analysis scripts.
"""

import pandas as pd


def filter_by_evalue(df: pd.DataFrame, max_evalue: float) -> pd.DataFrame:
    """Keep hits with e-value <= threshold."""
    return df[df["evalue"] <= max_evalue].copy()


def filter_by_pident(df: pd.DataFrame, min_pident: float) -> pd.DataFrame:
    """Keep hits with percent identity >= threshold."""
    return df[df["pident"] >= min_pident].copy()


def filter_by_length(df: pd.DataFrame, min_length: int) -> pd.DataFrame:
    """Keep hits with alignment length >= threshold."""
    return df[df["length"] >= min_length].copy()


def apply_filters(
    df: pd.DataFrame,
    max_evalue: float | None = None,
    min_pident: float | None = None,
    min_length: int | None = None,
) -> pd.DataFrame:
    """Apply multiple quality filters in sequence.

    Args:
        df: BLAST results DataFrame.
        max_evalue: Maximum e-value (inclusive).
        min_pident: Minimum percent identity (inclusive).
        min_length: Minimum alignment length (inclusive).

    Returns:
        Filtered copy of the DataFrame.
    """
    result = df.copy()

    if max_evalue is not None:
        result = result[result["evalue"] <= max_evalue]
    if min_pident is not None:
        result = result[result["pident"] >= min_pident]
    if min_length is not None:
        result = result[result["length"] >= min_length]

    return result
