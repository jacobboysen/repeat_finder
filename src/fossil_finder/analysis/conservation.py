"""Conservation / phyloP scoring for TE fossil hits.

Converts BLAST hits to genomic BED coordinates, scores them against a
phyloP bigWig file using bigWigAverageOverBed, and provides summary
statistics (group-level aggregation, pident-conservation correlation).

Ported from WS6 logic in scripts/run_full_analysis.py.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats


def hits_to_genomic_bed(
    df: pd.DataFrame,
    regions: pd.DataFrame,
) -> pd.DataFrame:
    """Convert BLAST hits to genomic BED coordinates.

    Maps query-relative hit positions (qstart, qend) onto the genome
    using the region metadata (chrom, start, end, strand).

    Args:
        df: BLAST results with at least qseqid, qstart, qend columns.
        regions: Region metadata with region_id, chrom, start, end, strand.

    Returns:
        DataFrame with original columns plus bed_chrom, g_start, g_end,
        and hit_id.  Invalid coordinates (g_start < 0 or g_end <= g_start)
        are removed.
    """
    if df.empty:
        out = df.copy()
        for col in ("bed_chrom", "g_start", "g_end", "hit_id"):
            out[col] = pd.Series(dtype="object" if col in ("bed_chrom", "hit_id") else "int64")
        return out

    # Build a lookup keyed by region_id
    region_df = regions.set_index("region_id")[["chrom", "start", "end", "strand"]].rename(
        columns={"chrom": "r_chrom", "start": "r_start", "end": "r_end", "strand": "r_strand"},
    )

    hits = df.merge(region_df, left_on="qseqid", right_index=True, how="inner")

    if hits.empty:
        out = df.iloc[:0].copy()
        for col in ("bed_chrom", "g_start", "g_end", "hit_id"):
            out[col] = pd.Series(dtype="object" if col in ("bed_chrom", "hit_id") else "int64")
        return out

    # Vectorised coordinate conversion
    plus_mask = hits["r_strand"] == "+"

    # + strand: bed_start = region_start + qstart - 2  (0-based)
    #           bed_end   = region_start + qend   - 1  (exclusive)
    hits.loc[plus_mask, "g_start"] = (
        hits.loc[plus_mask, "r_start"] + hits.loc[plus_mask, "qstart"] - 2
    )
    hits.loc[plus_mask, "g_end"] = (
        hits.loc[plus_mask, "r_start"] + hits.loc[plus_mask, "qend"] - 1
    )

    # - strand: bed_start = region_end - qend      (0-based)
    #           bed_end   = region_end - qstart + 1 (exclusive)
    hits.loc[~plus_mask, "g_start"] = (
        hits.loc[~plus_mask, "r_end"] - hits.loc[~plus_mask, "qend"]
    )
    hits.loc[~plus_mask, "g_end"] = (
        hits.loc[~plus_mask, "r_end"] - hits.loc[~plus_mask, "qstart"] + 1
    )

    hits["g_start"] = hits["g_start"].astype(int)
    hits["g_end"] = hits["g_end"].astype(int)

    # UCSC BED requires "chr" prefix
    hits["bed_chrom"] = "chr" + hits["r_chrom"].astype(str)

    # Unique hit identifiers for round-trip merging
    hits["hit_id"] = [f"hit_{i}" for i in range(len(hits))]

    # Drop intermediate region columns
    hits = hits.drop(columns=["r_chrom", "r_start", "r_end", "r_strand"])

    # Filter invalid coordinates
    hits = hits[(hits["g_start"] >= 0) & (hits["g_end"] > hits["g_start"])]
    hits = hits.reset_index(drop=True)

    return hits


def score_with_bigwig(
    bed_df: pd.DataFrame,
    bigwig_path: str | Path,
    tool_path: str | Path,
) -> pd.DataFrame:
    """Score BED intervals with bigWigAverageOverBed.

    Writes a temporary BED file, invokes the external
    ``bigWigAverageOverBed`` binary, and parses its tab-delimited output.

    Args:
        bed_df: DataFrame with bed_chrom, g_start, g_end, hit_id columns.
        bigwig_path: Path to phyloP bigWig file.
        tool_path: Path to bigWigAverageOverBed binary.

    Returns:
        DataFrame with columns: hit_id, size, covered, sum, mean0, mean.

    Raises:
        FileNotFoundError: If *bigwig_path* or *tool_path* does not exist.
        RuntimeError: If the subprocess exits with a non-zero return code.
    """
    import subprocess

    bigwig_path = Path(bigwig_path)
    tool_path = Path(tool_path)

    if not bigwig_path.exists():
        raise FileNotFoundError(f"bigWig file not found: {bigwig_path}")
    if not tool_path.exists():
        raise FileNotFoundError(f"bigWigAverageOverBed not found: {tool_path}")

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False,
    ) as bed_fh, tempfile.NamedTemporaryFile(
        mode="r", suffix=".tab", delete=False,
    ) as out_fh:
        bed_path = Path(bed_fh.name)
        out_path = Path(out_fh.name)

    # Write sorted BED (no header)
    sorted_bed = bed_df[["bed_chrom", "g_start", "g_end", "hit_id"]].sort_values(
        ["bed_chrom", "g_start"],
    )
    sorted_bed.to_csv(bed_path, sep="\t", index=False, header=False)

    try:
        cmd = [str(tool_path), str(bigwig_path), str(bed_path), str(out_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            raise RuntimeError(
                f"bigWigAverageOverBed failed (rc={result.returncode}): "
                f"{result.stderr[:500]}"
            )

        scores = pd.read_csv(
            out_path,
            sep="\t",
            header=None,
            names=["hit_id", "size", "covered", "sum", "mean0", "mean"],
        )
    finally:
        bed_path.unlink(missing_ok=True)
        out_path.unlink(missing_ok=True)

    return scores


def summarize_conservation_by_group(
    scores_df: pd.DataFrame,
    group_col: str,
) -> pd.DataFrame:
    """Summarise phyloP scores by a categorical grouping column.

    Args:
        scores_df: DataFrame containing *group_col* and ``phylop_mean``.
        group_col: Column name to group by (e.g. ``"tier"``,
            ``"in_repeatmasker"``).

    Returns:
        DataFrame indexed by group value with columns:
        n_hits, mean_phylop, median_phylop, std_phylop,
        pct_positive, pct_negative.
    """
    if scores_df.empty or group_col not in scores_df.columns:
        return pd.DataFrame(
            columns=["n_hits", "mean_phylop", "median_phylop",
                      "std_phylop", "pct_positive", "pct_negative"],
        )

    rows: list[dict] = []
    for group_val, subset in scores_df.groupby(group_col):
        phylop = subset["phylop_mean"].dropna()
        if phylop.empty:
            continue
        rows.append({
            group_col: group_val,
            "n_hits": len(phylop),
            "mean_phylop": float(phylop.mean()),
            "median_phylop": float(phylop.median()),
            "std_phylop": float(phylop.std()),
            "pct_positive": float((phylop > 0).mean() * 100),
            "pct_negative": float((phylop < 0).mean() * 100),
        })

    if not rows:
        return pd.DataFrame(
            columns=["n_hits", "mean_phylop", "median_phylop",
                      "std_phylop", "pct_positive", "pct_negative"],
        )

    result = pd.DataFrame(rows).set_index(group_col)
    return result


def compute_pident_conservation_correlation(
    df: pd.DataFrame,
) -> Optional[dict]:
    """Spearman correlation between percent-identity and phyloP score.

    Args:
        df: DataFrame with ``pident`` and ``phylop_mean`` columns.

    Returns:
        Dict with rho, p_value, n — or ``None`` if fewer than 10 valid
        data points are available.
    """
    if df.empty:
        return None

    valid = df[["pident", "phylop_mean"]].dropna()
    if len(valid) < 10:
        return None

    rho, p_value = stats.spearmanr(valid["pident"], valid["phylop_mean"])
    return {
        "rho": float(rho),
        "p_value": float(p_value),
        "n": len(valid),
    }
