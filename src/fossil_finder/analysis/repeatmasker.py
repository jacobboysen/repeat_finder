"""RepeatMasker .out file parsing and overlap classification.

Parses RepeatMasker output, detects overlaps with query regions,
and classifies BLAST hits as known (in RM annotation) vs novel.
"""

from pathlib import Path

import pandas as pd


def parse_repeatmasker_out(path: str | Path) -> pd.DataFrame:
    """Parse RepeatMasker .out file into a DataFrame.

    Args:
        path: Path to RepeatMasker .out file.

    Returns:
        DataFrame with columns: score, divergence, chrom, start, end,
        strand, repeat_name, repeat_class.

    Raises:
        FileNotFoundError: If file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"RepeatMasker file not found: {path}")

    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue

            parts = line.split()
            if len(parts) < 15:
                continue

            try:
                score = int(parts[0])
                divergence = float(parts[1])
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                strand_raw = parts[8]
                strand = "-" if strand_raw == "C" else "+"
                repeat_name = parts[9]
                repeat_class = parts[10]
            except (ValueError, IndexError):
                continue

            records.append({
                "score": score,
                "divergence": divergence,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "repeat_name": repeat_name,
                "repeat_class": repeat_class,
            })

    return pd.DataFrame(records)


def find_overlaps(
    rm_regions: pd.DataFrame,
    query_regions: pd.DataFrame,
) -> pd.DataFrame:
    """Find overlaps between RepeatMasker regions and query regions.

    Uses chromosome-grouped merge + vectorized interval intersection
    instead of O(N*M) nested loop. Both inputs use 1-based inclusive
    coordinates (GFF3/RM convention).

    Args:
        rm_regions: RepeatMasker regions with chrom, start, end,
            repeat_name, repeat_class, strand, divergence.
        query_regions: Query regions with region_id, chrom, start, end.

    Returns:
        DataFrame of overlaps with region_id, overlap_bp, and RM metadata.
    """
    if rm_regions.empty or query_regions.empty:
        return pd.DataFrame()

    # Inner join on chromosome to reduce comparison space
    merged = rm_regions.merge(
        query_regions,
        on="chrom",
        suffixes=("_rm", "_qr"),
    )

    if merged.empty:
        return pd.DataFrame()

    # Vectorized overlap detection
    overlap_start = merged[["start_rm", "start_qr"]].max(axis=1)
    overlap_end = merged[["end_rm", "end_qr"]].min(axis=1)
    has_overlap = overlap_start <= overlap_end

    merged = merged[has_overlap].copy()
    if merged.empty:
        return pd.DataFrame()

    overlap_start = merged[["start_rm", "start_qr"]].max(axis=1)
    overlap_end = merged[["end_rm", "end_qr"]].min(axis=1)

    merged["overlap_start"] = overlap_start
    merged["overlap_end"] = overlap_end
    merged["overlap_bp"] = overlap_end - overlap_start + 1

    # Compute RM interval in query-relative 1-based coordinates
    query_len = merged["end_qr"] - merged["start_qr"] + 1
    is_minus = merged.get("strand_qr", pd.Series("+" * len(merged))) == "-"

    # Plus strand: relative = overlap - query_start + 1
    rm_start_rel = overlap_start - merged["start_qr"] + 1
    rm_end_rel = overlap_end - merged["start_qr"] + 1

    # Minus strand: positions are reversed
    if is_minus.any():
        rm_start_rel_minus = merged["end_qr"] - overlap_end + 1
        rm_end_rel_minus = merged["end_qr"] - overlap_start + 1
        rm_start_rel = rm_start_rel.where(~is_minus, rm_start_rel_minus)
        rm_end_rel = rm_end_rel.where(~is_minus, rm_end_rel_minus)

    rm_start_rel = rm_start_rel.clip(lower=1)
    rm_end_rel = rm_end_rel.clip(upper=query_len)

    result = pd.DataFrame({
        "region_id": merged["region_id"].values,
        "chrom": merged["chrom"].values,
        "overlap_start": merged["overlap_start"].values,
        "overlap_end": merged["overlap_end"].values,
        "overlap_bp": merged["overlap_bp"].values,
        "repeat_name": merged["repeat_name"].values,
        "repeat_class": merged["repeat_class"].values,
        "divergence": merged["divergence"].values if "divergence" in merged.columns else None,
        "rm_start_in_query": rm_start_rel.values,
        "rm_end_in_query": rm_end_rel.values,
    })

    return result


def classify_hits(
    blast_hits: pd.DataFrame,
    rm_overlaps: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classify BLAST hits as known (RepeatMasker) or novel.

    A hit is "known" if its query range [qstart, qend] overlaps with
    any RepeatMasker region mapped to the same query.

    Uses vectorized merge + interval comparison instead of row-by-row iteration.

    Args:
        blast_hits: BLAST results with qseqid, qstart, qend columns.
        rm_overlaps: RM overlap data with region_id, rm_start_in_query,
            rm_end_in_query, repeat_name, repeat_class.

    Returns:
        Tuple of (known_hits, novel_hits) DataFrames.
    """
    if blast_hits.empty:
        return pd.DataFrame(), pd.DataFrame()

    if rm_overlaps.empty:
        return pd.DataFrame(), blast_hits.copy()

    # Add a unique hit index for tracking
    hits = blast_hits.copy()
    hits["_hit_idx"] = range(len(hits))

    # Merge BLAST hits with RM overlaps on query ID
    rm_sub = rm_overlaps[["region_id", "rm_start_in_query", "rm_end_in_query",
                           "repeat_name", "repeat_class"]].copy()

    merged = hits.merge(
        rm_sub,
        left_on="qseqid",
        right_on="region_id",
        how="inner",
    )

    if merged.empty:
        return pd.DataFrame(), blast_hits.copy()

    # Vectorized overlap check: hit [qstart, qend] overlaps RM [rm_start, rm_end]
    has_overlap = (
        (merged["qstart"] <= merged["rm_end_in_query"]) &
        (merged["qend"] >= merged["rm_start_in_query"])
    )

    # Get unique hit indices that are known (take first matching RM annotation)
    known_merged = merged[has_overlap].drop_duplicates(subset=["_hit_idx"], keep="first")
    known_idx = set(known_merged["_hit_idx"].values)

    # Build known DataFrame with RM metadata
    known = hits[hits["_hit_idx"].isin(known_idx)].copy()
    # Attach RM metadata via the merge result
    rm_meta = known_merged.set_index("_hit_idx")[["repeat_name", "repeat_class"]]
    rm_meta.columns = ["rm_repeat_name", "rm_repeat_class"]
    known = known.join(rm_meta, on="_hit_idx")

    novel = hits[~hits["_hit_idx"].isin(known_idx)].copy()

    # Clean up temp column
    known = known.drop(columns=["_hit_idx"])
    novel = novel.drop(columns=["_hit_idx"])

    return known, novel
