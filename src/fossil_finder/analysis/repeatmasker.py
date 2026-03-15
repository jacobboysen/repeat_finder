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

    Uses simple interval intersection. Both inputs use 1-based inclusive
    coordinates (GFF3/RM convention).

    Args:
        rm_regions: RepeatMasker regions with chrom, start, end,
            repeat_name, repeat_class, strand, divergence.
        query_regions: Query regions with region_id, chrom, start, end.

    Returns:
        DataFrame of overlaps with region_id, overlap_bp, and RM metadata.
    """
    overlaps = []

    for _, rm in rm_regions.iterrows():
        for _, qr in query_regions.iterrows():
            if rm["chrom"] != qr["chrom"]:
                continue

            overlap_start = max(rm["start"], qr["start"])
            overlap_end = min(rm["end"], qr["end"])

            if overlap_start <= overlap_end:
                overlap_bp = overlap_end - overlap_start + 1
                # Compute RM interval in query-relative 1-based coordinates
                # to match BLAST qstart/qend (1-based inclusive)
                query_len = qr["end"] - qr["start"] + 1
                if "strand" in qr and qr["strand"] == "-":
                    # For minus-strand queries, positions are reversed
                    rm_start_rel = qr["end"] - overlap_end + 1
                    rm_end_rel = qr["end"] - overlap_start + 1
                else:
                    rm_start_rel = overlap_start - qr["start"] + 1
                    rm_end_rel = overlap_end - qr["start"] + 1
                rm_start_rel = max(1, rm_start_rel)
                rm_end_rel = min(query_len, rm_end_rel)
                overlaps.append({
                    "region_id": qr["region_id"],
                    "chrom": rm["chrom"],
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_bp": overlap_bp,
                    "repeat_name": rm["repeat_name"],
                    "repeat_class": rm["repeat_class"],
                    "divergence": rm.get("divergence", None),
                    "rm_start_in_query": rm_start_rel,
                    "rm_end_in_query": rm_end_rel,
                })

    return pd.DataFrame(overlaps)


def classify_hits(
    blast_hits: pd.DataFrame,
    rm_overlaps: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classify BLAST hits as known (RepeatMasker) or novel.

    A hit is "known" if its query range [qstart, qend] overlaps with
    any RepeatMasker region mapped to the same query.

    Args:
        blast_hits: BLAST results with qseqid, qstart, qend columns.
        rm_overlaps: RM overlap data with region_id, rm_start_in_query,
            rm_end_in_query, repeat_name, repeat_class.

    Returns:
        Tuple of (known_hits, novel_hits) DataFrames.
    """
    # Build lookup: region_id -> list of RM intervals in query space
    rm_by_region: dict[str, list[dict]] = {}
    for _, row in rm_overlaps.iterrows():
        rid = row["region_id"]
        if rid not in rm_by_region:
            rm_by_region[rid] = []
        rm_by_region[rid].append({
            "start": row["rm_start_in_query"],
            "end": row["rm_end_in_query"],
            "repeat_name": row["repeat_name"],
            "repeat_class": row["repeat_class"],
        })

    known_rows = []
    novel_rows = []

    for _, hit in blast_hits.iterrows():
        qid = hit["qseqid"]
        regions = rm_by_region.get(qid, [])

        is_known = False
        for rm in regions:
            if hit["qstart"] <= rm["end"] and hit["qend"] >= rm["start"]:
                known_row = hit.to_dict()
                known_row["rm_repeat_name"] = rm["repeat_name"]
                known_row["rm_repeat_class"] = rm["repeat_class"]
                known_rows.append(known_row)
                is_known = True
                break

        if not is_known:
            novel_rows.append(hit.to_dict())

    known = pd.DataFrame(known_rows) if known_rows else pd.DataFrame()
    novel = pd.DataFrame(novel_rows) if novel_rows else pd.DataFrame()

    return known, novel
