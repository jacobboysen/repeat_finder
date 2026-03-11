"""BLAST I/O utilities for fossil_finder.

Provides consistent parsing of BLAST tabular output.
All TSV files use the 17-column format (see BLAST_COLUMNS).
"""

from pathlib import Path

import pandas as pd


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq", "strand",
]

BLAST_COLUMNS_NO_STRAND = BLAST_COLUMNS[:-1]


def classify_strand(sstart: int, send: int) -> str:
    """Classify hit strand based on subject coordinates.

    When sstart > send in BLAST output, the hit is on the minus strand.
    """
    return "plus" if sstart < send else "minus"


def load_blast_results(
    results_file: str | Path,
    add_strand: bool = True,
    min_length: int | None = None,
    min_pident: float | None = None,
    max_evalue: float | None = None,
) -> pd.DataFrame:
    """Load BLAST results from TSV file.

    Auto-detects 16 vs 17 column format. Adds strand classification
    if missing. Returns empty DataFrame (with correct columns) for
    missing or empty files.
    """
    results_file = Path(results_file)

    if not results_file.exists() or results_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    with open(results_file) as f:
        first_line = f.readline()
    num_cols = len(first_line.strip().split("\t"))

    if num_cols == 17:
        df = pd.read_csv(results_file, sep="\t", names=BLAST_COLUMNS)
    elif num_cols == 16:
        df = pd.read_csv(results_file, sep="\t", names=BLAST_COLUMNS_NO_STRAND)
        if add_strand:
            df["strand"] = df.apply(
                lambda row: classify_strand(row["sstart"], row["send"]), axis=1
            )
    else:
        basic_cols = BLAST_COLUMNS[: min(num_cols, len(BLAST_COLUMNS))]
        df = pd.read_csv(results_file, sep="\t", names=basic_cols)
        if add_strand and "sstart" in df.columns and "send" in df.columns:
            df["strand"] = df.apply(
                lambda row: classify_strand(row["sstart"], row["send"]), axis=1
            )

    if min_length is not None and "length" in df.columns:
        df = df[df["length"] >= min_length]
    if min_pident is not None and "pident" in df.columns:
        df = df[df["pident"] >= min_pident]
    if max_evalue is not None and "evalue" in df.columns:
        df = df[df["evalue"] <= max_evalue]

    return df


def parse_blast_line(line: str) -> dict:
    """Parse a single BLAST TSV line into a dictionary.

    Useful for streaming large files without loading into memory.
    """
    parts = line.strip().split("\t")
    result = {}

    for i, col in enumerate(BLAST_COLUMNS):
        if i >= len(parts):
            result[col] = None
            continue

        value = parts[i]
        if col in ("pident", "evalue", "bitscore"):
            result[col] = float(value) if value else 0.0
        elif col in (
            "length", "mismatch", "gapopen", "qstart", "qend",
            "sstart", "send", "qlen", "slen",
        ):
            result[col] = int(value) if value else 0
        else:
            result[col] = value

    if "strand" not in result or result["strand"] is None:
        if result.get("sstart") is not None and result.get("send") is not None:
            result["strand"] = classify_strand(result["sstart"], result["send"])

    return result


def iter_blast_results(results_file: str | Path):
    """Iterate over BLAST results line by line.

    Yields dict per hit. Useful for large files that don't fit in memory.
    Raises FileNotFoundError if file doesn't exist (fail-fast).
    """
    results_file = Path(results_file)
    if not results_file.exists():
        raise FileNotFoundError(f"BLAST results file not found: {results_file}")

    with open(results_file) as f:
        for line in f:
            if line.strip():
                yield parse_blast_line(line)


def summarize_blast_results(df: pd.DataFrame) -> dict:
    """Generate summary statistics for BLAST results DataFrame."""
    if df.empty:
        return {
            "total_hits": 0, "unique_queries": 0, "unique_subjects": 0,
            "mean_pident": 0, "mean_length": 0, "mean_evalue": 0,
            "strand_plus": 0, "strand_minus": 0,
        }

    summary = {
        "total_hits": len(df),
        "unique_queries": df["qseqid"].nunique() if "qseqid" in df.columns else 0,
        "unique_subjects": df["sseqid"].nunique() if "sseqid" in df.columns else 0,
        "mean_pident": df["pident"].mean() if "pident" in df.columns else 0,
        "mean_length": df["length"].mean() if "length" in df.columns else 0,
        "mean_evalue": df["evalue"].mean() if "evalue" in df.columns else 0,
    }

    if "strand" in df.columns:
        strand_counts = df["strand"].value_counts()
        summary["strand_plus"] = strand_counts.get("plus", 0)
        summary["strand_minus"] = strand_counts.get("minus", 0)
    else:
        summary["strand_plus"] = 0
        summary["strand_minus"] = 0

    return summary
