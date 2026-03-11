"""Generic FASTA I/O for fossil_finder.

Parses any FASTA file without organism-specific assumptions.
Handles key=value header attributes flexibly.
"""

import re
from pathlib import Path
from collections.abc import Iterator


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Parse FASTA file into {sequence_id: sequence} dict.

    Sequence ID is the first whitespace-delimited token after '>'.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def parse_fasta_headers(path: str | Path) -> dict[str, dict]:
    """Parse FASTA headers, extracting key=value attributes.

    Returns {seq_id: {"description": "...", "key1": "val1", ...}}.
    Handles both 'key=value;' and 'key=value' formats.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    headers: dict[str, dict] = {}

    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            parts = line.strip()[1:].split(None, 1)
            seq_id = parts[0]
            raw_desc = parts[1] if len(parts) > 1 else ""

            attrs: dict[str, str] = {"description": raw_desc}

            for match in re.finditer(r"(\w+)=([^;\s]+)", raw_desc):
                attrs[match.group(1)] = match.group(2)

            headers[seq_id] = attrs

    return headers


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Iterate over FASTA records yielding (seq_id, sequence) tuples.

    Memory-efficient for large files.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    current_id: str | None = None
    current_seq: list[str] = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            yield current_id, "".join(current_seq)


def write_fasta(
    sequences: dict[str, str],
    path: str | Path,
    line_width: int = 80,
) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: {seq_id: sequence} dict.
        path: Output file path.
        line_width: Characters per sequence line (default 80, FASTA standard).
    """
    path = Path(path)
    with open(path, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")
