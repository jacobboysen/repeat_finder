"""Low-level sequence extraction and manipulation primitives.

Provides coordinate conversion (GFF3 1-based inclusive -> Python 0-based
half-open), reverse complement without BioPython, and MD5-based sequence
deduplication.
"""

import hashlib
from collections import defaultdict

_COMPLEMENT = str.maketrans("ATCGatcgNn", "TAGCtagcNn")


def gff_to_python_coords(gff_start: int, gff_end: int) -> tuple[int, int]:
    """Convert GFF3 coordinates to Python slice coordinates.

    GFF3 uses 1-based inclusive coordinates: feature at bases 1-10 means
    positions 1,2,3,...,10. Python slicing is 0-based half-open: [0:10].

    Args:
        gff_start: GFF3 start (1-based inclusive).
        gff_end: GFF3 end (1-based inclusive).

    Returns:
        (python_start, python_end) for use in seq[start:end].
    """
    return (gff_start - 1, gff_end)


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence.

    Handles uppercase, lowercase, and N bases. No BioPython dependency.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def extract_subsequence(
    genome: dict[str, str],
    chrom: str,
    gff_start: int,
    gff_end: int,
    strand: str,
) -> str:
    """Extract a subsequence from a loaded genome.

    Args:
        genome: Dict mapping chromosome names to sequences.
        chrom: Chromosome name (must exist in genome).
        gff_start: GFF3 start coordinate (1-based inclusive).
        gff_end: GFF3 end coordinate (1-based inclusive).
        strand: '+', '-', or '.' (dot treated as plus).

    Returns:
        Extracted sequence (reverse-complemented if minus strand).

    Raises:
        KeyError: If chrom not found in genome.
    """
    chrom_seq = genome[chrom]
    py_start, py_end = gff_to_python_coords(gff_start, gff_end)

    # Clamp to chromosome bounds
    py_start = max(0, py_start)
    py_end = min(len(chrom_seq), py_end)

    seq = chrom_seq[py_start:py_end]

    if strand == "-":
        seq = reverse_complement(seq)

    return seq


class SequenceDeduplicator:
    """MD5-based sequence deduplication.

    Tracks unique sequences and records which IDs map to the same
    sequence. Used to collapse identical isoforms/overlapping regions.
    """

    def __init__(self):
        self._hashes: dict[str, dict] = {}  # hash -> first region dict
        self._id_map: dict[str, list[str]] = defaultdict(list)  # hash -> all region IDs
        self.duplicates_skipped: int = 0

    def add(self, region_id: str, sequence: str, metadata: dict | None = None) -> bool:
        """Add a sequence. Returns True if unique, False if duplicate.

        Args:
            region_id: Unique identifier for this region.
            sequence: The DNA sequence.
            metadata: Optional metadata dict to store with the first occurrence.

        Returns:
            True if this is a new unique sequence, False if duplicate.
        """
        seq_hash = hashlib.md5(sequence.upper().encode()).hexdigest()
        self._id_map[seq_hash].append(region_id)

        if seq_hash in self._hashes:
            self.duplicates_skipped += 1
            return False

        self._hashes[seq_hash] = {
            "region_id": region_id,
            "sequence": sequence,
            "metadata": metadata or {},
        }
        return True

    def get_unique(self) -> list[dict]:
        """Return list of unique region dicts (region_id, sequence, metadata)."""
        return list(self._hashes.values())

    def get_isoform_map(self) -> dict[str, list[str]]:
        """Return mapping of primary region_id -> all duplicate IDs."""
        result = {}
        for seq_hash, entry in self._hashes.items():
            result[entry["region_id"]] = self._id_map[seq_hash]
        return result

    @property
    def stats(self) -> dict:
        return {
            "unique_sequences": len(self._hashes),
            "duplicates_skipped": self.duplicates_skipped,
            "total_input": sum(len(v) for v in self._id_map.values()),
        }
