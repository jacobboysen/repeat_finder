"""Genomic region extraction from GFF3 annotations + genome FASTA."""

from .extractor import RegionExtractor
from .sequence import (
    SequenceDeduplicator,
    extract_subsequence,
    gff_to_python_coords,
    reverse_complement,
)

__all__ = [
    "RegionExtractor",
    "SequenceDeduplicator",
    "extract_subsequence",
    "gff_to_python_coords",
    "reverse_complement",
]
