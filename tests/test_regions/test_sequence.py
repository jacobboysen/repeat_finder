"""Tests for sequence extraction primitives."""

import pytest

from fossil_finder.regions.sequence import (
    SequenceDeduplicator,
    extract_subsequence,
    gff_to_python_coords,
    reverse_complement,
)


class TestGffToPythonCoords:
    """GFF3 uses 1-based inclusive coordinates. Python uses 0-based half-open."""

    def test_simple_conversion(self):
        start, end = gff_to_python_coords(1, 10)
        assert start == 0
        assert end == 10

    def test_internal_feature(self):
        start, end = gff_to_python_coords(100, 200)
        assert start == 99
        assert end == 200

    def test_single_base(self):
        """GFF feature spanning exactly one base."""
        start, end = gff_to_python_coords(5, 5)
        assert start == 4
        assert end == 5


class TestReverseComplement:
    def test_simple_sequence(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        assert reverse_complement("ATAT") == "ATAT"

    def test_lowercase_input(self):
        assert reverse_complement("atcg") == "cgat"

    def test_n_bases_preserved(self):
        assert reverse_complement("ANCG") == "CGNT"

    def test_empty_string(self):
        assert reverse_complement("") == ""


class TestExtractSubsequence:
    def test_plus_strand(self):
        genome = {"chr1": "AAACCCGGGTTT"}
        # GFF coords 4-9 -> Python [3:9] -> "CCCGGG"
        seq = extract_subsequence(genome, "chr1", 4, 9, "+")
        assert seq == "CCCGGG"

    def test_minus_strand_nonpalindrome(self):
        genome = {"chr1": "XXXXXATCXXXXX"}
        seq_plus = extract_subsequence(genome, "chr1", 6, 8, "+")
        assert seq_plus == "ATC"
        seq_minus = extract_subsequence(genome, "chr1", 6, 8, "-")
        assert seq_minus == "GAT"

    def test_unknown_chrom_raises(self):
        genome = {"chr1": "ATCG"}
        with pytest.raises(KeyError):
            extract_subsequence(genome, "chrZ", 1, 4, "+")

    def test_out_of_bounds_clamps(self):
        """Coordinates beyond chromosome end should be clamped, not crash."""
        genome = {"chr1": "ATCG"}
        seq = extract_subsequence(genome, "chr1", 1, 10, "+")
        assert seq == "ATCG"

    def test_dot_strand_treated_as_plus(self):
        """GFF '.' strand should be treated as plus."""
        genome = {"chr1": "ATCGATCG"}
        seq = extract_subsequence(genome, "chr1", 1, 4, ".")
        assert seq == "ATCG"


class TestSequenceDeduplicator:
    def test_unique_sequences_kept(self):
        dedup = SequenceDeduplicator()
        assert dedup.add("r1", "ATCG") is True
        assert dedup.add("r2", "GGCC") is True
        assert len(dedup.get_unique()) == 2

    def test_duplicate_rejected(self):
        dedup = SequenceDeduplicator()
        assert dedup.add("r1", "ATCG") is True
        assert dedup.add("r2", "ATCG") is False
        assert len(dedup.get_unique()) == 1

    def test_case_insensitive_dedup(self):
        dedup = SequenceDeduplicator()
        assert dedup.add("r1", "ATCG") is True
        assert dedup.add("r2", "atcg") is False

    def test_stats_tracking(self):
        dedup = SequenceDeduplicator()
        dedup.add("r1", "AAAA")
        dedup.add("r2", "CCCC")
        dedup.add("r3", "AAAA")  # duplicate
        dedup.add("r4", "CCCC")  # duplicate
        assert dedup.stats["unique_sequences"] == 2
        assert dedup.stats["duplicates_skipped"] == 2
        assert dedup.stats["total_input"] == 4

    def test_isoform_map(self):
        dedup = SequenceDeduplicator()
        dedup.add("iso1", "ATCG")
        dedup.add("iso2", "ATCG")
        dedup.add("iso3", "ATCG")
        imap = dedup.get_isoform_map()
        assert "iso1" in imap
        assert set(imap["iso1"]) == {"iso1", "iso2", "iso3"}

    def test_metadata_stored_for_first_occurrence(self):
        dedup = SequenceDeduplicator()
        dedup.add("r1", "ATCG", metadata={"gene": "geneA"})
        dedup.add("r2", "ATCG", metadata={"gene": "geneB"})
        unique = dedup.get_unique()
        assert unique[0]["metadata"]["gene"] == "geneA"

    def test_empty_deduplicator(self):
        dedup = SequenceDeduplicator()
        assert dedup.get_unique() == []
        assert dedup.stats["unique_sequences"] == 0
