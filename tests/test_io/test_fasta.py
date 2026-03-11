"""Tests for FASTA I/O module."""

import pytest

from fossil_finder.io.fasta import (
    parse_fasta,
    parse_fasta_headers,
    iter_fasta,
    write_fasta,
)


class TestParseFasta:
    def test_parse_sequences(self, mini_genome_fasta):
        seqs = parse_fasta(mini_genome_fasta)
        assert len(seqs) == 2
        assert "chr1" in seqs
        assert "chr2" in seqs

    def test_sequences_are_strings(self, mini_genome_fasta):
        seqs = parse_fasta(mini_genome_fasta)
        for seq in seqs.values():
            assert isinstance(seq, str)
            assert "\n" not in seq

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_fasta(tmp_path / "nope.fasta")

    def test_parse_te_fasta(self, mini_tes_fasta):
        seqs = parse_fasta(mini_tes_fasta)
        assert len(seqs) == 3
        assert "TE_LTR1" in seqs
        assert "TE_DNA1" in seqs
        assert "TE_LINE1" in seqs


class TestParseFastaHeaders:
    def test_extract_key_value_attributes(self, mini_tes_fasta):
        headers = parse_fasta_headers(mini_tes_fasta)
        assert headers["TE_LTR1"]["class"] == "LTR"
        assert headers["TE_LTR1"]["family"] == "gypsy"

    def test_raw_description_preserved(self, mini_tes_fasta):
        headers = parse_fasta_headers(mini_tes_fasta)
        assert "description" in headers["TE_LTR1"]

    def test_plain_headers(self, mini_genome_fasta):
        headers = parse_fasta_headers(mini_genome_fasta)
        assert headers["chr1"]["description"] == "Synthetic chromosome 1"


class TestIterFasta:
    def test_yields_id_and_sequence(self, mini_genome_fasta):
        records = list(iter_fasta(mini_genome_fasta))
        assert len(records) == 2
        seq_id, seq = records[0]
        assert seq_id == "chr1"
        assert isinstance(seq, str)


class TestWriteFasta:
    def test_roundtrip(self, mini_tes_fasta, tmp_path):
        seqs = parse_fasta(mini_tes_fasta)
        out_path = tmp_path / "out.fasta"
        write_fasta(seqs, out_path)
        reloaded = parse_fasta(out_path)
        assert seqs == reloaded

    def test_line_wrapping(self, tmp_path):
        seqs = {"test": "A" * 200}
        out_path = tmp_path / "wrapped.fasta"
        write_fasta(seqs, out_path, line_width=80)
        text = out_path.read_text()
        lines = text.strip().split("\n")
        assert lines[0] == ">test"
        assert len(lines[1]) == 80
