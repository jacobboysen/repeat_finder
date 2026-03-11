"""Tests for BLAST I/O module."""

import pandas as pd
import pytest

from fossil_finder.io.blast import (
    BLAST_COLUMNS,
    classify_strand,
    load_blast_results,
    parse_blast_line,
    iter_blast_results,
    summarize_blast_results,
)


class TestClassifyStrand:
    def test_plus_strand(self):
        assert classify_strand(1, 100) == "plus"

    def test_minus_strand(self):
        assert classify_strand(100, 1) == "minus"

    def test_equal_positions(self):
        assert classify_strand(50, 50) == "minus"


class TestBlastColumns:
    def test_column_count(self):
        assert len(BLAST_COLUMNS) == 17

    def test_required_columns_present(self):
        for col in ["qseqid", "sseqid", "pident", "evalue", "bitscore"]:
            assert col in BLAST_COLUMNS


class TestLoadBlastResults:
    def test_load_17col_file(self, mini_blast_results):
        df = load_blast_results(mini_blast_results)
        assert len(df) == 5
        assert list(df.columns) == BLAST_COLUMNS

    def test_filter_by_evalue(self, mini_blast_results):
        df = load_blast_results(mini_blast_results, max_evalue=1e-10)
        assert all(df["evalue"] <= 1e-10)

    def test_filter_by_min_length(self, mini_blast_results):
        df = load_blast_results(mini_blast_results, min_length=100)
        assert all(df["length"] >= 100)

    def test_nonexistent_file_returns_empty(self, tmp_path):
        df = load_blast_results(tmp_path / "nonexistent.tsv")
        assert df.empty
        assert list(df.columns) == BLAST_COLUMNS

    def test_empty_file_returns_empty(self, tmp_path):
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        df = load_blast_results(empty)
        assert df.empty


class TestParseBlastLine:
    def test_parse_complete_line(self):
        line = "q1\ts1\t80.0\t100\t15\t2\t1\t100\t1\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG\tplus"
        result = parse_blast_line(line)
        assert result["qseqid"] == "q1"
        assert result["pident"] == 80.0
        assert result["length"] == 100
        assert result["strand"] == "plus"

    def test_parse_16col_line_adds_strand(self):
        line = "q1\ts1\t80.0\t100\t15\t2\t1\t100\t200\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG"
        result = parse_blast_line(line)
        assert result["strand"] == "minus"


class TestIterBlastResults:
    def test_iterates_all_lines(self, mini_blast_results):
        hits = list(iter_blast_results(mini_blast_results))
        assert len(hits) == 5

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            list(iter_blast_results(tmp_path / "nope.tsv"))


class TestSummarizeBlastResults:
    def test_summary_stats(self, mini_blast_results):
        df = load_blast_results(mini_blast_results)
        summary = summarize_blast_results(df)
        assert summary["total_hits"] == 5
        assert summary["unique_queries"] == 3
        assert summary["unique_subjects"] == 4
        assert summary["strand_plus"] == 4
        assert summary["strand_minus"] == 1

    def test_empty_dataframe_summary(self):
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        summary = summarize_blast_results(df)
        assert summary["total_hits"] == 0
