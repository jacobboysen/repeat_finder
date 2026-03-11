"""Tests for BLAST hit filtering."""

import pandas as pd
import pytest

from fossil_finder.blast.filter import apply_filters, filter_by_evalue, filter_by_pident, filter_by_length


@pytest.fixture
def sample_hits():
    """DataFrame with varied quality hits for filter testing."""
    return pd.DataFrame({
        "qseqid": ["q1", "q2", "q3", "q4", "q5"],
        "sseqid": ["s1", "s1", "s2", "s2", "s3"],
        "pident": [95.0, 65.0, 80.0, 45.0, 90.0],
        "length": [200, 50, 100, 30, 150],
        "evalue": [1e-20, 0.5, 1e-5, 10.0, 1e-10],
        "bitscore": [150.0, 20.0, 80.0, 10.0, 120.0],
        "strand": ["plus", "minus", "plus", "plus", "minus"],
    })


class TestFilterByEvalue:
    def test_strict_evalue(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=1e-5)
        assert len(result) == 3  # q1, q3, q5
        assert all(result["evalue"] <= 1e-5)

    def test_permissive_evalue(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=10.0)
        assert len(result) == 5  # all pass

    def test_no_hits_pass(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=1e-30)
        assert len(result) == 0


class TestFilterByPident:
    def test_high_identity(self, sample_hits):
        result = filter_by_pident(sample_hits, min_pident=80.0)
        assert len(result) == 3  # q1, q3, q5

    def test_low_threshold_passes_all(self, sample_hits):
        result = filter_by_pident(sample_hits, min_pident=0.0)
        assert len(result) == 5


class TestFilterByLength:
    def test_min_length(self, sample_hits):
        result = filter_by_length(sample_hits, min_length=100)
        assert len(result) == 3  # q1, q3, q5

    def test_length_zero_passes_all(self, sample_hits):
        result = filter_by_length(sample_hits, min_length=0)
        assert len(result) == 5


class TestApplyFilters:
    def test_combined_filters(self, sample_hits):
        result = apply_filters(
            sample_hits, max_evalue=1e-5, min_pident=80.0, min_length=100,
        )
        assert len(result) == 3  # q1, q3, q5

    def test_no_filters_returns_copy(self, sample_hits):
        result = apply_filters(sample_hits)
        assert len(result) == 5
        assert result is not sample_hits  # returns copy, not same object

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["qseqid", "evalue", "pident", "length"])
        result = apply_filters(df, max_evalue=1e-5)
        assert len(result) == 0

    def test_filters_are_composable(self, sample_hits):
        """Applying filters sequentially equals applying all at once."""
        step1 = filter_by_evalue(sample_hits, max_evalue=0.5)
        step2 = filter_by_pident(step1, min_pident=80.0)

        combined = apply_filters(sample_hits, max_evalue=0.5, min_pident=80.0)
        assert list(combined["qseqid"]) == list(step2["qseqid"])
