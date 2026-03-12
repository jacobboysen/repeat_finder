"""Tests for positional profiling of TE fossil hits."""

import numpy as np
import pandas as pd
import pytest

from fossil_finder.analysis.positional import (
    compute_end_bias,
    compute_positional_profile,
    compute_te_position,
    compute_utr_position,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_blast_df(**overrides):
    """Return a small BLAST-like DataFrame with sensible defaults."""
    defaults = {
        "qseqid": ["tr1", "tr2", "tr3"],
        "sseqid": ["TE_gypsy", "TE_copia", "TE_jockey"],
        "qstart": [100, 200, 50],
        "qend": [200, 400, 150],
        "qlen": [1000, 1000, 500],
        "sstart": [50, 100, 25],
        "send": [150, 300, 75],
        "slen": [500, 600, 300],
        "strand": ["plus", "plus", "minus"],
    }
    defaults.update(overrides)
    return pd.DataFrame(defaults)


# ---------------------------------------------------------------------------
# TestComputeUtrPosition
# ---------------------------------------------------------------------------

class TestComputeUtrPosition:
    def test_basic_calculation(self):
        df = _make_blast_df()
        result = compute_utr_position(df)
        assert "normalized_utr_pos" in result.columns
        expected = [100 / 1000, 200 / 1000, 50 / 500]
        np.testing.assert_allclose(result["normalized_utr_pos"].values, expected)

    def test_does_not_modify_input(self):
        df = _make_blast_df()
        original_cols = set(df.columns)
        compute_utr_position(df)
        assert set(df.columns) == original_cols

    def test_qlen_zero_gives_nan(self):
        df = _make_blast_df(qlen=[0, 1000, 500])
        result = compute_utr_position(df)
        assert np.isnan(result["normalized_utr_pos"].iloc[0])
        assert not np.isnan(result["normalized_utr_pos"].iloc[1])

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["qseqid", "qstart", "qend", "qlen"])
        result = compute_utr_position(df)
        assert "normalized_utr_pos" in result.columns
        assert len(result) == 0


# ---------------------------------------------------------------------------
# TestComputeTePosition
# ---------------------------------------------------------------------------

class TestComputeTePosition:
    def test_plus_strand(self):
        df = _make_blast_df(strand=["plus", "plus", "plus"])
        result = compute_te_position(df)
        expected = [50 / 500, 100 / 600, 25 / 300]
        np.testing.assert_allclose(result["normalized_te_pos"].values, expected)

    def test_minus_strand(self):
        df = _make_blast_df(strand=["minus", "minus", "minus"])
        result = compute_te_position(df)
        # For minus: (slen - sstart) / slen
        expected = [
            (500 - 50) / 500,
            (600 - 100) / 600,
            (300 - 25) / 300,
        ]
        np.testing.assert_allclose(result["normalized_te_pos"].values, expected)

    def test_mixed_strands(self):
        df = _make_blast_df(strand=["plus", "minus", "plus"])
        result = compute_te_position(df)
        expected = [
            50 / 500,           # plus
            (600 - 100) / 600,  # minus
            25 / 300,           # plus
        ]
        np.testing.assert_allclose(result["normalized_te_pos"].values, expected)

    def test_slen_zero_gives_nan(self):
        df = _make_blast_df(slen=[0, 600, 300])
        result = compute_te_position(df)
        assert np.isnan(result["normalized_te_pos"].iloc[0])
        assert not np.isnan(result["normalized_te_pos"].iloc[1])

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["sstart", "send", "slen", "strand"])
        result = compute_te_position(df)
        assert "normalized_te_pos" in result.columns
        assert len(result) == 0

    def test_does_not_modify_input(self):
        df = _make_blast_df()
        original_cols = set(df.columns)
        compute_te_position(df)
        assert set(df.columns) == original_cols


# ---------------------------------------------------------------------------
# TestComputePositionalProfile
# ---------------------------------------------------------------------------

class TestComputePositionalProfile:
    def test_correct_bin_count(self):
        df = pd.DataFrame({"pos": [0.05, 0.15, 0.25, 0.55, 0.95]})
        profile = compute_positional_profile(df, "pos", n_bins=10)
        assert len(profile) == 10
        assert list(profile.columns) == ["bin_start", "bin_end", "n_hits", "pct_hits"]

    def test_hits_sum_to_total(self):
        df = pd.DataFrame({"pos": [0.1, 0.3, 0.5, 0.7, 0.9]})
        profile = compute_positional_profile(df, "pos", n_bins=5)
        assert profile["n_hits"].sum() == 5
        assert profile["pct_hits"].sum() == pytest.approx(100.0)

    def test_all_hits_in_one_bin(self):
        df = pd.DataFrame({"pos": [0.05, 0.05, 0.05]})
        profile = compute_positional_profile(df, "pos", n_bins=10)
        assert profile.iloc[0]["n_hits"] == 3
        assert profile.iloc[0]["pct_hits"] == pytest.approx(100.0)
        assert profile.iloc[1:]["n_hits"].sum() == 0

    def test_custom_n_bins(self):
        df = pd.DataFrame({"pos": [0.25, 0.75]})
        profile = compute_positional_profile(df, "pos", n_bins=4)
        assert len(profile) == 4
        assert profile["n_hits"].sum() == 2

    def test_empty_column(self):
        df = pd.DataFrame({"pos": pd.Series(dtype=float)})
        profile = compute_positional_profile(df, "pos", n_bins=5)
        assert len(profile) == 5
        assert profile["n_hits"].sum() == 0
        assert (profile["pct_hits"] == 0.0).all()

    def test_nan_values_excluded(self):
        df = pd.DataFrame({"pos": [0.1, np.nan, 0.9, np.nan]})
        profile = compute_positional_profile(df, "pos", n_bins=10)
        assert profile["n_hits"].sum() == 2


# ---------------------------------------------------------------------------
# TestComputeEndBias
# ---------------------------------------------------------------------------

class TestComputeEndBias:
    def test_balanced_hits(self):
        # Even distribution across entire UTR
        pos = np.linspace(0.0, 1.0, 100, endpoint=False)
        df = pd.DataFrame({"normalized_utr_pos": pos})
        result = compute_end_bias(df)
        assert result["n_total"] == 100
        assert result["five_prime_pct"] == pytest.approx(
            result["three_prime_pct"], abs=5.0
        )
        assert result["end_ratio"] == pytest.approx(1.0, abs=0.3)

    def test_three_prime_biased(self):
        # All hits at 3' end
        df = pd.DataFrame({"normalized_utr_pos": [0.85, 0.90, 0.95, 0.99]})
        result = compute_end_bias(df)
        assert result["three_prime_pct"] == pytest.approx(100.0)
        assert result["five_prime_pct"] == pytest.approx(0.0)
        assert result["end_ratio"] == float("inf")
        assert result["n_total"] == 4

    def test_five_prime_biased(self):
        # All hits at 5' end
        df = pd.DataFrame({"normalized_utr_pos": [0.01, 0.05, 0.10, 0.15]})
        result = compute_end_bias(df)
        assert result["five_prime_pct"] == pytest.approx(100.0)
        assert result["three_prime_pct"] == pytest.approx(0.0)
        assert result["end_ratio"] == pytest.approx(0.0)

    def test_empty_data(self):
        df = pd.DataFrame({"normalized_utr_pos": pd.Series(dtype=float)})
        result = compute_end_bias(df)
        assert result["n_total"] == 0
        assert result["five_prime_pct"] == 0.0
        assert result["three_prime_pct"] == 0.0
        assert result["end_ratio"] == 0.0

    def test_nan_values_excluded(self):
        df = pd.DataFrame({"normalized_utr_pos": [0.1, np.nan, 0.9, np.nan]})
        result = compute_end_bias(df)
        assert result["n_total"] == 2
        assert result["five_prime_pct"] == pytest.approx(50.0)
        assert result["three_prime_pct"] == pytest.approx(50.0)

    def test_required_keys(self):
        df = pd.DataFrame({"normalized_utr_pos": [0.5]})
        result = compute_end_bias(df)
        assert set(result.keys()) == {
            "five_prime_pct", "three_prime_pct", "end_ratio", "n_total"
        }
