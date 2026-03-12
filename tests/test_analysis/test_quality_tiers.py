"""Tests for quality tier assignment and edit statistics."""

import pandas as pd
import pytest

from fossil_finder.analysis.quality_tiers import (
    assign_quality_tiers,
    compute_edit_stats,
    compute_tier_edit_summary,
    summarize_tiers,
)


@pytest.fixture
def sample_hits():
    """BLAST hits spanning all three quality tiers."""
    return pd.DataFrame({
        "pident": [95.0, 80.0, 60.0, 90.0, 70.0],
        "length": [200, 80, 30, 150, 60],
        "evalue": [1e-10, 1e-5, 0.5, 1e-8, 0.01],
        "mismatch": [10, 16, 12, 15, 18],
        "gapopen": [1, 2, 1, 2, 3],
    })


@pytest.fixture
def empty_df():
    """Empty DataFrame with the required BLAST columns."""
    return pd.DataFrame({
        "pident": pd.Series(dtype=float),
        "length": pd.Series(dtype=int),
        "evalue": pd.Series(dtype=float),
        "mismatch": pd.Series(dtype=int),
        "gapopen": pd.Series(dtype=int),
    })


# ── assign_quality_tiers ─────────────────────────────────────────────────────


class TestAssignQualityTiers:
    def test_strict_assignment(self, sample_hits):
        """Hits with pident >= 85 AND length >= 100 are strict."""
        result = assign_quality_tiers(sample_hits)
        # Row 0: pident=95, length=200 -> strict
        assert result.iloc[0]["tier"] == "strict"
        # Row 3: pident=90, length=150 -> strict
        assert result.iloc[3]["tier"] == "strict"

    def test_moderate_assignment(self, sample_hits):
        """Hits with pident >= 75 AND length >= 50 (but not strict) are moderate."""
        result = assign_quality_tiers(sample_hits)
        # Row 1: pident=80, length=80 -> moderate
        assert result.iloc[1]["tier"] == "moderate"
        # Row 4: pident=70, length=60 -> relaxed (pident < 75)
        assert result.iloc[4]["tier"] == "relaxed"

    def test_relaxed_assignment(self, sample_hits):
        """Hits not meeting moderate thresholds are relaxed."""
        result = assign_quality_tiers(sample_hits)
        # Row 2: pident=60, length=30 -> relaxed
        assert result.iloc[2]["tier"] == "relaxed"
        # Row 4: pident=70, length=60 -> relaxed (pident < 75)
        assert result.iloc[4]["tier"] == "relaxed"

    def test_custom_thresholds(self, sample_hits):
        """Custom thresholds override defaults."""
        result = assign_quality_tiers(
            sample_hits,
            strict_pident=90,
            strict_length=150,
            moderate_pident=70,
            moderate_length=50,
        )
        # Row 0: pident=95 >= 90, length=200 >= 150 -> strict
        assert result.iloc[0]["tier"] == "strict"
        # Row 3: pident=90 >= 90, length=150 >= 150 -> strict
        assert result.iloc[3]["tier"] == "strict"
        # Row 1: pident=80 >= 70, length=80 >= 50 -> moderate (not strict)
        assert result.iloc[1]["tier"] == "moderate"
        # Row 4: pident=70 >= 70, length=60 >= 50 -> moderate
        assert result.iloc[4]["tier"] == "moderate"
        # Row 2: pident=60 < 70 -> relaxed
        assert result.iloc[2]["tier"] == "relaxed"

    def test_boundary_exact_strict(self):
        """Hit exactly on strict boundary is classified as strict."""
        df = pd.DataFrame({
            "pident": [85.0],
            "length": [100],
        })
        result = assign_quality_tiers(df)
        assert result.iloc[0]["tier"] == "strict"

    def test_boundary_just_below_strict(self):
        """Hit just below strict boundary on pident falls to moderate."""
        df = pd.DataFrame({
            "pident": [84.9],
            "length": [100],
        })
        result = assign_quality_tiers(df)
        assert result.iloc[0]["tier"] == "moderate"

    def test_boundary_exact_moderate(self):
        """Hit exactly on moderate boundary is classified as moderate."""
        df = pd.DataFrame({
            "pident": [75.0],
            "length": [50],
        })
        result = assign_quality_tiers(df)
        assert result.iloc[0]["tier"] == "moderate"

    def test_boundary_just_below_moderate(self):
        """Hit just below moderate boundary falls to relaxed."""
        df = pd.DataFrame({
            "pident": [74.9],
            "length": [50],
        })
        result = assign_quality_tiers(df)
        assert result.iloc[0]["tier"] == "relaxed"

    def test_empty_dataframe(self, empty_df):
        """Empty input returns empty DataFrame with tier column."""
        result = assign_quality_tiers(empty_df)
        assert "tier" in result.columns
        assert len(result) == 0

    def test_does_not_modify_original(self, sample_hits):
        """Function returns a copy; original DataFrame is not modified."""
        original_cols = set(sample_hits.columns)
        _ = assign_quality_tiers(sample_hits)
        assert "tier" not in sample_hits.columns
        assert set(sample_hits.columns) == original_cols

    def test_strict_requires_both_conditions(self):
        """High pident but short length should NOT be strict."""
        df = pd.DataFrame({
            "pident": [95.0],
            "length": [50],  # below strict_length=100
        })
        result = assign_quality_tiers(df)
        assert result.iloc[0]["tier"] == "moderate"  # meets moderate thresholds


# ── compute_edit_stats ────────────────────────────────────────────────────────


class TestComputeEditStats:
    def test_mismatch_rate(self, sample_hits):
        """mismatch_rate = mismatch / length."""
        result = compute_edit_stats(sample_hits)
        # Row 0: 10 / 200 = 0.05
        assert result.iloc[0]["mismatch_rate"] == pytest.approx(0.05)
        # Row 1: 16 / 80 = 0.20
        assert result.iloc[1]["mismatch_rate"] == pytest.approx(0.20)

    def test_gap_rate(self, sample_hits):
        """gap_rate = gapopen / length."""
        result = compute_edit_stats(sample_hits)
        # Row 0: 1 / 200 = 0.005
        assert result.iloc[0]["gap_rate"] == pytest.approx(0.005)
        # Row 4: 3 / 60 = 0.05
        assert result.iloc[4]["gap_rate"] == pytest.approx(0.05)

    def test_edit_distance(self, sample_hits):
        """edit_distance = mismatch + gapopen."""
        result = compute_edit_stats(sample_hits)
        # Row 0: 10 + 1 = 11
        assert result.iloc[0]["edit_distance"] == 11
        # Row 4: 18 + 3 = 21
        assert result.iloc[4]["edit_distance"] == 21

    def test_zero_length_protection(self):
        """Zero-length hits use clipped length of 1 to avoid division by zero."""
        df = pd.DataFrame({
            "pident": [80.0],
            "length": [0],
            "mismatch": [5],
            "gapopen": [2],
        })
        result = compute_edit_stats(df)
        # Clipped to length=1: mismatch_rate = 5/1 = 5.0
        assert result.iloc[0]["mismatch_rate"] == pytest.approx(5.0)
        assert result.iloc[0]["gap_rate"] == pytest.approx(2.0)
        assert result.iloc[0]["edit_distance"] == 7

    def test_empty_dataframe(self, empty_df):
        """Empty input returns empty DataFrame with edit columns."""
        result = compute_edit_stats(empty_df)
        assert "mismatch_rate" in result.columns
        assert "gap_rate" in result.columns
        assert "edit_distance" in result.columns
        assert len(result) == 0

    def test_does_not_modify_original(self, sample_hits):
        """Function returns a copy; original DataFrame is not modified."""
        original_cols = set(sample_hits.columns)
        _ = compute_edit_stats(sample_hits)
        assert "mismatch_rate" not in sample_hits.columns
        assert set(sample_hits.columns) == original_cols


# ── summarize_tiers ───────────────────────────────────────────────────────────


class TestSummarizeTiers:
    def test_correct_grouping(self, sample_hits):
        """Summary groups by tier with correct counts."""
        tiered = assign_quality_tiers(sample_hits)
        summary = summarize_tiers(tiered)
        # Expect: strict=2 (rows 0,3), moderate=1 (row 1), relaxed=2 (rows 2,4)
        assert summary.loc["strict", "n_hits"] == 2
        assert summary.loc["moderate", "n_hits"] == 1
        assert summary.loc["relaxed", "n_hits"] == 2

    def test_percentages(self, sample_hits):
        """Percentages should sum to 100."""
        tiered = assign_quality_tiers(sample_hits)
        summary = summarize_tiers(tiered)
        assert summary["pct"].sum() == pytest.approx(100.0)

    def test_mean_values(self, sample_hits):
        """Mean pident/length/evalue are correct per tier."""
        tiered = assign_quality_tiers(sample_hits)
        summary = summarize_tiers(tiered)
        # Strict tier: rows 0 (pident=95) and 3 (pident=90) -> mean = 92.5
        assert summary.loc["strict", "mean_pident"] == pytest.approx(92.5)
        # Strict tier: rows 0 (length=200) and 3 (length=150) -> mean = 175.0
        assert summary.loc["strict", "mean_length"] == pytest.approx(175.0)
        # Moderate tier: row 1 only (pident=80)
        assert summary.loc["moderate", "mean_pident"] == pytest.approx(80.0)

    def test_empty_dataframe(self, empty_df):
        """Empty input returns empty summary with correct columns."""
        result = summarize_tiers(empty_df)
        assert list(result.columns) == [
            "n_hits", "pct", "mean_pident", "mean_length", "mean_evalue",
        ]
        assert len(result) == 0

    def test_single_tier(self):
        """All hits in one tier still produces valid summary."""
        df = pd.DataFrame({
            "pident": [95.0, 90.0],
            "length": [200, 150],
            "evalue": [1e-10, 1e-8],
            "tier": ["strict", "strict"],
        })
        summary = summarize_tiers(df)
        assert len(summary) == 1
        assert summary.loc["strict", "n_hits"] == 2
        assert summary.loc["strict", "pct"] == pytest.approx(100.0)


# ── compute_tier_edit_summary ─────────────────────────────────────────────────


class TestComputeTierEditSummary:
    def test_per_tier_edit_summary(self, sample_hits):
        """Edit summary computes correct mean rates per tier."""
        tiered = assign_quality_tiers(sample_hits)
        with_edits = compute_edit_stats(tiered)
        summary = compute_tier_edit_summary(with_edits)

        # Strict tier: rows 0 and 3
        # Row 0: mismatch_rate = 10/200 = 0.05, gap_rate = 1/200 = 0.005
        # Row 3: mismatch_rate = 15/150 = 0.10, gap_rate = 2/150 = 0.01333
        expected_mismatch = (0.05 + 0.10) / 2
        expected_gap = (0.005 + 2 / 150) / 2
        assert summary.loc["strict", "mean_mismatch_rate"] == pytest.approx(expected_mismatch)
        assert summary.loc["strict", "mean_gap_rate"] == pytest.approx(expected_gap)

    def test_all_tiers_present(self, sample_hits):
        """All tiers that appear in data are present in the summary."""
        tiered = assign_quality_tiers(sample_hits)
        with_edits = compute_edit_stats(tiered)
        summary = compute_tier_edit_summary(with_edits)
        assert "strict" in summary.index
        assert "moderate" in summary.index
        assert "relaxed" in summary.index

    def test_empty_dataframe(self, empty_df):
        """Empty input returns empty summary with correct columns."""
        result = compute_tier_edit_summary(empty_df)
        assert list(result.columns) == ["mean_mismatch_rate", "mean_gap_rate"]
        assert len(result) == 0
