"""Tests for TE family distribution analysis."""

import pandas as pd
import pytest

from fossil_finder.analysis.families import (
    compute_family_stats,
    compute_class_distribution,
    compute_fold_enrichment,
)


@pytest.fixture
def blast_hits():
    return pd.DataFrame({
        "sseqid": ["TE_gypsy", "TE_gypsy", "TE_gypsy", "TE_copia", "TE_jockey"],
        "strand": ["plus", "plus", "minus", "plus", "minus"],
        "pident": [80.0, 75.0, 85.0, 70.0, 90.0],
        "evalue": [1e-10, 1e-8, 1e-12, 1e-5, 1e-15],
        "length": [100, 80, 120, 60, 150],
    })


@pytest.fixture
def te_metadata():
    return {
        "TE_gypsy": {"te_class": "LTR"},
        "TE_copia": {"te_class": "LTR"},
        "TE_jockey": {"te_class": "LINE"},
    }


class TestComputeFamilyStats:
    def test_hit_counts(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "hit_count"] == 3
        assert result.loc["TE_copia", "hit_count"] == 1
        assert result.loc["TE_jockey", "hit_count"] == 1

    def test_strand_counts(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "sense_hits"] == 2
        assert result.loc["TE_gypsy", "antisense_hits"] == 1

    def test_mean_metrics(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "mean_pident"] == pytest.approx(80.0)
        assert result.loc["TE_gypsy", "total_bp"] == 300  # 100+80+120

    def test_frequency(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "frequency"] == pytest.approx(3 / 5)

    def test_empty_input(self):
        df = pd.DataFrame(columns=["sseqid", "strand", "pident", "evalue", "length"])
        result = compute_family_stats(df)
        assert len(result) == 0


class TestComputeClassDistribution:
    def test_class_counts(self, blast_hits, te_metadata):
        result = compute_class_distribution(blast_hits, te_metadata)
        assert result["LTR"] == 4  # gypsy(3) + copia(1)
        assert result["LINE"] == 1

    def test_unknown_class(self, blast_hits):
        """TEs not in metadata get class 'Unknown'."""
        result = compute_class_distribution(blast_hits, {})
        assert result["Unknown"] == 5


class TestComputeFoldEnrichment:
    def test_enrichment_calculation(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1", "TE1", "TE2"],
            "strand": ["plus"] * 4,
            "pident": [80.0] * 4,
            "evalue": [1e-5] * 4,
            "length": [100] * 4,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE1", "TE2", "TE2", "TE2"],
            "strand": ["plus"] * 4,
            "pident": [80.0] * 4,
            "evalue": [1e-5] * 4,
            "length": [100] * 4,
        })
        result = compute_fold_enrichment(set_a, set_b)
        # TE1: freq_a=3/4=0.75, freq_b=1/4=0.25 -> enrichment=3.0
        assert result.loc["TE1", "fold_enrichment"] == pytest.approx(3.0)
        # TE2: freq_a=1/4=0.25, freq_b=3/4=0.75 -> enrichment=1/3
        assert result.loc["TE2", "fold_enrichment"] == pytest.approx(1 / 3)

    def test_log2_enrichment(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1"], "strand": ["plus"] * 2,
            "pident": [80.0] * 2, "evalue": [1e-5] * 2, "length": [100] * 2,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE1"], "strand": ["plus"],
            "pident": [80.0], "evalue": [1e-5], "length": [100],
        })
        result = compute_fold_enrichment(set_a, set_b)
        # freq_a=1.0, freq_b=1.0 -> enrichment=1.0, log2=0.0
        # Both sets have TE1 at 100% frequency -> fold=1.0
        assert result.loc["TE1", "log2_enrichment"] == pytest.approx(0.0)

    def test_family_unique_to_one_set(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1"], "strand": ["plus"] * 2,
            "pident": [80.0] * 2, "evalue": [1e-5] * 2, "length": [100] * 2,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE2"], "strand": ["plus"],
            "pident": [80.0], "evalue": [1e-5], "length": [100],
        })
        result = compute_fold_enrichment(set_a, set_b)
        # TE1 only in set_a -> fold_enrichment = inf, log2 = inf
        assert result.loc["TE1", "fold_enrichment"] == float("inf")
        assert result.loc["TE1", "log2_enrichment"] == float("inf")
        # TE2 only in set_b -> fold_enrichment = 0.0, log2 = -inf
        assert result.loc["TE2", "fold_enrichment"] == 0.0
        assert result.loc["TE2", "log2_enrichment"] == float("-inf")
