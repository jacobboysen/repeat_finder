"""Tests for statistical enrichment testing."""

import numpy as np
import pandas as pd
import pytest

from fossil_finder.analysis.enrichment import (
    fisher_exact_enrichment,
    mannwhitney_enrichment,
    correct_multiple_testing,
    test_gene_set_enrichment,
)


class TestFisherExactEnrichment:
    def test_enriched_set(self):
        """Set with more TE-positive genes than background -> OR > 1."""
        result = fisher_exact_enrichment(
            set_positive=8, set_negative=2,
            bg_positive=20, bg_negative=80,
        )
        assert result["odds_ratio"] > 1.0
        assert result["p_value"] < 0.05

    def test_depleted_set(self):
        """Set with fewer TE-positive genes -> OR < 1."""
        result = fisher_exact_enrichment(
            set_positive=1, set_negative=9,
            bg_positive=50, bg_negative=50,
        )
        assert result["odds_ratio"] < 1.0

    def test_no_difference(self):
        """Equal proportions -> OR ~ 1."""
        result = fisher_exact_enrichment(
            set_positive=5, set_negative=5,
            bg_positive=50, bg_negative=50,
        )
        assert result["odds_ratio"] == pytest.approx(1.0, abs=0.5)
        assert result["p_value"] > 0.05

    def test_returns_required_keys(self):
        result = fisher_exact_enrichment(5, 5, 50, 50)
        assert "odds_ratio" in result
        assert "p_value" in result


class TestMannWhitneyEnrichment:
    def test_higher_group(self):
        """Group with clearly higher values -> p < 0.05."""
        group_a = [10.0, 12.0, 15.0, 11.0, 14.0]
        group_b = [1.0, 2.0, 3.0, 1.5, 2.5]
        result = mannwhitney_enrichment(group_a, group_b)
        assert result["p_value"] < 0.05

    def test_equal_groups(self):
        """Identical distributions -> p > 0.05."""
        group = [5.0, 5.0, 5.0, 5.0, 5.0]
        result = mannwhitney_enrichment(group, group)
        assert result["p_value"] > 0.05

    def test_returns_required_keys(self):
        result = mannwhitney_enrichment([1.0, 2.0], [3.0, 4.0])
        assert "u_statistic" in result
        assert "p_value" in result

    def test_too_few_samples(self):
        """With < 2 samples, return NaN p-value."""
        result = mannwhitney_enrichment([1.0], [2.0, 3.0])
        assert np.isnan(result["p_value"])


class TestCorrectMultipleTesting:
    def test_bonferroni_correction(self):
        p_values = [0.01, 0.04, 0.06]
        q_values = correct_multiple_testing(p_values, method="bonferroni")
        assert len(q_values) == 3
        # Bonferroni: p * n_tests, capped at 1.0
        assert q_values[0] == pytest.approx(0.03)
        assert q_values[1] == pytest.approx(0.12)
        assert q_values[2] == pytest.approx(0.18)

    def test_fdr_correction(self):
        p_values = [0.01, 0.04, 0.06]
        q_values = correct_multiple_testing(p_values, method="fdr")
        assert len(q_values) == 3
        # BH: rank 1 -> 0.01*3/1=0.03, rank 2 -> 0.04*3/2=0.06, rank 3 -> 0.06*3/3=0.06
        # After monotonicity sweep: [0.03, 0.06, 0.06]
        assert q_values[0] == pytest.approx(0.03)
        assert q_values[1] == pytest.approx(0.06)
        assert q_values[2] == pytest.approx(0.06)

    def test_empty_input(self):
        assert correct_multiple_testing([], method="fdr") == []

    def test_single_pvalue(self):
        q = correct_multiple_testing([0.05], method="bonferroni")
        assert q[0] == pytest.approx(0.05)


class TestGeneSetEnrichment:
    def test_full_enrichment_pipeline(self):
        """Test the combined pipeline with gene set vs background."""
        gene_densities = pd.Series(
            [5.0, 8.0, 12.0, 3.0, 15.0],
            index=["g1", "g2", "g3", "g4", "g5"],
        )
        gene_set = {"g1", "g2", "g3"}  # genes of interest
        te_positive = {"g1", "g2", "g3", "g5"}  # genes with TE hits

        result = test_gene_set_enrichment(
            gene_set=gene_set,
            te_positive_genes=te_positive,
            gene_densities=gene_densities,
        )
        assert "fisher" in result
        assert "mannwhitney" in result
        assert result["fisher"]["odds_ratio"] is not None
        assert result["mannwhitney"]["p_value"] is not None
