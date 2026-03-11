"""Tests for strand bias analysis."""

import pandas as pd
import pytest

from fossil_finder.analysis.strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_te_strand_bias,
    compute_genome_strand_bias,
)


@pytest.fixture
def blast_hits_with_genes():
    """BLAST hits with gene_id column already populated."""
    return pd.DataFrame({
        "qseqid": ["tr1", "tr1", "tr1", "tr2", "tr2"],
        "sseqid": ["TE_gypsy", "TE_gypsy", "TE_copia", "TE_gypsy", "TE_jockey"],
        "gene_id": ["gA", "gA", "gA", "gB", "gB"],
        "strand": ["plus", "plus", "minus", "plus", "minus"],
        "length": [100, 80, 60, 120, 90],
    })


class TestClassifyBias:
    def test_strong_sense(self):
        assert classify_bias(0.85) == "strong_sense"

    def test_sense(self):
        assert classify_bias(0.60) == "sense"

    def test_balanced(self):
        assert classify_bias(0.50) == "balanced"

    def test_antisense(self):
        assert classify_bias(0.35) == "antisense"

    def test_strong_antisense(self):
        assert classify_bias(0.15) == "strong_antisense"


class TestComputeGeneStrandBias:
    def test_per_gene_counts(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "sense_hits"] == 2
        assert result.loc["gA", "antisense_hits"] == 1
        assert result.loc["gB", "sense_hits"] == 1
        assert result.loc["gB", "antisense_hits"] == 1

    def test_sense_pct(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "sense_pct"] == pytest.approx(2 / 3 * 100)
        assert result.loc["gB", "sense_pct"] == pytest.approx(50.0)

    def test_bias_classification(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "bias"] == "sense"
        assert result.loc["gB", "bias"] == "balanced"

    def test_empty_input(self):
        df = pd.DataFrame(columns=["gene_id", "strand", "length"])
        result = compute_gene_strand_bias(df)
        assert len(result) == 0


class TestComputeTEStrandBias:
    def test_per_te_counts(self, blast_hits_with_genes):
        result = compute_te_strand_bias(blast_hits_with_genes)
        assert result.loc["TE_gypsy", "sense_hits"] == 3
        assert result.loc["TE_gypsy", "antisense_hits"] == 0

    def test_min_hits_filter(self, blast_hits_with_genes):
        result = compute_te_strand_bias(blast_hits_with_genes, min_hits=3)
        # Only TE_gypsy has 3+ hits
        assert "TE_gypsy" in result.index
        assert "TE_copia" not in result.index
        assert "TE_jockey" not in result.index


class TestComputeGenomeStrandBias:
    def test_genome_totals(self, blast_hits_with_genes):
        result = compute_genome_strand_bias(blast_hits_with_genes)
        assert result["total_hits"] == 5
        assert result["sense_hits"] == 3
        assert result["antisense_hits"] == 2
        assert result["sense_pct"] == pytest.approx(60.0)
