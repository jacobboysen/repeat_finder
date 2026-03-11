"""Tests for per-gene hit aggregation."""

import pandas as pd
import pytest

from fossil_finder.analysis.aggregation import aggregate_by_gene, compute_density


@pytest.fixture
def blast_hits():
    """BLAST results with multiple transcripts mapping to same genes."""
    return pd.DataFrame({
        "qseqid": ["tr1", "tr1", "tr2", "tr2", "tr3"],
        "sseqid": ["TE1", "TE2", "TE1", "TE3", "TE1"],
        "pident": [80.0, 75.0, 90.0, 60.0, 85.0],
        "length": [100, 50, 200, 80, 150],
        "evalue": [1e-10, 1e-5, 1e-20, 0.01, 1e-12],
        "bitscore": [80.0, 40.0, 150.0, 30.0, 100.0],
        "qstart": [10, 200, 50, 300, 100],
        "qend": [110, 250, 250, 380, 250],
        "sstart": [1, 1, 1, 1, 1],
        "send": [100, 50, 200, 80, 150],
        "qlen": [500, 500, 800, 800, 600],
        "slen": [5000, 4000, 5000, 3000, 5000],
        "strand": ["plus", "minus", "plus", "plus", "plus"],
    })


@pytest.fixture
def query_to_gene():
    """Mapping from transcript/query ID to gene ID."""
    return {"tr1": "geneA", "tr2": "geneA", "tr3": "geneB"}


class TestAggregateByGene:
    def test_groups_by_gene(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert len(result) == 2  # geneA, geneB
        assert "geneA" in result.index
        assert "geneB" in result.index

    def test_hit_count(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "hit_count"] == 4  # tr1(2) + tr2(2)
        assert result.loc["geneB", "hit_count"] == 1

    def test_hit_bp(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        # geneA: 100 + 50 + 200 + 80 = 430
        assert result.loc["geneA", "hit_bp"] == 430
        # geneB: 150
        assert result.loc["geneB", "hit_bp"] == 150

    def test_sense_antisense_counts(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "sense_hits"] == 3
        assert result.loc["geneA", "antisense_hits"] == 1
        assert result.loc["geneB", "sense_hits"] == 1
        assert result.loc["geneB", "antisense_hits"] == 0

    def test_query_length_uses_max_per_gene(self, blast_hits, query_to_gene):
        """Gene query length = max qlen across its transcripts."""
        result = aggregate_by_gene(blast_hits, query_to_gene)
        # geneA has tr1(500) and tr2(800) -> max = 800
        assert result.loc["geneA", "query_len"] == 800
        assert result.loc["geneB", "query_len"] == 600

    def test_unmapped_queries_skipped(self, blast_hits):
        partial_map = {"tr1": "geneA"}  # tr2, tr3 not mapped
        result = aggregate_by_gene(blast_hits, partial_map)
        assert len(result) == 1
        assert result.loc["geneA", "hit_count"] == 2

    def test_empty_dataframe(self, query_to_gene):
        df = pd.DataFrame(columns=["qseqid", "sseqid", "length", "strand", "qlen"])
        result = aggregate_by_gene(df, query_to_gene)
        assert len(result) == 0

    def test_te_families_tracked(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "n_te_families"] == 3  # TE1, TE2, TE3
        assert result.loc["geneB", "n_te_families"] == 1


class TestComputeDensity:
    def test_density_calculation(self):
        """Density = hit_bp / query_len * 1000 (per kb)."""
        gene_stats = pd.DataFrame({
            "hit_bp": [500, 100],
            "query_len": [1000, 500],
        }, index=["geneA", "geneB"])
        result = compute_density(gene_stats)
        assert result.loc["geneA", "density"] == pytest.approx(500.0)
        assert result.loc["geneB", "density"] == pytest.approx(200.0)

    def test_density_zero_length_safe(self):
        gene_stats = pd.DataFrame({
            "hit_bp": [100],
            "query_len": [0],
        }, index=["geneA"])
        result = compute_density(gene_stats)
        assert result.loc["geneA", "density"] == 0.0

    def test_density_added_as_column(self, blast_hits, query_to_gene):
        agg = aggregate_by_gene(blast_hits, query_to_gene)
        result = compute_density(agg)
        assert "density" in result.columns
