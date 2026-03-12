"""Tests for hit multiplicity analysis."""

import pandas as pd
import pytest

from fossil_finder.analysis.multiplicity import (
    compute_hit_multiplicity,
    compute_query_hit_counts,
    compute_te_breadth,
)


@pytest.fixture
def blast_hits():
    """Multi-query, multi-TE BLAST results for multiplicity testing."""
    return pd.DataFrame({
        "qseqid": ["q1", "q1", "q1", "q2", "q2", "q3"],
        "sseqid": ["TE_A", "TE_B", "TE_A", "TE_A", "TE_C", "TE_B"],
        "pident": [80.0, 75.0, 85.0, 70.0, 90.0, 60.0],
    })


@pytest.fixture
def query_to_gene():
    """Mapping: q1 and q2 belong to gene g1, q3 to gene g2."""
    return {"q1": "g1", "q2": "g1", "q3": "g2"}


class TestComputeHitMultiplicity:
    def test_basic_stats(self, blast_hits):
        """Verify median, mean, max, n for hits_per_query and queries_per_te."""
        result = compute_hit_multiplicity(blast_hits)

        hpq = result["hits_per_query"]
        # q1: 3 hits, q2: 2 hits, q3: 1 hit
        assert hpq["n"] == 3
        assert hpq["max"] == 3
        assert hpq["median"] == pytest.approx(2.0)
        assert hpq["mean"] == pytest.approx(2.0)

        qpt = result["queries_per_te"]
        # TE_A: q1,q2 (2), TE_B: q1,q3 (2), TE_C: q2 (1)
        assert qpt["n"] == 3
        assert qpt["max"] == 2
        assert qpt["median"] == pytest.approx(2.0)
        assert qpt["mean"] == pytest.approx(5 / 3)

    def test_with_gene_mapping(self, blast_hits, query_to_gene):
        """Gene-level metrics are present when query_to_gene is provided."""
        result = compute_hit_multiplicity(blast_hits, query_to_gene=query_to_gene)

        assert "genes_per_te" in result
        assert "hits_per_gene" in result

        gpt = result["genes_per_te"]
        # TE_A: g1 (via q1,q2) -> 1 gene; TE_B: g1 (q1), g2 (q3) -> 2 genes;
        # TE_C: g1 (q2) -> 1 gene
        assert gpt["n"] == 3
        assert gpt["max"] == 2

        hpg = result["hits_per_gene"]
        # g1: 5 hits (q1:3 + q2:2), g2: 1 hit (q3:1)
        assert hpg["n"] == 2
        assert hpg["max"] == 5
        assert hpg["mean"] == pytest.approx(3.0)

    def test_without_gene_mapping(self, blast_hits):
        """No gene-level keys when query_to_gene is None."""
        result = compute_hit_multiplicity(blast_hits, query_to_gene=None)
        assert "genes_per_te" not in result
        assert "hits_per_gene" not in result

    def test_empty_dataframe(self):
        """Empty input returns zeroed-out summaries."""
        df = pd.DataFrame(columns=["qseqid", "sseqid", "pident"])
        result = compute_hit_multiplicity(df)

        assert result["hits_per_query"]["n"] == 0
        assert result["hits_per_query"]["max"] == 0
        assert result["queries_per_te"]["mean"] == pytest.approx(0.0)

    def test_empty_dataframe_with_gene_mapping(self):
        """Empty input with gene mapping still returns gene-level keys."""
        df = pd.DataFrame(columns=["qseqid", "sseqid", "pident"])
        result = compute_hit_multiplicity(df, query_to_gene={"q1": "g1"})

        assert "genes_per_te" in result
        assert "hits_per_gene" in result
        assert result["genes_per_te"]["n"] == 0
        assert result["hits_per_gene"]["n"] == 0

    def test_single_hit(self):
        """A single hit produces consistent multiplicity stats."""
        df = pd.DataFrame({
            "qseqid": ["q1"],
            "sseqid": ["TE_A"],
            "pident": [85.0],
        })
        result = compute_hit_multiplicity(df)

        for key in ("hits_per_query", "queries_per_te"):
            assert result[key]["n"] == 1
            assert result[key]["max"] == 1
            assert result[key]["median"] == pytest.approx(1.0)
            assert result[key]["mean"] == pytest.approx(1.0)

    def test_unmapped_queries_excluded_from_gene_metrics(self):
        """Queries not in query_to_gene are excluded from gene-level stats."""
        df = pd.DataFrame({
            "qseqid": ["q1", "q_unknown"],
            "sseqid": ["TE_A", "TE_B"],
            "pident": [80.0, 70.0],
        })
        result = compute_hit_multiplicity(df, query_to_gene={"q1": "g1"})

        # Only q1 maps, so hits_per_gene has just g1 with 1 hit
        assert result["hits_per_gene"]["n"] == 1
        assert result["hits_per_gene"]["max"] == 1


class TestComputeTeBreadth:
    def test_multiple_tes(self, blast_hits):
        """Verify per-TE breadth counts and identity stats."""
        result = compute_te_breadth(blast_hits)

        assert result.loc["TE_A", "n_hits"] == 3
        assert result.loc["TE_A", "n_queries"] == 2  # q1, q2
        assert result.loc["TE_B", "n_hits"] == 2
        assert result.loc["TE_B", "n_queries"] == 2  # q1, q3
        assert result.loc["TE_C", "n_hits"] == 1
        assert result.loc["TE_C", "n_queries"] == 1  # q2

        # pident for TE_A: 80, 85, 70 -> mean=78.33, median=80
        assert result.loc["TE_A", "mean_pident"] == pytest.approx(78.333, abs=0.01)
        assert result.loc["TE_A", "median_pident"] == pytest.approx(80.0)

    def test_with_gene_mapping(self, blast_hits, query_to_gene):
        """n_genes column present when query_to_gene is provided."""
        result = compute_te_breadth(blast_hits, query_to_gene=query_to_gene)

        assert "n_genes" in result.columns
        # TE_A hit by q1 (g1) and q2 (g1) -> 1 unique gene
        assert result.loc["TE_A", "n_genes"] == 1
        # TE_B hit by q1 (g1) and q3 (g2) -> 2 unique genes
        assert result.loc["TE_B", "n_genes"] == 2
        # TE_C hit by q2 (g1) -> 1 unique gene
        assert result.loc["TE_C", "n_genes"] == 1

    def test_without_gene_mapping(self, blast_hits):
        """n_genes column absent when query_to_gene is None."""
        result = compute_te_breadth(blast_hits, query_to_gene=None)
        assert "n_genes" not in result.columns

    def test_empty_dataframe(self):
        """Empty input returns empty DataFrame with correct columns."""
        df = pd.DataFrame(columns=["qseqid", "sseqid", "pident"])
        result = compute_te_breadth(df)
        assert len(result) == 0
        assert "n_hits" in result.columns
        assert "mean_pident" in result.columns

    def test_empty_dataframe_with_gene_mapping(self):
        """Empty input with gene mapping has n_genes column."""
        df = pd.DataFrame(columns=["qseqid", "sseqid", "pident"])
        result = compute_te_breadth(df, query_to_gene={"q1": "g1"})
        assert len(result) == 0
        assert "n_genes" in result.columns


class TestComputeQueryHitCounts:
    def test_basic_counting(self, blast_hits):
        """Each query gets its correct hit count."""
        result = compute_query_hit_counts(blast_hits)

        assert result.loc["q1", "hit_count"] == 3
        assert result.loc["q2", "hit_count"] == 2
        assert result.loc["q3", "hit_count"] == 1

    def test_single_query(self):
        """Single query returns one-row DataFrame."""
        df = pd.DataFrame({
            "qseqid": ["q1", "q1"],
            "sseqid": ["TE_A", "TE_B"],
        })
        result = compute_query_hit_counts(df)

        assert len(result) == 1
        assert result.loc["q1", "hit_count"] == 2

    def test_empty_dataframe(self):
        """Empty input returns empty DataFrame with hit_count column."""
        df = pd.DataFrame(columns=["qseqid", "sseqid"])
        result = compute_query_hit_counts(df)

        assert len(result) == 0
        assert "hit_count" in result.columns

    def test_index_name(self, blast_hits):
        """Index is named qseqid."""
        result = compute_query_hit_counts(blast_hits)
        assert result.index.name == "qseqid"
