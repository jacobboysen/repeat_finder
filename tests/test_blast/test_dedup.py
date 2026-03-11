"""Tests for BLAST hit deduplication."""

import pandas as pd
import pytest

from fossil_finder.blast.dedup import HitDeduplicator


DEFAULT_KEY = ["gene_id", "chrom", "genomic_start", "genomic_end",
               "te_id", "te_start", "te_end"]


@pytest.fixture
def annotated_hits():
    """DataFrame simulating annotated BLAST hits with genomic coordinates."""
    return pd.DataFrame({
        "gene_id": ["g1", "g1", "g1", "g2", "g2"],
        "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "genomic_start": [100, 100, 200, 500, 500],
        "genomic_end": [200, 200, 300, 600, 600],
        "te_id": ["TE_gypsy1", "TE_gypsy1", "TE_gypsy1", "TE_jockey1", "TE_jockey1"],
        "te_start": [1, 1, 50, 1, 1],
        "te_end": [100, 100, 150, 100, 100],
        "pident": [80.0, 80.0, 75.0, 90.0, 90.0],
        "evalue": [1e-10, 1e-10, 1e-5, 1e-15, 1e-15],
        "source_transcript": ["tr1", "tr2", "tr1", "tr3", "tr4"],
    })


class TestHitDeduplicator:
    def test_removes_exact_duplicates(self, annotated_hits):
        dedup = HitDeduplicator()
        result = dedup.deduplicate(annotated_hits)
        # Rows 0 and 1 are duplicates (same key, from different transcripts)
        # Rows 3 and 4 are duplicates
        assert len(result) == 3  # rows 0, 2, 3 (or 1, 2, 4)

    def test_keeps_different_te_regions(self, annotated_hits):
        dedup = HitDeduplicator()
        result = dedup.deduplicate(annotated_hits)
        # g1 has two distinct TE regions: (1,100) and (50,150)
        g1_hits = result[result["gene_id"] == "g1"]
        assert len(g1_hits) == 2

    def test_stats_tracking(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        assert dedup.stats["total_input"] == 5
        assert dedup.stats["unique_hits"] == 3
        assert dedup.stats["duplicates_removed"] == 2

    def test_custom_key_columns(self, annotated_hits):
        """Using a smaller key should produce more dedup."""
        dedup = HitDeduplicator(key_columns=["gene_id", "chrom", "te_id"])
        result = dedup.deduplicate(annotated_hits)
        # With this key: g1+chr1+TE_gypsy1 = 1 unique, g2+chr2+TE_jockey1 = 1 unique
        assert len(result) == 2

    def test_empty_dataframe(self):
        dedup = HitDeduplicator()
        df = pd.DataFrame(columns=DEFAULT_KEY + ["pident"])
        result = dedup.deduplicate(df)
        assert len(result) == 0
        assert dedup.stats["total_input"] == 0

    def test_no_duplicates(self):
        df = pd.DataFrame({
            "gene_id": ["g1", "g2"],
            "chrom": ["chr1", "chr2"],
            "genomic_start": [100, 200],
            "genomic_end": [200, 300],
            "te_id": ["TE1", "TE2"],
            "te_start": [1, 1],
            "te_end": [50, 50],
        })
        dedup = HitDeduplicator()
        result = dedup.deduplicate(df)
        assert len(result) == 2
        assert dedup.stats["duplicates_removed"] == 0


class TestTECoordinateNormalization:
    def test_normalizes_te_coords(self):
        """TE start/end should be normalized so start < end before keying."""
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 100],   # forward
            "te_end": [100, 1],     # reverse (same region, different strand)
            "pident": [80.0, 80.0],
        })
        dedup = HitDeduplicator(normalize_te_coords=True)
        result = dedup.deduplicate(df)
        # After normalization, both rows have te_start=1, te_end=100 → duplicate
        assert len(result) == 1

    def test_normalization_preserves_original_coords(self):
        """Original te_start/te_end values should not be mutated."""
        df = pd.DataFrame({
            "gene_id": ["g1", "g2"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 300],
            "genomic_end": [200, 400],
            "te_id": ["TE1", "TE2"],
            "te_start": [100, 200],   # reversed
            "te_end": [1, 50],        # reversed
        })
        dedup = HitDeduplicator(normalize_te_coords=True)
        result = dedup.deduplicate(df)
        # Original coords should be preserved in the output
        assert result.iloc[0]["te_start"] == 100
        assert result.iloc[0]["te_end"] == 1
        assert result.iloc[1]["te_start"] == 200
        assert result.iloc[1]["te_end"] == 50

    def test_without_normalization_keeps_both(self):
        """Without normalization, (1,100) and (100,1) are different keys."""
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 100],
            "te_end": [100, 1],
            "pident": [80.0, 80.0],
        })
        dedup = HitDeduplicator(normalize_te_coords=False)
        result = dedup.deduplicate(df)
        assert len(result) == 2


class TestDuplicationReport:
    def test_per_gene_stats(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        gene_stats = dedup.per_gene_stats(annotated_hits)
        assert "g1" in gene_stats
        assert gene_stats["g1"]["raw_hits"] == 3
        assert gene_stats["g1"]["unique_hits"] == 2
        assert gene_stats["g1"]["duplicates_removed"] == 1

    def test_per_gene_stats_empty(self):
        dedup = HitDeduplicator()
        df = pd.DataFrame(columns=DEFAULT_KEY)
        dedup.deduplicate(df)
        gene_stats = dedup.per_gene_stats(df)
        assert gene_stats == {}

    def test_duplication_rate(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        assert dedup.stats["duplication_rate"] == pytest.approx(2 / 5)
