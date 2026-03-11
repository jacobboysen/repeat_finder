"""Tests for individual pipeline steps."""

import json
from pathlib import Path

import pandas as pd
import pytest

from fossil_finder.pipeline.steps import (
    step_extract_regions,
    step_load_and_filter,
    step_deduplicate,
    step_aggregate,
    step_strand_analysis,
    step_family_analysis,
    step_enrichment_analysis,
    step_repeatmasker_overlap,
)


@pytest.fixture
def blast_df():
    """Minimal BLAST results DataFrame."""
    return pd.DataFrame({
        "qseqid": ["utr1", "utr1", "utr2", "utr2", "utr3"],
        "sseqid": ["TE_gypsy1", "TE_copia1", "TE_gypsy1", "TE_jockey1", "TE_gypsy1"],
        "pident": [80.0, 65.0, 85.0, 70.0, 90.0],
        "length": [120, 85, 200, 60, 150],
        "mismatch": [20, 25, 30, 15, 12],
        "gapopen": [3, 4, 5, 2, 1],
        "qstart": [50, 200, 10, 300, 100],
        "qend": [170, 284, 210, 360, 250],
        "sstart": [1, 300, 50, 100, 1],
        "send": [120, 216, 250, 160, 150],
        "evalue": [1.5e-10, 3.2e-5, 2.1e-15, 0.005, 8.3e-20],
        "bitscore": [85.2, 42.1, 120.5, 30.8, 155.0],
        "qlen": [500, 500, 800, 800, 600],
        "slen": [5000, 4500, 5000, 3000, 4000],
        "qseq": ["A"] * 5,
        "sseq": ["A"] * 5,
        "strand": ["plus", "minus", "plus", "plus", "plus"],
    })


@pytest.fixture
def query_to_gene():
    return {"utr1": "gene1", "utr2": "gene2", "utr3": "gene3"}


class TestStepExtractRegions:
    def test_extracts_to_fasta_and_metadata(self, custom_config, tmp_path):
        fasta_path = tmp_path / "regions.fa"
        meta_path = tmp_path / "regions.tsv"
        regions = step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
        )
        assert fasta_path.exists()
        assert meta_path.exists()
        assert len(regions) > 0
        assert "sequence" in regions[0]

    def test_returns_region_dicts(self, custom_config, tmp_path):
        regions = step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=tmp_path / "r.fa",
            metadata_out=tmp_path / "r.tsv",
        )
        for r in regions:
            assert "region_id" in r
            assert "chrom" in r
            assert "start" in r

    def test_skip_if_exists(self, custom_config, tmp_path):
        fasta_path = tmp_path / "regions.fa"
        meta_path = tmp_path / "regions.tsv"
        # First run
        step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
        )
        mtime = fasta_path.stat().st_mtime
        # Second run with force=False should skip
        step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
            force=False,
        )
        assert fasta_path.stat().st_mtime == mtime


class TestStepLoadAndFilter:
    def test_loads_from_path(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results)
        assert len(df) == 5
        assert "strand" in df.columns

    def test_applies_evalue_filter(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results, max_evalue=1e-10)
        assert all(df["evalue"] <= 1e-10)

    def test_applies_pident_filter(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results, min_pident=80.0)
        assert all(df["pident"] >= 80.0)


class TestStepDeduplicate:
    def test_removes_duplicates(self):
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 1],
            "te_end": [100, 100],
            "source_transcript": ["tr1", "tr2"],
        })
        result, stats = step_deduplicate(df)
        assert len(result) == 1
        assert stats["duplicates_removed"] == 1

    def test_no_duplicates_passthrough(self, blast_df):
        # blast_df lacks default dedup key columns, so specify BLAST-level keys
        result, stats = step_deduplicate(
            blast_df,
            key_columns=["qseqid", "sseqid", "qstart", "qend"],
        )
        assert stats["duplicates_removed"] == 0


class TestStepAggregate:
    def test_produces_gene_stats(self, blast_df, query_to_gene):
        result = step_aggregate(blast_df, query_to_gene)
        assert "gene1" in result.index
        assert "density" in result.columns
        assert result.loc["gene1", "hit_count"] == 2

    def test_empty_input(self):
        df = pd.DataFrame(columns=["qseqid", "sseqid", "length", "strand", "qlen"])
        result = step_aggregate(df, {})
        assert len(result) == 0


class TestStepStrandAnalysis:
    def test_returns_three_levels(self, blast_df):
        blast_df["gene_id"] = blast_df["qseqid"].map(
            {"utr1": "g1", "utr2": "g2", "utr3": "g3"}
        )
        result = step_strand_analysis(blast_df)
        assert "gene" in result
        assert "te_family" in result
        assert "genome" in result


class TestStepFamilyAnalysis:
    def test_returns_family_stats(self, blast_df):
        result = step_family_analysis(blast_df)
        assert "family_stats" in result
        assert "TE_gypsy1" in result["family_stats"].index

    def test_class_distribution_with_metadata(self, blast_df):
        te_meta = {"TE_gypsy1": {"te_class": "LTR"}, "TE_copia1": {"te_class": "LTR"}}
        result = step_family_analysis(blast_df, te_metadata=te_meta)
        assert "class_distribution" in result
        assert "LTR" in result["class_distribution"]


class TestStepEnrichmentAnalysis:
    def test_enrichment_for_gene_set(self):
        """Use a 4-gene setup so Mann-Whitney gets >= 2 values per group."""
        gene_stats = pd.DataFrame({
            "hit_count": [5, 3, 1, 2],
            "hit_bp": [500, 300, 100, 200],
            "query_len": [1000, 800, 600, 700],
            "density": [500.0, 375.0, 166.7, 285.7],
        }, index=["gene1", "gene2", "gene3", "gene4"])
        te_positive = {"gene1", "gene2", "gene3", "gene4"}
        result = step_enrichment_analysis(
            gene_set={"gene1", "gene2"},
            te_positive_genes=te_positive,
            gene_densities=gene_stats["density"],
        )
        assert "fisher" in result
        assert "mannwhitney" in result
        # Both groups have 2 values, so Mann-Whitney should produce a real p-value
        assert result["mannwhitney"]["p_value"] is not None


class TestStepRepeatMaskerOverlap:
    def test_classifies_hits(self, blast_df, mini_repeatmasker):
        # Need qseqid-based region mapping for RM overlap
        query_regions = pd.DataFrame({
            "region_id": ["utr1", "utr2", "utr3"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [1, 200, 50],
            "end": [500, 600, 400],
        })
        result = step_repeatmasker_overlap(
            blast_hits=blast_df,
            repeatmasker_path=mini_repeatmasker,
            query_regions=query_regions,
        )
        assert "known" in result
        assert "novel" in result
        assert "rm_stats" in result
