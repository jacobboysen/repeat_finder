"""Tests for individual pipeline steps."""

import json
from pathlib import Path

import pandas as pd
import pytest

from fossil_finder.pipeline.steps import (
    step_conservation_analysis,
    step_extract_regions,
    step_load_and_filter,
    step_deduplicate,
    step_aggregate,
    step_multiplicity_analysis,
    step_positional_analysis,
    step_quality_tiers,
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


class TestStepQualityTiers:
    def test_returns_df_and_summary(self, blast_df):
        annotated, summary = step_quality_tiers(blast_df)
        assert isinstance(annotated, pd.DataFrame)
        assert isinstance(summary, pd.DataFrame)
        assert len(annotated) == len(blast_df)

    def test_tier_column_values(self, blast_df):
        annotated, _ = step_quality_tiers(blast_df)
        assert "tier" in annotated.columns
        valid_tiers = {"strict", "moderate", "relaxed"}
        assert set(annotated["tier"].unique()).issubset(valid_tiers)

    def test_tier_assignment_correctness(self):
        """Verify tier thresholds are applied correctly."""
        df = pd.DataFrame({
            "qseqid": ["q1", "q2", "q3"],
            "sseqid": ["TE1", "TE2", "TE3"],
            "pident": [90.0, 78.0, 60.0],
            "length": [150, 80, 30],
            "mismatch": [10, 15, 20],
            "gapopen": [1, 2, 3],
            "qstart": [1, 1, 1],
            "qend": [150, 80, 30],
            "sstart": [1, 1, 1],
            "send": [150, 80, 30],
            "evalue": [1e-20, 1e-5, 0.1],
            "bitscore": [200.0, 50.0, 20.0],
            "qlen": [500, 500, 500],
            "slen": [3000, 3000, 3000],
            "qseq": ["A"] * 3,
            "sseq": ["A"] * 3,
            "strand": ["plus"] * 3,
        })
        annotated, _ = step_quality_tiers(df)
        tiers = annotated.set_index("qseqid")["tier"]
        # q1: pident=90, length=150 -> strict
        assert tiers["q1"] == "strict"
        # q2: pident=78, length=80 -> moderate
        assert tiers["q2"] == "moderate"
        # q3: pident=60, length=30 -> relaxed
        assert tiers["q3"] == "relaxed"

    def test_edit_stats_columns_added(self, blast_df):
        annotated, _ = step_quality_tiers(blast_df)
        for col in ("mismatch_rate", "gap_rate", "edit_distance"):
            assert col in annotated.columns

    def test_summary_has_expected_columns(self, blast_df):
        _, summary = step_quality_tiers(blast_df)
        for col in ("n_hits", "pct", "mean_pident", "mean_length", "mean_evalue"):
            assert col in summary.columns

    def test_empty_input(self):
        df = pd.DataFrame(columns=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen", "qseq", "sseq", "strand",
        ])
        annotated, summary = step_quality_tiers(df)
        assert len(annotated) == 0
        assert "tier" in annotated.columns


class TestStepPositionalAnalysis:
    def test_returns_dict_with_expected_keys(self, blast_df):
        result = step_positional_analysis(blast_df)
        assert isinstance(result, dict)
        assert "utr_profile" in result
        assert "te_profile" in result
        assert "end_bias" in result

    def test_utr_profile_is_dataframe(self, blast_df):
        result = step_positional_analysis(blast_df)
        profile = result["utr_profile"]
        assert isinstance(profile, pd.DataFrame)
        assert "bin_start" in profile.columns
        assert "bin_end" in profile.columns
        assert "n_hits" in profile.columns
        assert "pct_hits" in profile.columns

    def test_te_profile_is_dataframe(self, blast_df):
        result = step_positional_analysis(blast_df)
        profile = result["te_profile"]
        assert isinstance(profile, pd.DataFrame)
        assert len(profile) == 10  # default n_bins

    def test_end_bias_is_dict_with_expected_keys(self, blast_df):
        result = step_positional_analysis(blast_df)
        end_bias = result["end_bias"]
        assert isinstance(end_bias, dict)
        for key in ("five_prime_pct", "three_prime_pct", "end_ratio", "n_total"):
            assert key in end_bias

    def test_end_bias_n_total_matches_input(self, blast_df):
        result = step_positional_analysis(blast_df)
        # n_total should equal number of non-NaN positions (all 5 rows have qlen > 0)
        assert result["end_bias"]["n_total"] == len(blast_df)

    def test_profile_bins_sum_to_total(self, blast_df):
        result = step_positional_analysis(blast_df)
        total_hits = result["utr_profile"]["n_hits"].sum()
        assert total_hits == len(blast_df)


class TestStepMultiplicityAnalysis:
    def test_returns_expected_keys(self, blast_df):
        result = step_multiplicity_analysis(blast_df)
        assert "multiplicity" in result
        assert "te_breadth" in result

    def test_multiplicity_without_query_to_gene(self, blast_df):
        result = step_multiplicity_analysis(blast_df)
        mult = result["multiplicity"]
        assert "hits_per_query" in mult
        assert "queries_per_te" in mult
        # Without query_to_gene, gene-level keys should be absent
        assert "genes_per_te" not in mult
        assert "hits_per_gene" not in mult

    def test_multiplicity_with_query_to_gene(self, blast_df, query_to_gene):
        result = step_multiplicity_analysis(blast_df, query_to_gene)
        mult = result["multiplicity"]
        assert "hits_per_query" in mult
        assert "queries_per_te" in mult
        assert "genes_per_te" in mult
        assert "hits_per_gene" in mult

    def test_multiplicity_summary_stats_keys(self, blast_df):
        result = step_multiplicity_analysis(blast_df)
        for stat_key in ("median", "mean", "max", "n"):
            assert stat_key in result["multiplicity"]["hits_per_query"]

    def test_te_breadth_is_dataframe(self, blast_df, query_to_gene):
        result = step_multiplicity_analysis(blast_df, query_to_gene)
        breadth = result["te_breadth"]
        assert isinstance(breadth, pd.DataFrame)
        assert "n_hits" in breadth.columns
        assert "n_queries" in breadth.columns
        assert "n_genes" in breadth.columns

    def test_te_breadth_without_query_to_gene(self, blast_df):
        result = step_multiplicity_analysis(blast_df)
        breadth = result["te_breadth"]
        assert isinstance(breadth, pd.DataFrame)
        assert "n_hits" in breadth.columns
        assert "n_genes" not in breadth.columns

    def test_empty_input(self):
        df = pd.DataFrame(columns=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen", "qseq", "sseq", "strand",
        ])
        result = step_multiplicity_analysis(df, query_to_gene={"a": "b"})
        assert result["multiplicity"]["hits_per_query"]["n"] == 0
        assert len(result["te_breadth"]) == 0


class TestStepConservationAnalysis:
    def test_missing_bigwig_raises_file_not_found(self, blast_df):
        regions = pd.DataFrame({
            "region_id": ["utr1", "utr2", "utr3"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [1, 200, 50],
            "end": [500, 600, 400],
            "strand": ["+", "+", "-"],
        })
        with pytest.raises(FileNotFoundError, match="bigWig file not found"):
            step_conservation_analysis(
                df=blast_df,
                regions=regions,
                bigwig_path="/nonexistent/path/phylop.bw",
                tool_path="/nonexistent/path/bigWigAverageOverBed",
            )

    def test_missing_tool_raises_file_not_found(self, blast_df, tmp_path):
        # Create a fake bigwig file so that check passes, but tool is missing
        fake_bw = tmp_path / "phylop.bw"
        fake_bw.write_bytes(b"\x00" * 16)
        regions = pd.DataFrame({
            "region_id": ["utr1", "utr2", "utr3"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [1, 200, 50],
            "end": [500, 600, 400],
            "strand": ["+", "+", "-"],
        })
        with pytest.raises(FileNotFoundError, match="bigWigAverageOverBed not found"):
            step_conservation_analysis(
                df=blast_df,
                regions=regions,
                bigwig_path=fake_bw,
                tool_path="/nonexistent/path/bigWigAverageOverBed",
            )

    def test_empty_df_returns_empty_result(self):
        df = pd.DataFrame(columns=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen", "qseq", "sseq", "strand",
        ])
        regions = pd.DataFrame(columns=["region_id", "chrom", "start", "end", "strand"])
        result = step_conservation_analysis(
            df=df,
            regions=regions,
            bigwig_path="/fake/path.bw",
            tool_path="/fake/tool",
        )
        assert "scores" in result
        assert "scored_df" in result
        assert len(result["scores"]) == 0

    def test_selects_tiers_for_scoring(self, blast_df, tmp_path):
        """Verify that tier filtering logic runs without error when tier column present."""
        # Add tier column to simulate post-quality-tier data
        blast_df["tier"] = ["strict", "moderate", "relaxed", "relaxed", "strict"]
        regions = pd.DataFrame({
            "region_id": ["utr1", "utr2", "utr3"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [1, 200, 50],
            "end": [500, 600, 400],
            "strand": ["+", "+", "-"],
        })
        # Create fake files so we get past the file existence checks
        fake_bw = tmp_path / "phylop.bw"
        fake_bw.write_bytes(b"\x00" * 16)
        fake_tool = tmp_path / "bigWigAverageOverBed"
        fake_tool.write_text("#!/bin/sh\nexit 1\n")
        fake_tool.chmod(0o755)
        # The subprocess call will fail, but we're testing that tier filtering
        # selects the right hits before reaching that point
        with pytest.raises(RuntimeError):
            step_conservation_analysis(
                df=blast_df,
                regions=regions,
                bigwig_path=fake_bw,
                tool_path=fake_tool,
            )
