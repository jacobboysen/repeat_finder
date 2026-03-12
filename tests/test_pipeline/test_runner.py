"""Tests for PipelineRunner orchestration."""

import json
from pathlib import Path

import pandas as pd
import pytest

from fossil_finder.pipeline.runner import PipelineResult, PipelineRunner


@pytest.fixture
def pipeline_config(test_data_dir):
    """Config with blast results pre-computed (no BLAST+ needed)."""
    from fossil_finder.config.schema import GenomeConfig

    return GenomeConfig(
        genome={
            "species": "Testus synthetica",
            "assembly": "test_v1",
            "chromosomes": ["chr1", "chr2"],
        },
        source={
            "adapter": "custom",
            "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
            "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            "te_consensus": str(test_data_dir / "mini_tes.fasta"),
        },
        blast={"word_size": 7, "dust": False},
    )


class TestPipelineResult:
    def test_default_fields(self):
        result = PipelineResult()
        assert result.gene_stats is None
        assert result.strand_bias is None
        assert result.family_stats is None
        assert result.enrichment == {}
        assert result.rm_overlap is None
        assert result.blast_hits is None
        assert result.tier_summary is None
        assert result.positional is None
        assert result.multiplicity is None
        assert result.conservation is None

    def test_summary_dict(self):
        result = PipelineResult()
        result.gene_stats = pd.DataFrame(
            {"hit_count": [5, 3]}, index=["g1", "g2"]
        )
        summary = result.summary()
        assert "n_genes_analyzed" in summary
        assert summary["n_genes_analyzed"] == 2


class TestPipelineRunner:
    def test_init_creates_output_dir(self, pipeline_config, tmp_path):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        assert (tmp_path / "output").exists()

    def test_extract_step(self, pipeline_config, tmp_path):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        regions = runner.extract(feature_type="three_prime_UTR")
        assert len(regions) > 0
        assert (tmp_path / "output" / "regions.fa").exists()
        assert (tmp_path / "output" / "regions.tsv").exists()

    def test_analyze_with_precomputed_blast(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        """Full analysis using pre-computed BLAST results (no BLAST+ needed)."""
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert isinstance(result, PipelineResult)
        assert result.gene_stats is not None
        assert len(result.gene_stats) > 0
        assert result.strand_bias is not None
        assert result.family_stats is not None

    def test_analyze_empty_blast(self, pipeline_config, tmp_path):
        """Empty BLAST results should produce empty but valid result."""
        empty_blast = tmp_path / "empty.tsv"
        empty_blast.write_text("")
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=empty_blast,
            query_to_gene={},
        )
        assert isinstance(result, PipelineResult)
        assert result.gene_stats is not None
        assert len(result.gene_stats) == 0

    def test_analyze_with_gene_sets(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        gene_sets = {
            "test_set": {"gene1", "gene2"},
        }
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
            gene_sets=gene_sets,
        )
        assert "test_set" in result.enrichment

    def test_analyze_with_repeatmasker(
        self, pipeline_config, mini_blast_results, mini_repeatmasker, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        query_regions = pd.DataFrame({
            "region_id": ["gene1_utr", "gene2_utr", "gene3_utr"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [281, 350, 381],
            "end": [300, 369, 400],
        })
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
            repeatmasker_path=mini_repeatmasker,
            query_regions=query_regions,
        )
        assert result.rm_overlap is not None
        assert "known" in result.rm_overlap
        assert "novel" in result.rm_overlap

    def test_saves_summary_json(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        summary_path = tmp_path / "output" / "summary.json"
        assert summary_path.exists()
        summary = json.loads(summary_path.read_text())
        assert "n_genes_analyzed" in summary

    def test_rerun_produces_consistent_results(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        q2g = {"gene1_utr": "gene1", "gene2_utr": "gene2",
               "gene3_utr": "gene3"}
        result1 = runner.analyze(blast_results=mini_blast_results,
                                 query_to_gene=q2g)
        result2 = runner.analyze(blast_results=mini_blast_results,
                                 query_to_gene=q2g)
        assert isinstance(result2, PipelineResult)
        assert result2.gene_stats is not None
        assert len(result1.gene_stats) == len(result2.gene_stats)

    def test_analyze_produces_tier_summary(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert result.tier_summary is not None
        assert isinstance(result.tier_summary, pd.DataFrame)
        assert "n_hits" in result.tier_summary.columns

    def test_analyze_produces_positional(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert result.positional is not None
        assert "utr_profile" in result.positional
        assert "te_profile" in result.positional
        assert "end_bias" in result.positional

    def test_analyze_produces_multiplicity(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert result.multiplicity is not None
        assert "multiplicity" in result.multiplicity
        assert "te_breadth" in result.multiplicity

    def test_analyze_produces_motifs(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert result.motifs is not None
        assert "kmer_counts" in result.motifs
        assert "top_motifs" in result.motifs

    def test_motif_summary_in_summary(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        summary = result.summary()
        assert "n_top_motifs" in summary

    def test_motif_files_saved(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert (tmp_path / "output" / "motif_summary.json").exists()
        motif_data = json.loads(
            (tmp_path / "output" / "motif_summary.json").read_text()
        )
        assert "top_motifs" in motif_data
        assert "n_unique_kmers" in motif_data
