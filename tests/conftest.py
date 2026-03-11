"""Shared test fixtures for fossil_finder."""

from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir() -> Path:
    """Path to test fixture data."""
    return Path(__file__).parent / "data"


@pytest.fixture
def mini_genome_fasta(test_data_dir) -> Path:
    """Path to synthetic mini genome FASTA."""
    return test_data_dir / "mini_genome.fasta"


@pytest.fixture
def mini_tes_fasta(test_data_dir) -> Path:
    """Path to synthetic mini TE consensus FASTA."""
    return test_data_dir / "mini_tes.fasta"


@pytest.fixture
def mini_annotation_gff(test_data_dir) -> Path:
    """Path to synthetic mini genome GFF3 annotation."""
    return test_data_dir / "mini_annotation.gff3"


@pytest.fixture
def mini_blast_results(test_data_dir) -> Path:
    """Path to pre-computed BLAST results for mini genome."""
    return test_data_dir / "mini_blast_results.tsv"


@pytest.fixture
def mini_genome_config(test_data_dir) -> Path:
    """Path to genome config YAML for mini test genome."""
    return test_data_dir / "mini_genome_config.yaml"


# --- Phase 2 fixtures ---


@pytest.fixture
def mini_gene_list(test_data_dir) -> Path:
    """Path to synthetic mini gene list (FBgn IDs)."""
    return test_data_dir / "mini_gene_list.txt"


@pytest.fixture
def mini_symbol_map(test_data_dir) -> Path:
    """Path to synthetic gene ID -> symbol mapping TSV."""
    return test_data_dir / "mini_symbol_map.tsv"


@pytest.fixture
def mini_expression(test_data_dir) -> Path:
    """Path to synthetic expression matrix TSV."""
    return test_data_dir / "mini_expression.tsv"


@pytest.fixture
def mini_gaf(test_data_dir) -> Path:
    """Path to synthetic GAF file (GO annotations)."""
    return test_data_dir / "mini_gene_association.gaf"


@pytest.fixture
def mini_gene_groups(test_data_dir) -> Path:
    """Path to synthetic gene group data TSV."""
    return test_data_dir / "mini_gene_groups.tsv"


@pytest.fixture
def mini_flyfish(test_data_dir) -> Path:
    """Path to synthetic FlyFISH localization CSV."""
    return test_data_dir / "mini_flyfish.csv"


@pytest.fixture
def mini_te_consensus_fb(test_data_dir) -> Path:
    """Path to synthetic TE consensus FASTA with FlyBase-style headers."""
    return test_data_dir / "mini_te_consensus_fb.fasta"


@pytest.fixture
def mini_repeatmasker(test_data_dir) -> Path:
    """Path to synthetic RepeatMasker .out file."""
    return test_data_dir / "mini_repeatmasker.out"


@pytest.fixture
def flybase_config(test_data_dir):
    """A GenomeConfig wired to FlyBase adapter with test fixture paths."""
    from fossil_finder.config.schema import GenomeConfig

    return GenomeConfig(
        genome={
            "species": "Testus synthetica",
            "assembly": "test_v1",
            "chromosomes": ["chr1", "chr2"],
        },
        source={
            "adapter": "flybase",
            "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
            "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            "gene_id_prefix": "FBgn",
            "transcript_id_prefix": "FBtr",
            "te_consensus": str(test_data_dir / "mini_te_consensus_fb.fasta"),
            "gene_symbol_map": str(test_data_dir / "mini_symbol_map.tsv"),
            "functional": {
                "expression": str(test_data_dir / "mini_expression.tsv"),
                "go_annotations": str(test_data_dir / "mini_gene_association.gaf"),
                "gene_groups": str(test_data_dir / "mini_gene_groups.tsv"),
                "localization": str(test_data_dir / "mini_flyfish.csv"),
            },
        },
    )


@pytest.fixture
def custom_config(test_data_dir):
    """A GenomeConfig wired to custom adapter (no functional annotations)."""
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
            "gene_id_prefix": "GENE",
            "te_consensus": str(test_data_dir / "mini_tes.fasta"),
        },
    )
