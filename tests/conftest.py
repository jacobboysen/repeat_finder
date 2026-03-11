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
