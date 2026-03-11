"""Tests for genome configuration schema."""

from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from fossil_finder.config.schema import (
    GenomeConfig,
    load_genome_config,
)


MINIMAL_CONFIG = {
    "genome": {
        "species": "Testus synthetica",
        "assembly": "test_v1",
    },
    "source": {
        "adapter": "custom",
        "genome_fasta": "/path/to/genome.fasta",
        "annotation_gff": "/path/to/annotation.gff3",
    },
}


class TestGenomeConfig:
    def test_minimal_config_loads(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.genome.species == "Testus synthetica"
        assert config.genome.assembly == "test_v1"
        assert config.source.adapter == "custom"

    def test_missing_required_fields_raises(self):
        with pytest.raises(ValidationError):
            GenomeConfig(**{"genome": {"species": "Test"}})

    def test_chromosomes_default_to_auto(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.genome.chromosomes is None

    def test_optional_fields_default_to_none(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.source.te_instances is None
        assert config.source.phylop_bigwig is None
        assert config.source.synteny is None

    def test_blast_params_have_defaults(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.blast.word_size == 7
        assert config.blast.dust is False
        assert config.blast.evalue == 0.001

    def test_gene_sets_optional(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.gene_sets == {}

    def test_full_config_with_gene_sets(self):
        full = {
            **MINIMAL_CONFIG,
            "gene_sets": {
                "target_genes": {
                    "description": "Genes of interest",
                    "genes": "/path/to/gene_ids.txt",
                }
            },
        }
        config = GenomeConfig(**full)
        assert "target_genes" in config.gene_sets
        assert config.gene_sets["target_genes"].description == "Genes of interest"

    def test_adapter_must_be_valid_choice(self):
        bad = {**MINIMAL_CONFIG}
        bad["source"] = {**bad["source"], "adapter": "nonexistent_db"}
        with pytest.raises(ValidationError):
            GenomeConfig(**bad)


class TestBlastSpecExtended:
    def test_max_hsps_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.max_hsps == 10

    def test_program_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.program == "blastn"

    def test_program_override(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec(program="tblastx")
        assert spec.program == "tblastx"

    def test_soft_masking_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.soft_masking is True

    def test_blast_spec_in_config(self, mini_genome_config):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        assert config.blast.program == "blastn"
        assert config.blast.max_hsps == 10


class TestLoadGenomeConfig:
    def test_load_from_yaml_file(self, tmp_path):
        config_path = tmp_path / "test_config.yaml"
        config_path.write_text(yaml.dump(MINIMAL_CONFIG))
        config = load_genome_config(config_path)
        assert config.genome.species == "Testus synthetica"

    def test_load_nonexistent_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_genome_config("/nonexistent/config.yaml")


class TestDmelConfig:
    def test_dmel_config_validates(self):
        dmel_path = (
            Path(__file__).parent.parent.parent
            / "src" / "fossil_finder" / "config" / "genomes" / "dmel_r6.66.yaml"
        )
        config = load_genome_config(dmel_path)
        assert config.genome.species == "Drosophila melanogaster"
        assert config.genome.assembly == "dm6"
        assert config.source.adapter == "flybase"
        assert "2L" in config.genome.chromosomes
        assert config.blast.word_size == 7
        assert config.blast.dust is False


class TestIntegration:
    """Verify config + I/O modules work together."""

    def test_config_points_to_readable_fixtures(self, mini_genome_config, test_data_dir):
        config = load_genome_config(mini_genome_config)

        from fossil_finder.io.fasta import parse_fasta
        from fossil_finder.io.gff import parse_gff3

        # Resolve paths relative to test data dir
        genome_path = test_data_dir / config.source.genome_fasta
        gff_path = test_data_dir / config.source.annotation_gff
        te_path = test_data_dir / config.source.te_consensus

        seqs = parse_fasta(genome_path)
        assert len(seqs) == len(config.genome.chromosomes)

        features = parse_gff3(gff_path)
        assert len(features) > 0

        tes = parse_fasta(te_path)
        assert len(tes) == 3
