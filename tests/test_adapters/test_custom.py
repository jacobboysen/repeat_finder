"""Tests for CustomAdapter (generic genome, no database)."""

import pytest

from fossil_finder.adapters.custom import CustomAdapter


class TestCustomGeneIDs:
    def test_load_gene_ids_with_custom_prefix(self, custom_config, tmp_path):
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("GENE001\nGENE002\nOTHER003\n")
        adapter = CustomAdapter(custom_config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["GENE001", "GENE002"]

    def test_load_gene_ids_no_prefix_loads_all(self, test_data_dir, tmp_path):
        """When no gene_id_prefix is configured, load all non-comment lines."""
        from fossil_finder.config.schema import GenomeConfig

        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("# comment\nGENE001\nGENE002\n")
        adapter = CustomAdapter(config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["GENE001", "GENE002"]


class TestCustomOptionalMethods:
    def test_symbol_map_returns_empty(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_gene_id_symbol_map() == {}

    def test_expression_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_expression() is None

    def test_go_annotations_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_go_annotations() is None

    def test_localization_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_localization() is None


class TestCustomTEMetadata:
    def test_load_te_metadata_from_plain_fasta(self, custom_config, mini_tes_fasta):
        adapter = CustomAdapter(custom_config)
        te_meta = adapter.load_te_metadata(mini_tes_fasta)
        assert len(te_meta) == 3
        assert te_meta["TE_LTR1"]["te_class"] == "LTR"
        assert te_meta["TE_DNA1"]["te_class"] == "DNA"
        assert te_meta["TE_LINE1"]["te_class"] == "LINE"

    def test_parse_fasta_metadata_key_value(self, custom_config):
        adapter = CustomAdapter(custom_config)
        header = ">TE001 class=LTR; family=gypsy; length=200"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["class"] == "LTR"
        assert meta["family"] == "gypsy"

    def test_parse_fasta_metadata_dfam_notation(self, custom_config):
        """Dfam-style headers: >Gypsy-2_DM#LTR/Gypsy"""
        adapter = CustomAdapter(custom_config)
        header = ">Gypsy-2_DM#LTR/Gypsy"
        meta = adapter.parse_fasta_metadata(header)
        assert meta.get("class") == "LTR"
        assert meta.get("family") == "Gypsy"
