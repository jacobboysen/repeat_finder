"""Tests for FlyBaseAdapter."""

import pytest

from fossil_finder.adapters.flybase import FlyBaseAdapter


class TestFlyBaseGeneIDs:
    def test_load_gene_ids_simple_list(self, flybase_config, mini_gene_list):
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(mini_gene_list)
        assert ids == [
            "FBgn0000001", "FBgn0000002", "FBgn0000003",
            "FBgn0000004", "FBgn0000005",
        ]

    def test_load_gene_ids_skips_comments_and_blanks(self, flybase_config, tmp_path):
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("# comment\n\nFBgn0000001\nFBgn0000002\n\n")
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["FBgn0000001", "FBgn0000002"]

    def test_load_gene_ids_extracts_from_tsv(self, flybase_config, tmp_path):
        gene_file = tmp_path / "genes.tsv"
        gene_file.write_text("geneA\tFBgn0000001\nFBgn0000002\tgeneB\n")
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(gene_file)
        assert "FBgn0000001" in ids
        assert "FBgn0000002" in ids


class TestFlyBaseSymbolMap:
    def test_load_symbol_map(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        mapping = adapter.load_gene_id_symbol_map()
        assert mapping["FBgn0000001"] == "geneA"
        assert mapping["FBgn0000005"] == "geneE"
        assert len(mapping) == 5

    def test_symbol_map_returns_empty_when_unconfigured(self, custom_config):
        from fossil_finder.adapters.custom import CustomAdapter

        adapter = CustomAdapter(custom_config)
        mapping = adapter.load_gene_id_symbol_map()
        assert mapping == {}


class TestFlyBaseExpression:
    def test_load_expression_matrix_format(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        expr = adapter.load_expression()
        assert expr is not None
        assert "FBgn0000001" in expr
        assert abs(expr["FBgn0000001"]["ovary"] - 125.4) < 0.01
        assert abs(expr["FBgn0000002"]["testis"] - 200.3) < 0.01

    def test_load_expression_gene_count(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        expr = adapter.load_expression()
        assert len(expr) == 5

    def test_load_expression_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig

        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_expression() is None


class TestFlyBaseGO:
    def test_load_go_annotations(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        go = adapter.load_go_annotations()
        assert go is not None
        # FBgn0000001 has 2 annotations (F and P)
        assert len(go["FBgn0000001"]) == 2
        aspects = {a["aspect"] for a in go["FBgn0000001"]}
        assert aspects == {"F", "P"}

    def test_go_skips_negative_annotations(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        go = adapter.load_go_annotations()
        # FBgn0000002 has 1 positive + 1 NOT annotation -> only 1 kept
        assert len(go["FBgn0000002"]) == 1

    def test_go_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig

        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_go_annotations() is None


class TestFlyBaseGeneGroups:
    def test_load_gene_groups(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        groups = adapter.load_gene_groups()
        assert groups is not None
        assert "RIBOSOMAL" in groups
        assert len(groups["RIBOSOMAL"]) == 3
        assert "KINASE" in groups
        assert len(groups["KINASE"]) == 2

    def test_gene_groups_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig

        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_gene_groups() is None


class TestFlyBaseLocalization:
    def test_load_localization(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        loc = adapter.load_localization()
        assert loc is not None
        # geneA has "posterior" and "pole plasm"
        assert "geneA" in loc
        assert "posterior" in loc["geneA"]
        assert "pole plasm" in loc["geneA"]

    def test_localization_skips_na(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        loc = adapter.load_localization()
        # geneD has "NA" -- should be skipped
        assert "geneD" not in loc

    def test_localization_returns_none_when_unconfigured(self, custom_config):
        from fossil_finder.adapters.custom import CustomAdapter

        adapter = CustomAdapter(custom_config)
        assert adapter.load_localization() is None


class TestFlyBaseTEMetadata:
    def test_load_te_metadata(self, flybase_config, mini_te_consensus_fb):
        adapter = FlyBaseAdapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)
        assert len(te_meta) == 3
        assert te_meta["FBte0000001"]["name"] == "gypsy12"
        assert te_meta["FBte0000002"]["name"] == "mariner2"
        assert te_meta["FBte0000003"]["name"] == "jockey3"

    def test_te_metadata_has_class_and_length(self, flybase_config, mini_te_consensus_fb):
        adapter = FlyBaseAdapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)
        for te_id, info in te_meta.items():
            assert "name" in info
            assert "te_class" in info
            assert "length" in info
            assert info["length"] > 0


class TestFlyBaseFastaMetadata:
    def test_parse_flybase_header(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        header = ">FBte0000001 name=gypsy12; class=LTR; family=Gypsy; length=200"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["name"] == "gypsy12"
        assert meta["class"] == "LTR"
        assert meta["family"] == "Gypsy"

    def test_parse_header_with_parent(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        header = ">FBtr0070000 type=three_prime_UTR; loc=2L:100..200; parent=FBgn0031081;"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["type"] == "three_prime_UTR"
        assert meta["parent"] == "FBgn0031081"


class TestFlyBaseIntegration:
    """End-to-end: config -> adapter -> all data loading methods."""

    def test_full_pipeline_load(self, flybase_config):
        from fossil_finder.adapters import get_adapter

        adapter = get_adapter(flybase_config)

        symbols = adapter.load_gene_id_symbol_map()
        assert len(symbols) > 0

        expr = adapter.load_expression()
        assert expr is not None and len(expr) > 0

        go = adapter.load_go_annotations()
        assert go is not None and len(go) > 0

        groups = adapter.load_gene_groups()
        assert groups is not None and len(groups) > 0

        loc = adapter.load_localization()
        assert loc is not None and len(loc) > 0

    def test_te_metadata_feeds_classifier(self, flybase_config, mini_te_consensus_fb):
        from fossil_finder.adapters import get_adapter
        from fossil_finder.te.classifier import classify_te_domain

        adapter = get_adapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)

        for te_id, info in te_meta.items():
            result = classify_te_domain(
                te_id=te_id,
                sstart=1,
                send=info["length"] // 2,
                te_length=info["length"],
                te_class=info["te_class"],
            )
            assert result["te_class"] in ("LTR", "LINE", "DNA", "Helitron", "Unknown")
            assert result["domain"] != ""
