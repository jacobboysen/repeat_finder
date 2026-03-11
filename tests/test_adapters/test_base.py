"""Tests for GenomeAdapter ABC contract."""

import inspect

import pytest

from fossil_finder.adapters.base import GenomeAdapter


class TestGenomeAdapterABC:
    def test_cannot_instantiate_abc_directly(self):
        """ABC should not be instantiable without implementing all methods."""
        with pytest.raises(TypeError, match="abstract"):
            GenomeAdapter(config=None)

    def test_concrete_subclass_must_implement_all_methods(self):
        """Partial implementation should also fail."""

        class PartialAdapter(GenomeAdapter):
            def load_gene_ids(self, path):
                return []

        with pytest.raises(TypeError, match="abstract"):
            PartialAdapter(config=None)

    def test_complete_subclass_can_instantiate(self):
        """A complete implementation should work."""

        class FullAdapter(GenomeAdapter):
            def load_gene_ids(self, path):
                return []

            def load_gene_id_symbol_map(self):
                return {}

            def load_expression(self):
                return None

            def load_go_annotations(self):
                return None

            def load_gene_groups(self):
                return None

            def load_localization(self):
                return None

            def load_te_metadata(self, path):
                return {}

            def parse_fasta_metadata(self, header):
                return {}

        adapter = FullAdapter(config=None)
        assert adapter is not None
        assert adapter.config is None

    def test_abstract_methods_list_is_complete(self):
        """Verify all 8 abstract methods are declared."""
        abstract_methods = {
            name
            for name, _ in inspect.getmembers(GenomeAdapter)
            if getattr(
                getattr(GenomeAdapter, name, None), "__isabstractmethod__", False
            )
        }
        expected = {
            "load_gene_ids",
            "load_gene_id_symbol_map",
            "load_expression",
            "load_go_annotations",
            "load_gene_groups",
            "load_localization",
            "load_te_metadata",
            "parse_fasta_metadata",
        }
        assert abstract_methods == expected


class TestGetAdapterFactory:
    def test_flybase_adapter_from_config(self, flybase_config):
        from fossil_finder.adapters import get_adapter
        from fossil_finder.adapters.flybase import FlyBaseAdapter

        adapter = get_adapter(flybase_config)
        assert isinstance(adapter, FlyBaseAdapter)

    def test_custom_adapter_from_config(self, custom_config):
        from fossil_finder.adapters import get_adapter
        from fossil_finder.adapters.custom import CustomAdapter

        adapter = get_adapter(custom_config)
        assert isinstance(adapter, CustomAdapter)

    def test_unknown_adapter_raises(self, test_data_dir):
        from fossil_finder.adapters import get_adapter
        from fossil_finder.config.schema import GenomeConfig

        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        # Monkey-patch the adapter name to test error handling
        config.source.adapter = "nonexistent"
        with pytest.raises(ValueError, match="Unknown adapter"):
            get_adapter(config)


class TestAdapterConfigIntegration:
    def test_both_adapters_share_same_interface(self, flybase_config, custom_config):
        """Both adapters return compatible types from all methods."""
        from fossil_finder.adapters import get_adapter

        fb = get_adapter(flybase_config)
        custom = get_adapter(custom_config)

        # Both return dict from symbol map (possibly empty)
        assert isinstance(fb.load_gene_id_symbol_map(), dict)
        assert isinstance(custom.load_gene_id_symbol_map(), dict)

        # Both return None or dict from optional methods
        for method_name in [
            "load_expression", "load_go_annotations",
            "load_gene_groups", "load_localization",
        ]:
            fb_result = getattr(fb, method_name)()
            custom_result = getattr(custom, method_name)()
            assert fb_result is None or isinstance(fb_result, dict)
            assert custom_result is None or isinstance(custom_result, dict)
