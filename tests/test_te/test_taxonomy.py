"""Tests for TE taxonomy / class inference."""

import pytest

from fossil_finder.te.taxonomy import (
    DEFAULT_FAMILY_LISTS,
    FLYBASE_OVERRIDES,
    infer_te_class,
    parse_consensus_fasta,
    parse_instance_fasta,
    strip_instance_suffix,
)


class TestDfamNotation:
    def test_ltr_gypsy_dfam(self):
        assert infer_te_class("Gypsy-2_DM#LTR/Gypsy") == "LTR"

    def test_line_jockey_dfam(self):
        assert infer_te_class("Jockey-1_DM#LINE/Jockey") == "LINE"

    def test_dna_mariner_dfam(self):
        assert infer_te_class("Mariner-N1_DM#DNA/Tc1-Mariner") == "DNA"

    def test_sine_dfam(self):
        assert infer_te_class("SINE1A#SINE/tRNA") == "SINE"

    def test_helitron_dfam(self):
        assert infer_te_class("Helitron-1_DM#RC/Helitron") == "Helitron"


class TestFamilyNameMatching:
    def test_gypsy_without_dfam_tag(self):
        assert infer_te_class("gypsy12") == "LTR"

    def test_jockey_without_dfam_tag(self):
        assert infer_te_class("jockey3") == "LINE"

    def test_mariner_without_dfam_tag(self):
        assert infer_te_class("mariner2") == "DNA"

    def test_helitron_without_dfam_tag(self):
        assert infer_te_class("DINE-1") == "Helitron"

    def test_unknown_te(self):
        assert infer_te_class("completely_novel_element") == "Unknown"


class TestCustomFamilyLists:
    def test_override_with_custom_families(self):
        """Custom family lists should work for non-standard genomes."""
        custom = {
            "LTR": {"metavirus", "chromovirus"},
            "DNA": {"mutator", "cacta"},
        }
        assert infer_te_class("Metavirus-1", family_lists=custom) == "LTR"
        assert infer_te_class("CACTA_Os", family_lists=custom) == "DNA"
        # Default families should NOT match when custom is provided
        assert infer_te_class("gypsy12", family_lists=custom) == "Unknown"


class TestDefaultFamilyLists:
    def test_all_classes_present(self):
        assert "LTR" in DEFAULT_FAMILY_LISTS
        assert "LINE" in DEFAULT_FAMILY_LISTS
        assert "DNA" in DEFAULT_FAMILY_LISTS
        assert "Helitron" in DEFAULT_FAMILY_LISTS

    def test_common_families_included(self):
        """Cross-species families should be in defaults."""
        assert "gypsy" in DEFAULT_FAMILY_LISTS["LTR"]
        assert "copia" in DEFAULT_FAMILY_LISTS["LTR"]
        assert "jockey" in DEFAULT_FAMILY_LISTS["LINE"]
        assert "mariner" in DEFAULT_FAMILY_LISTS["DNA"]


# --- New tests for consensus FASTA-based classification ---


class TestParseConsensusFasta:
    def test_basic_parsing(self, tmp_path):
        """Parse standard RepBase-format headers."""
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(
            ">gypsy#LTR/Gypsy\nATCG\n"
            ">hobo#DNA/hAT\nGCTA\n"
            ">jockey#LINE/Jockey\nAAAA\n"
        )
        result = parse_consensus_fasta(fasta)
        assert result == {"gypsy": "LTR", "hobo": "DNA", "jockey": "LINE"}

    def test_rc_normalization(self, tmp_path):
        """RC class should be normalized to Helitron."""
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(">DINE-1#RC/Helitron\nATCG\n")
        result = parse_consensus_fasta(fasta)
        assert result["dine-1"] == "Helitron"

    def test_missing_file_raises(self, tmp_path):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_consensus_fasta(tmp_path / "nonexistent.fa")

    def test_empty_file(self, tmp_path):
        """Empty file should return empty dict."""
        fasta = tmp_path / "empty.fa"
        fasta.write_text("")
        assert parse_consensus_fasta(fasta) == {}

    def test_headers_without_hash_ignored(self, tmp_path):
        """Headers missing # delimiter should be skipped."""
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(
            ">no_class_info\nATCG\n"
            ">gypsy#LTR/Gypsy\nGCTA\n"
        )
        result = parse_consensus_fasta(fasta)
        assert len(result) == 1
        assert "gypsy" in result

    def test_case_normalization(self, tmp_path):
        """Family names should be lowercased."""
        fasta = tmp_path / "consensus.fa"
        fasta.write_text(">Gypsy12#LTR/Gypsy\nATCG\n")
        result = parse_consensus_fasta(fasta)
        assert "gypsy12" in result
        assert result["gypsy12"] == "LTR"


class TestParseInstanceFasta:
    def test_basic_parsing(self, tmp_path):
        """Parse FlyBase TE instance headers with name= field."""
        fasta = tmp_path / "te_instances.fa"
        fasta.write_text(
            ">FBti0019256 type=TE; loc=2R:5000..6000; name=invader2{}555; md5=abc\n"
            "ATCGATCG\n"
            ">FBti0019300 type=TE; loc=3L:1000..2000; name=gypsy{}123; md5=def\n"
            "GCTAGCTA\n"
        )
        result = parse_instance_fasta(fasta)
        assert result == {
            "FBti0019256": "invader2{}555",
            "FBti0019300": "gypsy{}123",
        }

    def test_missing_file_raises(self, tmp_path):
        """Should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_instance_fasta(tmp_path / "nonexistent.fa")

    def test_header_without_name_skipped(self, tmp_path):
        """Headers without name= field should be skipped."""
        fasta = tmp_path / "te_instances.fa"
        fasta.write_text(
            ">FBti0019256 type=TE; loc=2R:5000..6000\nATCG\n"
            ">FBti0019300 type=TE; loc=3L:1000..2000; name=gypsy{}1; md5=x\n"
            "GCTA\n"
        )
        result = parse_instance_fasta(fasta)
        assert len(result) == 1
        assert result["FBti0019300"] == "gypsy{}1"

    def test_te_id_extraction(self, tmp_path):
        """TE ID should be extracted from first whitespace-delimited token."""
        fasta = tmp_path / "te_instances.fa"
        fasta.write_text(
            ">FBti9999999 extra info name=roo{}42; end\nACGT\n"
        )
        result = parse_instance_fasta(fasta)
        assert "FBti9999999" in result
        assert result["FBti9999999"] == "roo{}42"


class TestStripInstanceSuffix:
    def test_basic_strip(self):
        assert strip_instance_suffix("invader2{}555") == "invader2"

    def test_no_suffix(self):
        assert strip_instance_suffix("gypsy") == "gypsy"

    def test_complex_suffix(self):
        """Handle FlyBase suffixes like 1360{}Eph[1520]."""
        assert strip_instance_suffix("1360{}Eph[1520]") == "1360"

    def test_empty_after_braces(self):
        assert strip_instance_suffix("roo{}") == "roo"


class TestInferTeClassWithConsensus:
    @pytest.fixture()
    def consensus_map(self):
        return {
            "gypsy": "LTR",
            "hobo": "DNA",
            "jockey": "LINE",
            "doc2-element": "LINE",
            "dine-1": "Helitron",
        }

    def test_exact_match(self, consensus_map):
        """Exact lowercase match in consensus map."""
        assert infer_te_class("gypsy", consensus_map=consensus_map) == "LTR"
        assert infer_te_class("Hobo", consensus_map=consensus_map) == "DNA"

    def test_instance_suffix_stripped(self, consensus_map):
        """Instance suffixes should be stripped before consensus lookup."""
        assert infer_te_class("gypsy{}123", consensus_map=consensus_map) == "LTR"
        assert infer_te_class("jockey{}42", consensus_map=consensus_map) == "LINE"

    def test_element_suffix_match(self, consensus_map):
        """Names should match consensus entries with -element suffix."""
        assert infer_te_class("Doc2", consensus_map=consensus_map) == "LINE"

    def test_flybase_overrides(self, consensus_map):
        """FLYBASE_OVERRIDES should fire when consensus has no match."""
        assert infer_te_class("H", consensus_map=consensus_map) == "DNA"
        assert infer_te_class("het-tag", consensus_map=consensus_map) == "LINE"
        assert infer_te_class("antonia", consensus_map=consensus_map) == "LTR"
        assert infer_te_class("ninja-dsim-like", consensus_map=consensus_map) == "LTR"
        assert infer_te_class("Y", consensus_map=consensus_map) == "LINE"

    def test_fallback_to_pattern_matching(self, consensus_map):
        """Unknown names should fall through to existing pattern strategies."""
        # 'blood' is in DEFAULT_FAMILY_LISTS["LTR"] but not in our
        # consensus_map fixture, so it should fall through to Strategy 2.
        assert infer_te_class("blood5", consensus_map=consensus_map) == "LTR"

    def test_backwards_compatible_without_consensus(self):
        """Without consensus_map, behavior should be identical to before."""
        assert infer_te_class("gypsy12") == "LTR"
        assert infer_te_class("jockey3") == "LINE"
        assert infer_te_class("completely_novel_element") == "Unknown"
        assert infer_te_class("Gypsy-2_DM#LTR/Gypsy") == "LTR"

    def test_consensus_takes_priority_over_pattern(self):
        """Consensus map should win even if pattern matching would disagree."""
        # Contrived: map "gypsy" to DNA instead of LTR
        weird_map = {"gypsy": "DNA"}
        assert infer_te_class("gypsy", consensus_map=weird_map) == "DNA"

    def test_empty_consensus_map_falls_through(self):
        """An empty consensus map should act like no consensus map."""
        assert infer_te_class("gypsy12", consensus_map={}) == "LTR"


class TestFlybaseOverridesConstant:
    def test_expected_entries(self):
        """All five Drosophila-specific overrides should be present."""
        assert FLYBASE_OVERRIDES == {
            "h": "DNA",
            "het-tag": "LINE",
            "antonia": "LTR",
            "ninja-dsim-like": "LTR",
            "y": "LINE",
        }
