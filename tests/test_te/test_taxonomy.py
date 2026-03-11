"""Tests for TE taxonomy / class inference."""

import pytest

from fossil_finder.te.taxonomy import infer_te_class, DEFAULT_FAMILY_LISTS


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
