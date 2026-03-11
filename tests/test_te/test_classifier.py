"""Tests for TE domain position classifier."""

import pytest

from fossil_finder.te.classifier import (
    classify_te_domain,
    get_relative_position,
    get_domain_description,
    is_coding_domain,
    is_regulatory_domain,
)


class TestRelativePosition:
    def test_start_of_element(self):
        start, end = get_relative_position(1, 50, 1000)
        assert abs(start - 0.001) < 0.01
        assert abs(end - 0.05) < 0.01

    def test_middle_of_element(self):
        start, end = get_relative_position(400, 600, 1000)
        assert abs(start - 0.4) < 0.01
        assert abs(end - 0.6) < 0.01

    def test_reversed_orientation(self):
        """sstart > send means minus strand; positions should still be ordered."""
        start, end = get_relative_position(600, 400, 1000)
        assert start < end


class TestClassifyTEDomain:
    def test_ltr_gag_region(self):
        """Hit at 20-30% of an LTR element should be gag."""
        result = classify_te_domain("gypsy1", 200, 300, 1000, te_class="LTR")
        assert result["te_class"] == "LTR"
        assert result["domain"] == "gag"

    def test_ltr_pol_with_subdomain(self):
        """Hit at 50% of an LTR element should be pol/rt."""
        result = classify_te_domain("copia1", 450, 550, 1000, te_class="LTR")
        assert result["domain"] == "pol"
        assert result["domain_detail"] == "rt"

    def test_line_orf2_rt(self):
        """Hit at 50% of a LINE should be orf2_rt."""
        result = classify_te_domain("jockey1", 450, 550, 1000, te_class="LINE")
        assert result["domain"] == "orf2_rt"

    def test_dna_transposase(self):
        """Hit at 50% of a DNA transposon should be transposase."""
        result = classify_te_domain("mariner1", 400, 600, 1000, te_class="DNA")
        assert result["domain"] == "transposase"

    def test_infers_class_from_name(self):
        """When te_class is not provided, infer from name."""
        result = classify_te_domain("gypsy12", 200, 300, 1000)
        assert result["te_class"] == "LTR"

    def test_zero_length_te_does_not_crash(self):
        """te_length=0 should not raise ZeroDivisionError."""
        result = classify_te_domain("gypsy1", 0, 0, 0, te_class="LTR")
        assert result["te_class"] == "LTR"
        # Midpoint 0.0 falls in 5_ltr (0.0-0.12); key thing is no ZeroDivisionError
        assert result["domain"] is not None

    def test_sine_has_no_domain_map(self):
        """SINE elements have no domain map, so domain should be 'unknown'."""
        result = classify_te_domain("SINE1#SINE/tRNA", 50, 100, 200, te_class="SINE")
        assert result["te_class"] == "SINE"
        assert result["domain"] == "unknown"


class TestRelativePositionEdgeCases:
    def test_zero_length_returns_zeros(self):
        start, end = get_relative_position(10, 20, 0)
        assert start == 0.0
        assert end == 0.0

    def test_negative_length_returns_zeros(self):
        start, end = get_relative_position(10, 20, -5)
        assert start == 0.0
        assert end == 0.0


class TestDomainHelpers:
    def test_coding_domain_gag(self):
        assert is_coding_domain("gag") is True

    def test_regulatory_domain_ltr(self):
        assert is_regulatory_domain("5_ltr") is True

    def test_domain_description_exists(self):
        desc = get_domain_description("rt")
        assert "Reverse transcriptase" in desc
