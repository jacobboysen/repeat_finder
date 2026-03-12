"""Tests for GFF3 I/O module."""

import pytest

from fossil_finder.io.gff import (
    parse_gff3,
    iter_gff3,
    get_features_by_type,
    get_children,
    get_gene_to_transcripts,
)


class TestParseGff3:
    def test_loads_features(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        assert len(features) > 0

    def test_feature_has_required_fields(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        f = features[0]
        assert "seqid" in f
        assert "type" in f
        assert "start" in f
        assert "end" in f
        assert "strand" in f
        assert "attributes" in f

    def test_coordinates_are_integers(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        for f in features:
            assert isinstance(f["start"], int)
            assert isinstance(f["end"], int)

    def test_skips_comment_lines(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        for f in features:
            assert not f["seqid"].startswith("#")

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_gff3(tmp_path / "nope.gff3")


class TestGetFeaturesByType:
    def test_get_genes(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        genes = get_features_by_type(features, "gene")
        assert len(genes) == 3

    def test_get_utrs(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        utr3s = get_features_by_type(features, "three_prime_UTR")
        assert len(utr3s) == 3

    def test_get_nonexistent_type(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        enhancers = get_features_by_type(features, "enhancer")
        assert enhancers == []


class TestParseGff3FeatureTypes:
    def test_filter_to_genes_only(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff, feature_types={"gene"})
        assert len(features) == 3
        assert all(f["type"] == "gene" for f in features)

    def test_filter_to_utrs(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff, feature_types={"three_prime_UTR"})
        assert len(features) == 3
        assert all(f["type"] == "three_prime_UTR" for f in features)

    def test_filter_to_multiple_types(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff, feature_types={"gene", "mRNA"})
        assert all(f["type"] in {"gene", "mRNA"} for f in features)
        assert len(features) == 6  # 3 genes + 3 mRNAs

    def test_no_filter_returns_all(self, mini_annotation_gff):
        all_features = parse_gff3(mini_annotation_gff)
        filtered = parse_gff3(mini_annotation_gff, feature_types=None)
        assert len(all_features) == len(filtered)


class TestIterGff3:
    def test_yields_all_features(self, mini_annotation_gff):
        features = list(iter_gff3(mini_annotation_gff))
        expected = parse_gff3(mini_annotation_gff)
        assert len(features) == len(expected)

    def test_features_match_parse(self, mini_annotation_gff):
        iter_features = list(iter_gff3(mini_annotation_gff))
        parse_features = parse_gff3(mini_annotation_gff)
        for a, b in zip(iter_features, parse_features):
            assert a == b

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            list(iter_gff3(tmp_path / "nope.gff3"))


class TestGetChildren:
    def test_get_exons_of_mrna(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        children = get_children(features, "mRNA001")
        child_types = {c["type"] for c in children}
        assert "exon" in child_types
        assert "CDS" in child_types

    def test_get_mrnas_of_gene(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        children = get_children(features, "gene001")
        assert len(children) == 1
        assert children[0]["type"] == "mRNA"


class TestGetGeneToTranscripts:
    def test_maps_genes_to_transcripts(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        mapping = get_gene_to_transcripts(features)
        assert "gene001" in mapping
        assert "mRNA001" in mapping["gene001"]
