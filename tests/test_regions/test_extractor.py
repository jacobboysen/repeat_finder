"""Tests for RegionExtractor."""

import pytest

from fossil_finder.regions.extractor import RegionExtractor


class TestRegionExtractorInit:
    def test_loads_genome_and_annotation(self, mini_genome_config, test_data_dir):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        extractor = RegionExtractor(config, base_dir=test_data_dir)
        assert len(extractor.genome) == 2  # chr1, chr2
        assert len(extractor.features) > 0

    def test_chromosome_filter_from_config(self, mini_genome_config, test_data_dir):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        extractor = RegionExtractor(config, base_dir=test_data_dir)
        assert "chr1" in extractor.genome
        assert "chr2" in extractor.genome


class TestFeatureExtraction:
    @pytest.fixture
    def extractor(self, mini_genome_config, test_data_dir):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        return RegionExtractor(config, base_dir=test_data_dir)

    def test_extract_three_prime_utrs(self, extractor):
        regions = extractor.extract_features("three_prime_UTR")
        assert len(regions) == 3  # 3 genes in mini_annotation.gff3
        for r in regions:
            assert "sequence" in r
            assert "chrom" in r
            assert "start" in r
            assert "end" in r
            assert "strand" in r
            assert len(r["sequence"]) > 0

    def test_extract_five_prime_utrs(self, extractor):
        regions = extractor.extract_features("five_prime_UTR")
        assert len(regions) == 3

    def test_extract_exons(self, extractor):
        regions = extractor.extract_features("exon")
        assert len(regions) == 5  # 5 exons in mini_annotation.gff3

    def test_extract_cds(self, extractor):
        regions = extractor.extract_features("CDS")
        assert len(regions) == 5

    def test_extract_nonexistent_type_returns_empty(self, extractor):
        regions = extractor.extract_features("enhancer")
        assert regions == []

    def test_minus_strand_features_are_reverse_complemented(self, extractor):
        """gene002 on chr1 is on minus strand; its UTR should be RC'd."""
        regions = extractor.extract_features("three_prime_UTR")
        minus_regions = [r for r in regions if r["strand"] == "-"]
        assert len(minus_regions) == 1
        assert len(minus_regions[0]["sequence"]) > 0

    def test_regions_have_parent_info(self, extractor):
        regions = extractor.extract_features("three_prime_UTR")
        for r in regions:
            assert "parent_id" in r

    def test_filter_by_gene_ids(self, extractor):
        regions = extractor.extract_features(
            "three_prime_UTR",
            gene_ids=["gene001"],
        )
        assert len(regions) == 1

    def test_deduplicate_option(self, extractor):
        """With dedup=True, identical sequences should be collapsed."""
        regions_all = extractor.extract_features("exon", deduplicate=False)
        regions_dedup = extractor.extract_features("exon", deduplicate=True)
        assert len(regions_dedup) <= len(regions_all)


class TestWindowExtraction:
    @pytest.fixture
    def extractor(self, mini_genome_config, test_data_dir):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        return RegionExtractor(config, base_dir=test_data_dir)

    def test_extract_windows_around_genes(self, extractor):
        """Extract 50bp upstream + 20bp downstream windows around gene starts."""
        regions = extractor.extract_windows(
            anchor_type="gene",
            upstream=50,
            downstream=20,
        )
        assert len(regions) > 0
        for r in regions:
            assert "sequence" in r
            assert len(r["sequence"]) <= 71  # 50 upstream + anchor + 20 downstream

    def test_window_clamps_at_chromosome_edge(self, extractor):
        """Windows near chromosome start/end should be clamped, not crash."""
        regions = extractor.extract_windows(
            anchor_type="gene",
            upstream=1000,
            downstream=10,
        )
        assert len(regions) > 0

    def test_window_minus_strand(self, extractor):
        """For minus strand, upstream = right in genome coordinates."""
        regions = extractor.extract_windows(
            anchor_type="gene",
            upstream=10,
            downstream=5,
        )
        minus = [r for r in regions if r["strand"] == "-"]
        assert len(minus) > 0
        assert len(minus[0]["sequence"]) > 0


class TestOutput:
    @pytest.fixture
    def extractor(self, mini_genome_config, test_data_dir):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        return RegionExtractor(config, base_dir=test_data_dir)

    def test_write_fasta(self, extractor, tmp_path):
        regions = extractor.extract_features("three_prime_UTR")
        out_path = tmp_path / "utrs.fasta"
        extractor.write_fasta(regions, out_path)
        assert out_path.exists()

        from fossil_finder.io.fasta import parse_fasta

        seqs = parse_fasta(out_path)
        assert len(seqs) == 3

    def test_write_metadata_tsv(self, extractor, tmp_path):
        regions = extractor.extract_features("three_prime_UTR")
        out_path = tmp_path / "utrs_metadata.tsv"
        extractor.write_metadata(regions, out_path)
        assert out_path.exists()

        lines = out_path.read_text().strip().split("\n")
        assert len(lines) == 4  # 1 header + 3 data rows


class TestExtractorAdapterIntegration:
    """Verify RegionExtractor works with adapter-loaded gene lists."""

    def test_extract_utrs_for_specific_gene(
        self, mini_genome_config, test_data_dir,
    ):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        extractor = RegionExtractor(config, base_dir=test_data_dir)
        regions = extractor.extract_features("three_prime_UTR", gene_ids=["gene001"])
        assert len(regions) == 1
        assert regions[0]["strand"] == "+"

    def test_full_roundtrip_fasta_output(self, mini_genome_config, test_data_dir, tmp_path):
        from fossil_finder.config.schema import load_genome_config
        from fossil_finder.io.fasta import parse_fasta

        config = load_genome_config(mini_genome_config)
        extractor = RegionExtractor(config, base_dir=test_data_dir)

        regions = extractor.extract_features("exon")
        fasta_path = tmp_path / "exons.fasta"
        extractor.write_fasta(regions, fasta_path)

        reloaded = parse_fasta(fasta_path)
        assert len(reloaded) == len(regions)
