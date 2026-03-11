"""Tests for RepeatMasker parsing and overlap detection."""

import pandas as pd
import pytest

from fossil_finder.analysis.repeatmasker import (
    parse_repeatmasker_out,
    find_overlaps,
    classify_hits,
)


class TestParseRepeatMaskerOut:
    def test_loads_all_records(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        assert len(df) == 5

    def test_columns_present(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        for col in ["chrom", "start", "end", "strand", "repeat_name", "repeat_class"]:
            assert col in df.columns

    def test_coordinate_values(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        first = df.iloc[0]
        assert first["chrom"] == "chr1"
        assert first["start"] == 50
        assert first["end"] == 200
        assert first["repeat_name"] == "TE_gypsy1"
        assert first["repeat_class"] == "LTR/Gypsy"

    def test_strand_parsing(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        assert df.iloc[0]["strand"] == "+"
        assert df.iloc[1]["strand"] == "-"  # C -> minus

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_repeatmasker_out(tmp_path / "nope.out")


class TestFindOverlaps:
    def test_overlap_detected(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR/Gypsy"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr1"],
            "start": [150], "end": [300],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 1
        assert overlaps.iloc[0]["region_id"] == "utr1"
        assert overlaps.iloc[0]["overlap_bp"] == 51  # 150-200 inclusive

    def test_no_overlap(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr1"],
            "start": [300], "end": [400],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 0

    def test_different_chromosomes_no_overlap(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr2"],
            "start": [100], "end": [200],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 0


class TestClassifyHits:
    def test_known_vs_novel(self):
        blast_hits = pd.DataFrame({
            "qseqid": ["utr1", "utr1", "utr2"],
            "qstart": [50, 250, 10],
            "qend": [150, 350, 60],
            "sseqid": ["TE1", "TE2", "TE3"],
        })
        rm_overlaps = pd.DataFrame({
            "region_id": ["utr1"],
            "rm_start_in_query": [31],  # 1-based relative to query start
            "rm_end_in_query": [181],   # 1-based relative to query start
            "repeat_name": ["TE_gypsy1"],
            "repeat_class": ["LTR/Gypsy"],
        })
        known, novel = classify_hits(blast_hits, rm_overlaps)
        # Hit 1 (qstart=50, qend=150) overlaps RM region [31,181] -> known
        assert len(known) == 1
        assert known.iloc[0]["qseqid"] == "utr1"
        assert known.iloc[0]["qstart"] == 50
        # Hit 2 (qstart=250, qend=350) doesn't overlap -> novel
        # Hit 3 (utr2) has no RM data -> novel
        assert len(novel) == 2
