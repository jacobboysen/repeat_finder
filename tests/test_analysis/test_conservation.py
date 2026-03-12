"""Tests for conservation / phyloP scoring module."""

import numpy as np
import pandas as pd
import pytest

from fossil_finder.analysis.conservation import (
    compute_pident_conservation_correlation,
    hits_to_genomic_bed,
    score_with_bigwig,
    summarize_conservation_by_group,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_regions(**overrides) -> pd.DataFrame:
    """Minimal regions DataFrame."""
    defaults = {
        "region_id": ["reg_A", "reg_B"],
        "chrom": ["2L", "3R"],
        "start": [1000, 5000],
        "end": [2000, 6000],
        "strand": ["+", "-"],
    }
    defaults.update(overrides)
    return pd.DataFrame(defaults)


def _make_blast_hits(**overrides) -> pd.DataFrame:
    """Minimal BLAST hits DataFrame."""
    defaults = {
        "qseqid": ["reg_A", "reg_B"],
        "qstart": [10, 20],
        "qend": [50, 80],
        "sseqid": ["TE_gypsy", "TE_copia"],
        "pident": [75.0, 60.0],
        "length": [41, 61],
    }
    defaults.update(overrides)
    return pd.DataFrame(defaults)


# ---------------------------------------------------------------------------
# TestHitsToGenomicBed
# ---------------------------------------------------------------------------

class TestHitsToGenomicBed:
    def test_plus_strand_conversion(self):
        """Plus-strand region: g_start = region_start + qstart - 2."""
        regions = _make_regions(
            region_id=["reg_A"], chrom=["2L"],
            start=[1000], end=[2000], strand=["+"],
        )
        hits = _make_blast_hits(
            qseqid=["reg_A"], qstart=[10], qend=[50],
            sseqid=["TE_gypsy"], pident=[75.0], length=[41],
        )
        result = hits_to_genomic_bed(hits, regions)

        assert len(result) == 1
        row = result.iloc[0]
        # g_start = 1000 + 10 - 2 = 1008
        assert row["g_start"] == 1008
        # g_end = 1000 + 50 - 1 = 1049
        assert row["g_end"] == 1049
        assert row["bed_chrom"] == "chr2L"

    def test_minus_strand_conversion(self):
        """Minus-strand region: g_start = region_end - qend."""
        regions = _make_regions(
            region_id=["reg_B"], chrom=["3R"],
            start=[5000], end=[6000], strand=["-"],
        )
        hits = _make_blast_hits(
            qseqid=["reg_B"], qstart=[20], qend=[80],
            sseqid=["TE_copia"], pident=[60.0], length=[61],
        )
        result = hits_to_genomic_bed(hits, regions)

        assert len(result) == 1
        row = result.iloc[0]
        # g_start = 6000 - 80 = 5920
        assert row["g_start"] == 5920
        # g_end = 6000 - 20 + 1 = 5981
        assert row["g_end"] == 5981
        assert row["bed_chrom"] == "chr3R"

    def test_mixed_strands(self):
        """Both plus and minus strand hits processed together."""
        regions = _make_regions()
        hits = _make_blast_hits()
        result = hits_to_genomic_bed(hits, regions)

        assert len(result) == 2
        # Check that both chrom prefixes are present
        assert set(result["bed_chrom"]) == {"chr2L", "chr3R"}

    def test_missing_region_no_match(self):
        """Hits referencing unknown regions are dropped (inner join)."""
        regions = _make_regions(
            region_id=["reg_X"], chrom=["X"],
            start=[100], end=[200], strand=["+"],
        )
        hits = _make_blast_hits(
            qseqid=["reg_UNKNOWN"], qstart=[1], qend=[10],
            sseqid=["TE"], pident=[50.0], length=[10],
        )
        result = hits_to_genomic_bed(hits, regions)
        assert len(result) == 0

    def test_invalid_coords_filtered(self):
        """Coordinates with g_start < 0 or g_end <= g_start are removed."""
        regions = _make_regions(
            region_id=["reg_A"], chrom=["2L"],
            start=[0], end=[10], strand=["+"],
        )
        # qstart=1, qend=1 -> g_start = 0+1-2 = -1 (invalid)
        hits = _make_blast_hits(
            qseqid=["reg_A"], qstart=[1], qend=[1],
            sseqid=["TE"], pident=[50.0], length=[1],
        )
        result = hits_to_genomic_bed(hits, regions)
        assert len(result) == 0

    def test_empty_dataframe(self):
        """Empty hits DataFrame returns empty with expected columns."""
        regions = _make_regions()
        empty = pd.DataFrame(columns=["qseqid", "qstart", "qend"])
        result = hits_to_genomic_bed(empty, regions)
        assert len(result) == 0
        for col in ("bed_chrom", "g_start", "g_end", "hit_id"):
            assert col in result.columns

    def test_hit_id_uniqueness(self):
        """Every row gets a unique hit_id."""
        regions = _make_regions(
            region_id=["reg_A", "reg_A", "reg_A"],
            chrom=["2L", "2L", "2L"],
            start=[1000, 1000, 1000],
            end=[2000, 2000, 2000],
            strand=["+", "+", "+"],
        )
        hits = _make_blast_hits(
            qseqid=["reg_A", "reg_A", "reg_A"],
            qstart=[10, 20, 30],
            qend=[50, 60, 70],
            sseqid=["TE1", "TE2", "TE3"],
            pident=[75.0, 80.0, 65.0],
            length=[41, 41, 41],
        )
        # regions has duplicate region_id — deduplicate for set_index
        regions_dedup = regions.drop_duplicates(subset="region_id")
        result = hits_to_genomic_bed(hits, regions_dedup)
        assert result["hit_id"].nunique() == len(result)


# ---------------------------------------------------------------------------
# TestScoreWithBigwig
# ---------------------------------------------------------------------------

class TestScoreWithBigwig:
    def test_missing_bigwig_raises(self, tmp_path):
        """FileNotFoundError when bigWig path does not exist."""
        bed = pd.DataFrame({
            "bed_chrom": ["chr2L"],
            "g_start": [100],
            "g_end": [200],
            "hit_id": ["hit_0"],
        })
        with pytest.raises(FileNotFoundError, match="bigWig"):
            score_with_bigwig(
                bed,
                bigwig_path=tmp_path / "nonexistent.bw",
                tool_path=tmp_path / "bigWigAverageOverBed",
            )

    def test_missing_tool_raises(self, tmp_path):
        """FileNotFoundError when tool binary does not exist."""
        # Create a fake bigwig so it passes the first check
        bw = tmp_path / "fake.bw"
        bw.touch()

        bed = pd.DataFrame({
            "bed_chrom": ["chr2L"],
            "g_start": [100],
            "g_end": [200],
            "hit_id": ["hit_0"],
        })
        with pytest.raises(FileNotFoundError, match="bigWigAverageOverBed"):
            score_with_bigwig(
                bed,
                bigwig_path=bw,
                tool_path=tmp_path / "nonexistent_tool",
            )


# ---------------------------------------------------------------------------
# TestSummarizeConservationByGroup
# ---------------------------------------------------------------------------

class TestSummarizeConservationByGroup:
    def test_basic_grouping(self):
        """Groups produce correct summary statistics."""
        scores = pd.DataFrame({
            "tier": ["strict", "strict", "moderate", "moderate"],
            "phylop_mean": [0.5, 1.5, -0.3, 0.2],
        })
        result = summarize_conservation_by_group(scores, "tier")

        assert "strict" in result.index
        assert "moderate" in result.index
        assert result.loc["strict", "n_hits"] == 2
        assert result.loc["strict", "mean_phylop"] == pytest.approx(1.0)
        assert result.loc["strict", "pct_positive"] == pytest.approx(100.0)
        assert result.loc["moderate", "pct_negative"] == pytest.approx(50.0)

    def test_single_group(self):
        """Works with only one group value."""
        scores = pd.DataFrame({
            "category": ["A", "A", "A"],
            "phylop_mean": [0.1, 0.2, 0.3],
        })
        result = summarize_conservation_by_group(scores, "category")

        assert len(result) == 1
        assert result.loc["A", "n_hits"] == 3
        assert result.loc["A", "median_phylop"] == pytest.approx(0.2)

    def test_empty_dataframe(self):
        """Empty input returns empty result with expected columns."""
        result = summarize_conservation_by_group(pd.DataFrame(), "tier")
        assert len(result) == 0
        assert "n_hits" in result.columns
        assert "mean_phylop" in result.columns

    def test_missing_group_column(self):
        """If group_col is absent from DataFrame, return empty."""
        scores = pd.DataFrame({"phylop_mean": [0.1, 0.2]})
        result = summarize_conservation_by_group(scores, "nonexistent_col")
        assert len(result) == 0


# ---------------------------------------------------------------------------
# TestComputePidentConservationCorrelation
# ---------------------------------------------------------------------------

class TestComputePidentConservationCorrelation:
    def test_positive_correlation(self):
        """Monotonically increasing pident/phylop -> rho close to +1."""
        rng = np.random.default_rng(42)
        n = 50
        pident = np.linspace(50, 100, n)
        phylop = pident * 0.01 + rng.normal(0, 0.01, n)
        df = pd.DataFrame({"pident": pident, "phylop_mean": phylop})
        result = compute_pident_conservation_correlation(df)

        assert result is not None
        assert result["rho"] > 0.8
        assert result["p_value"] < 0.01
        assert result["n"] == n

    def test_negative_correlation(self):
        """Monotonically decreasing relationship -> rho close to -1."""
        rng = np.random.default_rng(42)
        n = 50
        pident = np.linspace(50, 100, n)
        phylop = -pident * 0.01 + rng.normal(0, 0.01, n)
        df = pd.DataFrame({"pident": pident, "phylop_mean": phylop})
        result = compute_pident_conservation_correlation(df)

        assert result is not None
        assert result["rho"] < -0.8

    def test_too_few_points_returns_none(self):
        """Fewer than 10 data points -> None."""
        df = pd.DataFrame({
            "pident": [60.0, 70.0, 80.0],
            "phylop_mean": [0.1, 0.2, 0.3],
        })
        assert compute_pident_conservation_correlation(df) is None

    def test_empty_dataframe_returns_none(self):
        """Completely empty input -> None."""
        df = pd.DataFrame(columns=["pident", "phylop_mean"])
        assert compute_pident_conservation_correlation(df) is None

    def test_nan_values_excluded(self):
        """Rows with NaN in either column are dropped before counting."""
        n = 15
        df = pd.DataFrame({
            "pident": list(range(50, 50 + n)),
            "phylop_mean": [0.1 * i for i in range(n)],
        })
        # Inject enough NaNs to drop below 10
        df.loc[:5, "phylop_mean"] = np.nan  # 6 NaNs -> 9 valid < 10
        assert compute_pident_conservation_correlation(df) is None
