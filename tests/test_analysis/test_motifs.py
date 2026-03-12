"""Tests for k-mer motif enrichment analysis."""

import numpy as np
import pandas as pd
import pytest

from fossil_finder.analysis.motifs import (
    _benjamini_hochberg,
    extract_kmers,
    count_kmers_in_sequences,
    count_kmers_from_blast,
    compute_kmer_enrichment,
    find_motif_positions,
    compute_motif_positional_profile,
    compute_motif_gene_set_enrichment,
)


# ---------------------------------------------------------------------------
# TestExtractKmers
# ---------------------------------------------------------------------------


class TestExtractKmers:
    def test_basic_extraction(self):
        """Extract 6-mers from a clean sequence."""
        kmers = extract_kmers("ACGTACGTAC", k=6)
        assert len(kmers) == 5  # 10 - 6 + 1
        assert kmers[0] == "ACGTAC"
        # ACGTACGTAC -> indices 0-5: ACGTAC, 1-6: CGTACG, 2-7: GTACGT,
        #   3-8: TACGTA, 4-9: ACGTAC
        assert kmers[-1] == "ACGTAC"
        assert kmers[2] == "GTACGT"

    def test_gaps_removed(self):
        """Gaps are stripped before extracting k-mers."""
        kmers = extract_kmers("ACG-TAC-GTA-C", k=6)
        # After gap removal: ACGTACGTAC (length 10)
        assert len(kmers) == 5

    def test_ambiguous_bases_skipped(self):
        """K-mers containing non-ACGT bases are excluded."""
        kmers = extract_kmers("ACGTNACGTAC", k=6)
        # Clean: ACGTNACGTAC (11 chars, 6 possible 6-mers)
        # ACGTNA -> has N -> skip
        # CGTNAC -> has N -> skip
        # GTNACG -> has N -> skip
        # TNACGT -> has N -> skip
        # NACGTA -> has N -> skip
        # ACGTAC -> valid
        assert len(kmers) == 1
        assert kmers[0] == "ACGTAC"

    def test_short_sequence(self):
        """Sequence shorter than k returns empty list."""
        kmers = extract_kmers("ACGT", k=6)
        assert kmers == []

    def test_empty_string(self):
        """Empty string returns empty list."""
        kmers = extract_kmers("", k=6)
        assert kmers == []

    def test_custom_k(self):
        """Custom k-mer size works correctly."""
        kmers = extract_kmers("ACGTAC", k=3)
        # 6 - 3 + 1 = 4 kmers: ACG, CGT, GTA, TAC
        assert len(kmers) == 4
        assert kmers == ["ACG", "CGT", "GTA", "TAC"]

    def test_lowercase_uppercased(self):
        """Lowercase input is converted to uppercase."""
        kmers = extract_kmers("acgtac", k=4)
        assert all(kmer == kmer.upper() for kmer in kmers)
        assert kmers[0] == "ACGT"


# ---------------------------------------------------------------------------
# TestCountKmersInSequences
# ---------------------------------------------------------------------------


class TestCountKmersInSequences:
    def test_basic_counting(self):
        """Count k-mers from a single sequence."""
        counts = count_kmers_in_sequences(["ACGTACGT"], k=4)
        # ACGT, CGTA, GTAC, TACG, ACGT -> ACGT appears 2x
        assert counts["ACGT"] == 2
        assert counts["CGTA"] == 1

    def test_multiple_sequences(self):
        """Counts are aggregated across multiple sequences."""
        counts = count_kmers_in_sequences(["ACGTAC", "ACGTAC"], k=4)
        # Each sequence has ACGT once -> total 2
        assert counts["ACGT"] == 2

    def test_empty_list(self):
        """Empty input returns empty dict."""
        counts = count_kmers_in_sequences([], k=6)
        assert counts == {}


# ---------------------------------------------------------------------------
# TestCountKmersFromBlast
# ---------------------------------------------------------------------------


class TestCountKmersFromBlast:
    def test_count_from_dataframe(self):
        """Extract k-mers from BLAST result qseq column."""
        df = pd.DataFrame({"qseq": ["ACGTACGT", "TGCATGCA"]})
        counts = count_kmers_from_blast(df, k=4)
        assert "ACGT" in counts
        assert "TGCA" in counts

    def test_custom_seq_col(self):
        """Use a non-default sequence column."""
        df = pd.DataFrame({"sseq": ["AAAACCCC"], "qseq": ["TTTTGGGG"]})
        counts = count_kmers_from_blast(df, k=4, seq_col="sseq")
        assert "AAAA" in counts
        assert "TTTT" not in counts

    def test_empty_dataframe(self):
        """Empty DataFrame returns empty dict."""
        df = pd.DataFrame(columns=["qseq"])
        counts = count_kmers_from_blast(df, k=4)
        assert counts == {}

    def test_missing_column(self):
        """Missing seq_col returns empty dict."""
        df = pd.DataFrame({"other": ["ACGT"]})
        counts = count_kmers_from_blast(df, k=4, seq_col="qseq")
        assert counts == {}


# ---------------------------------------------------------------------------
# TestComputeKmerEnrichment
# ---------------------------------------------------------------------------


class TestComputeKmerEnrichment:
    def test_enriched_motif(self):
        """A motif enriched in real data should have high z-score and enrichment > 1."""
        real = {"ACGTAC": 100, "TGCATG": 50}
        shuf1 = {"ACGTAC": 10, "TGCATG": 48}
        shuf2 = {"ACGTAC": 12, "TGCATG": 52}
        shuf3 = {"ACGTAC": 11, "TGCATG": 50}

        df = compute_kmer_enrichment(real, [shuf1, shuf2, shuf3], min_count=10)
        acgt_row = df[df["motif"] == "ACGTAC"].iloc[0]
        assert acgt_row["enrichment"] > 1.0
        assert acgt_row["z_score"] > 2.0
        assert acgt_row["p_value"] < 0.05

    def test_depleted_motif(self):
        """A motif depleted in real data should have enrichment < 1."""
        real = {"AAAAAA": 5}
        shuf1 = {"AAAAAA": 50}
        shuf2 = {"AAAAAA": 55}
        shuf3 = {"AAAAAA": 45}

        df = compute_kmer_enrichment(real, [shuf1, shuf2, shuf3], min_count=10)
        row = df[df["motif"] == "AAAAAA"].iloc[0]
        assert row["enrichment"] < 1.0
        assert row["z_score"] < 0

    def test_min_count_filter(self):
        """K-mers below min_count threshold are excluded."""
        real = {"ACGTAC": 2}
        shuf = [{"ACGTAC": 1}]

        df = compute_kmer_enrichment(real, shuf, min_count=10)
        assert len(df) == 0

    def test_single_replicate_std_zero(self):
        """Single shuffled replicate -> std=0 -> z_score uses sign heuristic."""
        real = {"ACGTAC": 100}
        shuf = [{"ACGTAC": 10}]

        df = compute_kmer_enrichment(real, shuf, min_count=10)
        row = df[df["motif"] == "ACGTAC"].iloc[0]
        assert row["shuf_std"] == 0.0
        assert row["z_score"] == 10.0  # sign(positive) * 10

    def test_empty_inputs(self):
        """Empty real counts and empty shuffled list returns empty DataFrame."""
        df = compute_kmer_enrichment({}, [], min_count=0)
        assert len(df) == 0
        assert "motif" in df.columns
        assert "fdr" in df.columns

    def test_fdr_values_capped(self):
        """FDR values should all be <= 1.0."""
        real = {f"ACGT{chr(65 + i)}{chr(67 + i)}": 50 + i * 5 for i in range(10)}
        shuf = [{k: max(1, v - 20) for k, v in real.items()} for _ in range(3)]

        df = compute_kmer_enrichment(real, shuf, min_count=5)
        assert (df["fdr"] <= 1.0).all()
        assert (df["fdr"] >= 0.0).all()

    def test_output_columns(self):
        """Output DataFrame has all required columns."""
        real = {"ACGTAC": 50}
        shuf = [{"ACGTAC": 30}, {"ACGTAC": 35}]
        df = compute_kmer_enrichment(real, shuf, min_count=10)

        expected_cols = {
            "motif", "real_count", "shuf_mean", "shuf_std",
            "enrichment", "log2_enrichment", "z_score", "p_value", "fdr",
        }
        assert expected_cols == set(df.columns)

    def test_sorted_by_pvalue(self):
        """Results are sorted by ascending p_value."""
        real = {"AAAAAA": 100, "CCCCCC": 20, "GGGGGG": 60}
        shuf = [
            {"AAAAAA": 10, "CCCCCC": 18, "GGGGGG": 30},
            {"AAAAAA": 12, "CCCCCC": 22, "GGGGGG": 28},
        ]
        df = compute_kmer_enrichment(real, shuf, min_count=5)
        p_vals = df["p_value"].tolist()
        assert p_vals == sorted(p_vals)


# ---------------------------------------------------------------------------
# TestBenjaminiHochberg
# ---------------------------------------------------------------------------


class TestBenjaminiHochberg:
    def test_basic_correction(self):
        """Simple BH correction on three p-values."""
        fdr = _benjamini_hochberg([0.01, 0.04, 0.06])
        assert len(fdr) == 3
        # rank 1: 0.01*3/1=0.03, rank 2: 0.04*3/2=0.06, rank 3: 0.06*3/3=0.06
        assert fdr[0] == pytest.approx(0.03)
        assert fdr[1] == pytest.approx(0.06)
        assert fdr[2] == pytest.approx(0.06)

    def test_empty(self):
        """Empty input returns empty array."""
        result = _benjamini_hochberg([])
        assert len(result) == 0

    def test_monotonicity(self):
        """FDR values are non-decreasing when sorted by p-value."""
        pvals = [0.005, 0.01, 0.03, 0.5, 0.9]
        fdr = _benjamini_hochberg(pvals)
        sorted_fdr = fdr[np.argsort(pvals)]
        for i in range(len(sorted_fdr) - 1):
            assert sorted_fdr[i] <= sorted_fdr[i + 1] + 1e-12


# ---------------------------------------------------------------------------
# TestFindMotifPositions
# ---------------------------------------------------------------------------


class TestFindMotifPositions:
    def test_basic_positions(self):
        """Find a motif at known positions."""
        positions = find_motif_positions("ACGTACGT", "ACGT")
        assert positions == [0, 4]

    def test_overlapping_occurrences(self):
        """Overlapping matches are all found."""
        positions = find_motif_positions("AAAA", "AA")
        assert positions == [0, 1, 2]

    def test_case_insensitive(self):
        """Search is case-insensitive."""
        positions = find_motif_positions("acgtACGT", "ACGT")
        assert positions == [0, 4]

    def test_no_match(self):
        """Motif not present returns empty list."""
        positions = find_motif_positions("ACGTACGT", "TTTT")
        assert positions == []

    def test_gaps_removed(self):
        """Gaps are removed before searching."""
        positions = find_motif_positions("AC-GT-AC-GT", "ACGT")
        assert positions == [0, 4]


# ---------------------------------------------------------------------------
# TestComputeMotifPositionalProfile
# ---------------------------------------------------------------------------


class TestComputeMotifPositionalProfile:
    def test_basic_profile(self):
        """Positional profile is computed for a motif."""
        df = pd.DataFrame(
            {
                "qseq": ["ACGTACGTAC"],
                "qstart": [0],
                "qlen": [100],
            }
        )
        result = compute_motif_positional_profile(df, ["ACGT"], n_bins=10)
        assert "motif" in result.columns
        assert "bin_start" in result.columns
        assert "count" in result.columns
        assert "density" in result.columns
        assert len(result) == 10
        assert result["count"].sum() == 2  # ACGT at pos 0 and 4

    def test_binning_correct(self):
        """Bin boundaries are evenly spaced in [0, 1]."""
        df = pd.DataFrame(
            {
                "qseq": ["ACGTACGT"],
                "qstart": [0],
                "qlen": [100],
            }
        )
        result = compute_motif_positional_profile(df, ["ACGT"], n_bins=5)
        assert len(result) == 5
        assert result.iloc[0]["bin_start"] == pytest.approx(0.0)
        assert result.iloc[0]["bin_end"] == pytest.approx(0.2)
        assert result.iloc[-1]["bin_end"] == pytest.approx(1.0)

    def test_empty_dataframe(self):
        """Empty DataFrame returns zero-count bins."""
        df = pd.DataFrame(columns=["qseq", "qstart", "qlen"])
        result = compute_motif_positional_profile(df, ["ACGT"], n_bins=10)
        assert len(result) == 10
        assert result["count"].sum() == 0

    def test_density_sums_to_one(self):
        """Density values sum to approximately 1.0 when there are hits."""
        df = pd.DataFrame(
            {
                "qseq": ["ACGTACGTACGTACGT"] * 5,
                "qstart": [10, 20, 30, 40, 50],
                "qlen": [1000] * 5,
            }
        )
        result = compute_motif_positional_profile(df, ["ACGT"], n_bins=10)
        total_density = result["density"].sum()
        if result["count"].sum() > 0:
            assert total_density == pytest.approx(1.0, abs=0.01)


# ---------------------------------------------------------------------------
# TestComputeMotifGeneSetEnrichment
# ---------------------------------------------------------------------------


class TestComputeMotifGeneSetEnrichment:
    def _make_blast_df(self):
        """Create a small BLAST-like DataFrame."""
        return pd.DataFrame(
            {
                "qseqid": ["gene1_utr", "gene2_utr", "gene3_utr", "gene4_utr"],
                "qseq": ["ACGTACGTAC", "TGCATGCATG", "ACGTACGTAC", "GGGGGGGGGG"],
            }
        )

    def test_basic_enrichment(self):
        """Fisher's exact test produces valid output."""
        df = self._make_blast_df()
        query_to_gene = {
            "gene1_utr": "g1",
            "gene2_utr": "g2",
            "gene3_utr": "g3",
            "gene4_utr": "g4",
        }
        gene_sets = {"set_a": {"g1", "g3"}, "set_b": {"g2", "g4"}}

        result = compute_motif_gene_set_enrichment(
            df, ["ACGT"], query_to_gene, gene_sets
        )
        assert "motif" in result.columns
        assert "gene_set" in result.columns
        assert "odds_ratio" in result.columns
        assert "p_value" in result.columns
        assert len(result) == 2  # 1 motif x 2 sets

    def test_motif_present_in_set(self):
        """Motif enriched in gene set should have odds_ratio >= 1."""
        df = self._make_blast_df()
        query_to_gene = {
            "gene1_utr": "g1",
            "gene2_utr": "g2",
            "gene3_utr": "g3",
            "gene4_utr": "g4",
        }
        # ACGT is in gene1 and gene3 (set_a) but not gene2 or gene4 (set_b)
        gene_sets = {"set_a": {"g1", "g3"}}

        result = compute_motif_gene_set_enrichment(
            df, ["ACGT"], query_to_gene, gene_sets
        )
        row = result.iloc[0]
        assert row["n_with_motif"] == 2  # g1 and g3 have ACGT
        assert row["odds_ratio"] >= 1.0

    def test_empty_inputs(self):
        """Empty DataFrame produces NaN results."""
        df = pd.DataFrame(columns=["qseqid", "qseq"])
        result = compute_motif_gene_set_enrichment(
            df, ["ACGT"], {}, {"set_a": {"g1"}}
        )
        assert len(result) == 1
        assert np.isnan(result.iloc[0]["odds_ratio"])

    def test_multiple_motifs_and_sets(self):
        """Multiple motifs x multiple gene sets produces correct row count."""
        df = self._make_blast_df()
        query_to_gene = {
            "gene1_utr": "g1",
            "gene2_utr": "g2",
            "gene3_utr": "g3",
            "gene4_utr": "g4",
        }
        gene_sets = {"set_a": {"g1"}, "set_b": {"g2"}, "set_c": {"g3", "g4"}}
        motifs = ["ACGT", "TGCA"]

        result = compute_motif_gene_set_enrichment(
            df, motifs, query_to_gene, gene_sets
        )
        assert len(result) == 6  # 2 motifs x 3 sets
