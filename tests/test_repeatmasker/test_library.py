"""Tests for RepeatMasker-derived TE library construction."""

from fossil_finder.io.fasta import parse_fasta
from fossil_finder.repeatmasker.library import build_te_fasta_from_repeatmasker


def test_build_te_fasta_from_repeatmasker(mini_genome_fasta, mini_repeatmasker, tmp_path):
    out_fasta = tmp_path / "te_instances.fasta"
    result = build_te_fasta_from_repeatmasker(
        genome_fasta=mini_genome_fasta,
        repeatmasker_out=mini_repeatmasker,
        output_fasta=out_fasta,
        deduplicate=False,
    )
    assert result == out_fasta
    assert out_fasta.exists()
    sequences = parse_fasta(out_fasta)
    assert len(sequences) == 5
