"""Build TE FASTA libraries from RepeatMasker output.

Creates a simple TE instance FASTA by extracting RepeatMasker-annotated
intervals from the genome. This provides a fallback TE database when no
consensus/instance FASTA is supplied.
"""

from pathlib import Path

from fossil_finder.analysis.repeatmasker import parse_repeatmasker_out
from fossil_finder.io.fasta import parse_fasta, write_fasta
from fossil_finder.regions.sequence import SequenceDeduplicator, extract_subsequence


def build_te_fasta_from_repeatmasker(
    genome_fasta: str | Path,
    repeatmasker_out: str | Path,
    output_fasta: str | Path,
    deduplicate: bool = True,
) -> Path:
    """Create a TE instance FASTA from RepeatMasker .out annotations.

    Args:
        genome_fasta: Path to genome FASTA.
        repeatmasker_out: Path to RepeatMasker .out file.
        output_fasta: Output FASTA path.
        deduplicate: If True, deduplicate identical sequences.

    Returns:
        Path to the output FASTA.
    """
    genome_fasta = Path(genome_fasta)
    repeatmasker_out = Path(repeatmasker_out)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    genome = parse_fasta(genome_fasta)
    rm_df = parse_repeatmasker_out(repeatmasker_out)

    dedup = SequenceDeduplicator() if deduplicate else None
    sequences: dict[str, str] = {}

    for idx, row in rm_df.iterrows():
        chrom = row["chrom"]
        if chrom not in genome:
            continue
        start = int(row["start"])
        end = int(row["end"])
        strand = row["strand"] if row["strand"] in ("+", "-") else "+"

        seq = extract_subsequence(genome, chrom, start, end, strand)
        if not seq:
            continue

        record_id = f"RM_{idx + 1}"
        header = (
            f"{record_id} repeat_name={row['repeat_name']}; "
            f"class={row['repeat_class']}; "
            f"loc={chrom}:{start}-{end}:{strand}; "
            f"divergence={row.get('divergence', '')}"
        )

        if dedup:
            if not dedup.add(record_id, seq, metadata={"header": header}):
                continue

        sequences[header] = seq

    write_fasta(sequences, output_fasta)
    return output_fasta
