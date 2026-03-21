"""Build TE FASTA libraries from RepeatMasker output.

Creates a simple TE instance FASTA by extracting RepeatMasker-annotated
intervals from the genome. This provides a fallback TE database when no
consensus/instance FASTA is supplied.
"""

from pathlib import Path

from fossil_finder.analysis.repeatmasker import parse_repeatmasker_out
from fossil_finder.io.fasta import parse_fasta, write_fasta
from fossil_finder.regions.sequence import SequenceDeduplicator, extract_subsequence


_EXCLUDED_REPEAT_CLASSES = {"Simple_repeat", "Low_complexity", "Satellite"}


def build_te_fasta_from_repeatmasker(
    genome_fasta: str | Path,
    repeatmasker_out: str | Path,
    output_fasta: str | Path,
    deduplicate: bool = True,
    exclude_simple: bool = True,
) -> Path:
    """Create a TE instance FASTA from RepeatMasker .out annotations.

    Args:
        genome_fasta: Path to genome FASTA.
        repeatmasker_out: Path to RepeatMasker .out file.
        output_fasta: Output FASTA path.
        deduplicate: If True, deduplicate identical sequences.
        exclude_simple: If True (default), exclude Simple_repeat,
            Low_complexity, and Satellite entries. These short tandem
            repeats massively inflate BLAST hits without providing
            meaningful TE fossil signal.

    Returns:
        Path to the output FASTA.
    """
    genome_fasta = Path(genome_fasta)
    repeatmasker_out = Path(repeatmasker_out)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    genome = parse_fasta(genome_fasta)
    rm_df = parse_repeatmasker_out(repeatmasker_out)

    if exclude_simple:
        before = len(rm_df)
        # repeat_class can be "Simple_repeat" or "class/family" like "LINE/Jockey"
        top_class = rm_df["repeat_class"].str.split("/").str[0]
        rm_df = rm_df[~top_class.isin(_EXCLUDED_REPEAT_CLASSES)].copy()
        rm_df = rm_df.reset_index(drop=True)
        filtered = before - len(rm_df)
        if filtered:
            import logging
            logging.getLogger(__name__).info(
                "Excluded %d Simple_repeat/Low_complexity/Satellite entries "
                "(%d remaining)", filtered, len(rm_df),
            )

    # Build a mapping from RM chrom names to genome keys.
    # RM may use UCSC names (chr2L) while genome uses bare names (2L) or vice versa.
    genome_keys = set(genome.keys())
    rm_to_genome: dict[str, str] = {}
    for chrom in rm_df["chrom"].unique():
        if chrom in genome_keys:
            rm_to_genome[chrom] = chrom
        elif chrom.startswith("chr") and chrom[3:] in genome_keys:
            rm_to_genome[chrom] = chrom[3:]
        elif f"chr{chrom}" in genome_keys:
            rm_to_genome[chrom] = f"chr{chrom}"

    dedup = SequenceDeduplicator() if deduplicate else None
    sequences: dict[str, str] = {}

    for idx, row in rm_df.iterrows():
        chrom = rm_to_genome.get(row["chrom"])
        if chrom is None:
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
