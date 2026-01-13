#!/usr/bin/env python3
"""
Generate dinucleotide-shuffled control sequences.

Creates shuffled versions of 3'UTR sequences that preserve dinucleotide
frequencies (Altschul-Erickson algorithm). This is a more biologically
realistic null model than simple random shuffling.

Reference:
Altschul SF, Erickson BW. Significance of nucleotide sequence alignments:
a method for random sequence permutation that preserves dinucleotide
and codon usage. Mol Biol Evol. 1985;2(6):526-38.
"""

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def build_edge_graph(sequence):
    """
    Build directed graph where edges represent dinucleotides.

    Each node is a nucleotide, each edge is a dinucleotide occurrence.
    """
    graph = defaultdict(list)
    seq = str(sequence).upper()

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if 'N' not in dinuc:  # Skip ambiguous bases
            graph[dinuc[0]].append(dinuc[1])

    return graph


def dinucleotide_shuffle(sequence, rng=None):
    """
    Shuffle sequence preserving dinucleotide frequencies.

    Uses an improved Altschul-Erickson algorithm with Hierholzer's method
    for finding Eulerian paths, which handles disconnected components better.

    Args:
        sequence: Input sequence string
        rng: Random number generator (for reproducibility)

    Returns:
        Shuffled sequence string
    """
    if rng is None:
        rng = random.Random()

    seq = str(sequence).upper()

    if len(seq) < 2:
        return seq

    # Build edge list (each dinucleotide is an edge)
    edges = []
    for i in range(len(seq) - 1):
        if seq[i] not in 'ACGT' or seq[i+1] not in 'ACGT':
            continue
        edges.append((seq[i], seq[i+1]))

    if not edges:
        return seq

    # Shuffle edges
    rng.shuffle(edges)

    # Build adjacency list from shuffled edges
    graph = defaultdict(list)
    for src, dst in edges:
        graph[src].append(dst)

    # Find Eulerian path using Hierholzer's algorithm
    # Start from first nucleotide
    start = seq[0]
    if start not in graph or not graph[start]:
        # Find any node with outgoing edges
        for node in graph:
            if graph[node]:
                start = node
                break

    stack = [start]
    path = []

    while stack:
        current = stack[-1]
        if graph[current]:
            next_node = graph[current].pop()
            stack.append(next_node)
        else:
            path.append(stack.pop())

    # Reverse to get correct order
    path = path[::-1]

    return ''.join(path)


def verify_dinucleotide_composition(original, shuffled):
    """Verify that dinucleotide frequencies are preserved."""
    def count_dinucs(seq):
        counts = defaultdict(int)
        seq = str(seq).upper()
        for i in range(len(seq) - 1):
            dinuc = seq[i:i+2]
            if 'N' not in dinuc:
                counts[dinuc] += 1
        return dict(counts)

    orig_counts = count_dinucs(original)
    shuf_counts = count_dinucs(shuffled)

    return orig_counts == shuf_counts


def shuffle_fasta(input_file, output_file, seed=42, n_shuffles=1, verify=True):
    """
    Shuffle all sequences in a FASTA file.

    Args:
        input_file: Input FASTA path
        output_file: Output FASTA path
        seed: Random seed for reproducibility
        n_shuffles: Number of shuffled copies per sequence
        verify: Verify dinucleotide preservation
    """
    rng = random.Random(seed)
    shuffled_records = []
    verification_failures = 0

    print(f"Shuffling sequences from: {input_file}")
    print(f"  Random seed: {seed}")
    print(f"  Shuffles per sequence: {n_shuffles}")

    for record in SeqIO.parse(input_file, 'fasta'):
        for shuffle_idx in range(n_shuffles):
            shuffled_seq = dinucleotide_shuffle(record.seq, rng)

            # Verify if requested
            if verify:
                if not verify_dinucleotide_composition(record.seq, shuffled_seq):
                    verification_failures += 1
                    print(f"  Warning: Dinucleotide mismatch for {record.id}")

            # Create new record
            if n_shuffles > 1:
                new_id = f"{record.id}_shuffled_{shuffle_idx+1}"
            else:
                new_id = f"{record.id}_shuffled"

            shuffled_record = SeqRecord(
                Seq(shuffled_seq),
                id=new_id,
                description=f"dinucleotide_shuffled seed={seed} original={record.id}"
            )
            shuffled_records.append(shuffled_record)

    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(shuffled_records, output_file, 'fasta')

    print(f"\nWrote {len(shuffled_records)} shuffled sequences to: {output_file}")

    if verify and verification_failures > 0:
        print(f"  Warning: {verification_failures} sequences had verification failures")

    return len(shuffled_records), verification_failures


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--input',
        type=Path,
        default=Path('data/queries/germ_plasm/3UTR_sense.fasta'),
        help='Input FASTA file (default: data/queries/germ_plasm/3UTR_sense.fasta)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('data/queries/germ_plasm/3UTR_shuffled.fasta'),
        help='Output FASTA file (default: data/queries/germ_plasm/3UTR_shuffled.fasta)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    parser.add_argument(
        '--n-shuffles',
        type=int,
        default=1,
        help='Number of shuffled copies per sequence (default: 1)'
    )
    parser.add_argument(
        '--no-verify',
        action='store_true',
        help='Skip dinucleotide verification'
    )
    parser.add_argument(
        '--tier1-only',
        action='store_true',
        help='Also shuffle tier1 subset'
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        print("Run extract_germ_plasm_3utrs.py first", file=sys.stderr)
        return 1

    print("Dinucleotide Sequence Shuffling")
    print("=" * 60)

    # Shuffle main file
    n_seqs, failures = shuffle_fasta(
        args.input,
        args.output,
        seed=args.seed,
        n_shuffles=args.n_shuffles,
        verify=not args.no_verify
    )

    # Also shuffle tier1 if requested or if tier1 file exists
    tier1_input = args.input.parent / '3UTR_sense_tier1.fasta'
    if tier1_input.exists():
        tier1_output = args.output.parent / '3UTR_shuffled_tier1.fasta'
        print("\n" + "-" * 60)
        shuffle_fasta(
            tier1_input,
            tier1_output,
            seed=args.seed,
            n_shuffles=args.n_shuffles,
            verify=not args.no_verify
        )

    print("\n" + "=" * 60)
    print("Shuffling complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
