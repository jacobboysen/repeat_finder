#!/usr/bin/env python3
"""
Extract top BLAST hits and display alignments.

Shows the actual sequence matches between 3'UTRs and TEs,
formatted for easy inspection of potential TE fossils.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# BLAST output columns (17 columns with strand)
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'slen', 'qseq', 'sseq', 'strand'
]


def load_te_names(te_fasta):
    """Load TE names from FASTA file."""
    te_names = {}
    try:
        from Bio import SeqIO
        for record in SeqIO.parse(te_fasta, 'fasta'):
            # Extract name from description
            match = re.search(r'name=([^;]+)', record.description)
            if match:
                te_names[record.id] = match.group(1)
            else:
                te_names[record.id] = record.id
    except ImportError:
        pass
    return te_names


def format_alignment(row, te_names=None, width=60):
    """Format a single BLAST hit as a readable alignment."""
    lines = []

    # Header
    gene = row['qseqid'].split('_')[0]
    te_id = row['sseqid']
    te_name = te_names.get(te_id, te_id) if te_names else te_id

    lines.append("=" * 80)
    lines.append(f"Gene: {gene} ({row['qseqid']})")
    lines.append(f"TE: {te_name} ({te_id})")
    lines.append(f"")
    lines.append(f"Identity: {row['pident']:.1f}%  Length: {row['length']} bp  E-value: {row['evalue']:.2e}  Bitscore: {row['bitscore']:.1f}")
    lines.append(f"Query pos: {row['qstart']}-{row['qend']} of {row['qlen']} bp")
    lines.append(f"Subject pos: {row['sstart']}-{row['send']} of {row['slen']} bp")
    lines.append(f"TE strand: {row['strand']} ({'sense' if row['strand'] == 'plus' else 'antisense'})")
    lines.append("")

    # Alignment
    qseq = str(row['qseq'])
    sseq = str(row['sseq'])

    # Create match line
    match_line = ""
    for q, s in zip(qseq, sseq):
        if q == s:
            match_line += "|"
        elif q == '-' or s == '-':
            match_line += " "
        else:
            match_line += "."

    # Print alignment in chunks
    for i in range(0, len(qseq), width):
        chunk_q = qseq[i:i+width]
        chunk_m = match_line[i:i+width]
        chunk_s = sseq[i:i+width]

        q_pos = row['qstart'] + i
        s_pos = row['sstart'] + i if row['strand'] == 'plus' else row['sstart'] - i

        lines.append(f"Query  {q_pos:>6}  {chunk_q}  {q_pos + len(chunk_q) - 1}")
        lines.append(f"              {chunk_m}")
        lines.append(f"Subj   {s_pos:>6}  {chunk_s}  {s_pos + len(chunk_s) - 1 if row['strand'] == 'plus' else s_pos - len(chunk_s) + 1}")
        lines.append("")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'blast_file',
        type=Path,
        help='BLAST results TSV file'
    )
    parser.add_argument(
        '--te-fasta',
        type=Path,
        default=Path('data/references/dmel_te_flybase.fasta'),
        help='TE database FASTA for names'
    )
    parser.add_argument(
        '--output',
        '-o',
        type=Path,
        help='Output file (default: stdout)'
    )
    parser.add_argument(
        '--top',
        '-n',
        type=int,
        default=50,
        help='Number of top hits to show (default: 50)'
    )
    parser.add_argument(
        '--sort-by',
        choices=['evalue', 'bitscore', 'pident', 'length'],
        default='bitscore',
        help='Sort criterion (default: bitscore)'
    )
    parser.add_argument(
        '--min-pident',
        type=float,
        default=0,
        help='Minimum percent identity filter'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=0,
        help='Minimum alignment length filter'
    )
    parser.add_argument(
        '--gene',
        type=str,
        help='Filter by gene name (e.g., "nos", "osk")'
    )
    parser.add_argument(
        '--te-family',
        type=str,
        help='Filter by TE family/name pattern'
    )
    parser.add_argument(
        '--strand',
        choices=['plus', 'minus'],
        help='Filter by strand'
    )

    args = parser.parse_args()

    # Load BLAST results
    if not args.blast_file.exists():
        print(f"Error: File not found: {args.blast_file}", file=sys.stderr)
        return 1

    # Detect column count
    first_line = open(args.blast_file).readline()
    num_cols = len(first_line.strip().split('\t'))

    if num_cols == 17:
        df = pd.read_csv(args.blast_file, sep='\t', header=None, names=BLAST_COLUMNS)
    else:
        df = pd.read_csv(args.blast_file, sep='\t', header=None,
                        names=BLAST_COLUMNS[:-1])
        df['strand'] = df.apply(
            lambda row: 'plus' if row['sstart'] < row['send'] else 'minus',
            axis=1
        )

    print(f"Loaded {len(df)} BLAST hits", file=sys.stderr)

    # Load TE names
    te_names = {}
    if args.te_fasta.exists():
        te_names = load_te_names(args.te_fasta)
        print(f"Loaded {len(te_names)} TE names", file=sys.stderr)

    # Apply filters
    if args.min_pident > 0:
        df = df[df['pident'] >= args.min_pident]
        print(f"After pident filter: {len(df)} hits", file=sys.stderr)

    if args.min_length > 0:
        df = df[df['length'] >= args.min_length]
        print(f"After length filter: {len(df)} hits", file=sys.stderr)

    if args.gene:
        df = df[df['qseqid'].str.startswith(args.gene)]
        print(f"After gene filter ({args.gene}): {len(df)} hits", file=sys.stderr)

    if args.te_family:
        # Match against TE names
        te_ids_matching = [te_id for te_id, name in te_names.items()
                         if args.te_family.lower() in name.lower()]
        df = df[df['sseqid'].isin(te_ids_matching)]
        print(f"After TE filter ({args.te_family}): {len(df)} hits", file=sys.stderr)

    if args.strand:
        df = df[df['strand'] == args.strand]
        print(f"After strand filter ({args.strand}): {len(df)} hits", file=sys.stderr)

    # Sort
    ascending = args.sort_by == 'evalue'  # Lower is better for e-value
    df = df.sort_values(args.sort_by, ascending=ascending)

    # Take top N
    df = df.head(args.top)

    # Format output
    output_lines = []
    output_lines.append(f"Top {len(df)} BLAST Hits")
    output_lines.append(f"Source: {args.blast_file}")
    output_lines.append(f"Sorted by: {args.sort_by}")
    output_lines.append("")

    # Summary table
    output_lines.append("SUMMARY TABLE")
    output_lines.append("-" * 100)
    output_lines.append(f"{'Gene':<15} {'TE Family':<25} {'%ID':>6} {'Len':>5} {'E-value':>10} {'Strand':>6} {'Query Pos':<15}")
    output_lines.append("-" * 100)

    for _, row in df.iterrows():
        gene = row['qseqid'].split('_')[0]
        te_name = te_names.get(row['sseqid'], row['sseqid'])[:24]
        pos = f"{row['qstart']}-{row['qend']}"
        output_lines.append(
            f"{gene:<15} {te_name:<25} {row['pident']:>6.1f} {row['length']:>5} "
            f"{row['evalue']:>10.2e} {row['strand']:>6} {pos:<15}"
        )

    output_lines.append("")
    output_lines.append("")
    output_lines.append("DETAILED ALIGNMENTS")
    output_lines.append("")

    # Detailed alignments
    for _, row in df.iterrows():
        output_lines.append(format_alignment(row, te_names))

    output_text = "\n".join(output_lines)

    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
        print(f"Wrote {len(df)} alignments to {args.output}", file=sys.stderr)
    else:
        print(output_text)

    return 0


if __name__ == '__main__':
    sys.exit(main())
