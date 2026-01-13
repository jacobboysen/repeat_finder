#!/usr/bin/env python3
"""
Filter genic segment FASTA files by FlyBase Gene Group membership.

Reads gene group annotations from FlyBase and filters CDS, UTR, and intron
sequences to include only genes belonging to specified groups.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def parse_gene_groups(gene_groups_file):
    """
    Parse FlyBase gene_group_data file.

    Expected format (tab-separated):
    FB_group_id  FB_group_symbol  FB_group_name  Parent_FB_group_id  Parent_FB_group_symbol  Group_member_FB_gene_id  Group_member_FB_gene_symbol

    Returns:
        Dictionary mapping gene group names to sets of gene IDs (FBgn...)
    """
    gene_groups_file = Path(gene_groups_file)

    if not gene_groups_file.exists():
        raise FileNotFoundError(f"Gene groups file not found: {gene_groups_file}")

    groups = defaultdict(set)

    # Handle gzipped files
    open_func = gzip.open if gene_groups_file.suffix == '.gz' else open

    print(f"Parsing gene groups from: {gene_groups_file.name}")

    with open_func(gene_groups_file, 'rt') as f:
        # Skip header
        header = f.readline().strip()

        if not header.startswith('##'):
            # If no ## comment, treat first line as header
            if not header.startswith('FB_group'):
                f.seek(0)  # No header, start from beginning

        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            if len(fields) < 7:
                continue

            group_symbol = fields[1]
            group_name = fields[2]
            gene_id = fields[5]

            # Store by both symbol and name
            groups[group_symbol].add(gene_id)
            groups[group_name].add(gene_id)

    print(f"  Loaded {len(groups)} gene groups")
    return dict(groups)


def extract_gene_id(seq_id):
    """
    Extract FlyBase gene ID (FBgn...) from sequence identifier.

    Handles various FASTA header formats from FlyBase.

    Examples:
        FBgn0000008
        >FBgn0000008 type=gene; loc=X:19961297..19969323
        >FBtr0078104 type=transcript; parent=FBgn0000008
    """
    # Try to find FBgn pattern
    parts = seq_id.split()

    for part in parts:
        if 'FBgn' in part:
            # Extract FBgn ID
            start = part.find('FBgn')
            if start != -1:
                # Get FBgn followed by digits
                fbgn = 'FBgn'
                for char in part[start + 4:]:
                    if char.isdigit():
                        fbgn += char
                    else:
                        break
                return fbgn

    # If FBgn not found, try parsing parent= field for transcripts
    for part in parts:
        if 'parent=' in part.lower():
            parent = part.split('=')[1].rstrip(';,')
            if 'FBgn' in parent:
                return extract_gene_id(parent)

    return None


def filter_fasta_by_genes(input_fasta, output_fasta, gene_ids):
    """
    Filter FASTA file to include only sequences from specified genes.

    Args:
        input_fasta: Path to input FASTA file
        output_fasta: Path to output FASTA file
        gene_ids: Set of FlyBase gene IDs to include

    Returns:
        Number of sequences written
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)

    if not input_fasta.exists():
        print(f"  Warning: Input file not found: {input_fasta}")
        return 0

    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Handle gzipped input
    open_func = gzip.open if input_fasta.suffix == '.gz' else open

    seq_count = 0
    total_count = 0

    with open_func(input_fasta, 'rt') as in_f:
        with open(output_fasta, 'w') as out_f:
            for record in SeqIO.parse(in_f, 'fasta'):
                total_count += 1

                # Extract gene ID from sequence header
                gene_id = extract_gene_id(record.description)

                if gene_id and gene_id in gene_ids:
                    SeqIO.write(record, out_f, 'fasta')
                    seq_count += 1

    return seq_count


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--gene-groups-file',
        type=Path,
        default=Path('data/references/gene_groups.tsv'),
        help='FlyBase gene_group_data file (default: data/references/gene_groups.tsv)'
    )
    parser.add_argument(
        '--reference-dir',
        type=Path,
        default=Path('data/references'),
        help='Directory containing reference FASTA files (default: data/references)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/queries'),
        help='Output directory for filtered sequences (default: data/queries)'
    )
    parser.add_argument(
        '--gene-group',
        type=str,
        required=False,
        help='Gene group to filter by (use --list-groups to see available)'
    )
    parser.add_argument(
        '--segment-types',
        nargs='+',
        choices=['CDS', '5UTR', '3UTR', 'intron', 'all'],
        default=['all'],
        help='Segment types to filter (default: all)'
    )
    parser.add_argument(
        '--list-groups',
        action='store_true',
        help='List available gene groups and exit'
    )

    args = parser.parse_args()

    # Parse gene groups
    try:
        gene_groups = parse_gene_groups(args.gene_groups_file)
    except Exception as e:
        print(f"Error parsing gene groups: {e}", file=sys.stderr)
        return 1

    # List groups if requested
    if args.list_groups:
        print("\nAvailable Gene Groups:")
        print("-" * 60)

        # Sort by size
        sorted_groups = sorted(
            gene_groups.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )

        for group_name, gene_ids in sorted_groups:
            print(f"{len(gene_ids):5d} genes  {group_name}")

        return 0

    # Check that gene group is specified
    if not args.gene_group:
        print("Error: --gene-group required (use --list-groups to see available)",
              file=sys.stderr)
        return 1

    # Check that gene group exists
    if args.gene_group not in gene_groups:
        print(f"Error: Gene group '{args.gene_group}' not found", file=sys.stderr)
        print("Use --list-groups to see available groups", file=sys.stderr)
        return 1

    gene_ids = gene_groups[args.gene_group]
    print(f"\nFiltering by gene group: {args.gene_group}")
    print(f"  Genes in group: {len(gene_ids)}\n")

    # Determine segment types to process
    if 'all' in args.segment_types:
        segment_types = ['CDS', '5UTR', '3UTR', 'intron']
    else:
        segment_types = args.segment_types

    # File mapping
    segment_files = {
        'CDS': 'dmel_cds.fasta',
        '5UTR': 'dmel_5utr.fasta',
        '3UTR': 'dmel_3utr.fasta',
        'intron': 'dmel_intron.fasta'
    }

    # Create output directory for this gene group
    group_output_dir = args.output_dir / args.gene_group.replace(' ', '_')
    group_output_dir.mkdir(parents=True, exist_ok=True)

    # Filter each segment type
    total_seqs = 0

    for segment_type in segment_types:
        input_file = args.reference_dir / segment_files[segment_type]
        output_file = group_output_dir / f"{segment_type.lower()}.fasta"

        print(f"Processing {segment_type}...")
        print(f"  Input:  {input_file}")
        print(f"  Output: {output_file}")

        seq_count = filter_fasta_by_genes(input_file, output_file, gene_ids)

        print(f"  Wrote {seq_count} sequences\n")
        total_seqs += seq_count

    # Summary
    print("=" * 60)
    print(f"Filtering complete!")
    print(f"  Gene group: {args.gene_group}")
    print(f"  Output directory: {group_output_dir}")
    print(f"  Total sequences: {total_seqs}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
