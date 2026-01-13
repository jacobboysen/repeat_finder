#!/usr/bin/env python3
"""
List available FlyBase gene groups with gene counts.

Shows all gene groups from FlyBase with the number of genes in each group.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


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
    group_info = {}  # Store additional info about groups

    # Handle gzipped files
    open_func = gzip.open if gene_groups_file.suffix == '.gz' else open

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

            group_id = fields[0]
            group_symbol = fields[1]
            group_name = fields[2]
            parent_id = fields[3] if len(fields) > 3 and fields[3] else None
            parent_symbol = fields[4] if len(fields) > 4 and fields[4] else None
            gene_id = fields[5]
            gene_symbol = fields[6] if len(fields) > 6 else None

            # Store by both symbol and name
            groups[group_symbol].add(gene_id)
            groups[group_name].add(gene_id)

            # Store group info
            if group_symbol not in group_info:
                group_info[group_symbol] = {
                    'id': group_id,
                    'symbol': group_symbol,
                    'name': group_name,
                    'parent_id': parent_id,
                    'parent_symbol': parent_symbol
                }

    return dict(groups), group_info


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
        '--sort-by',
        choices=['name', 'size'],
        default='size',
        help='Sort groups by name or size (default: size)'
    )
    parser.add_argument(
        '--min-genes',
        type=int,
        default=0,
        help='Show only groups with at least N genes (default: 0)'
    )
    parser.add_argument(
        '--search',
        type=str,
        help='Search for groups containing this text (case-insensitive)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show additional information'
    )

    args = parser.parse_args()

    # Check if file exists
    if not args.gene_groups_file.exists():
        print(f"Error: Gene groups file not found: {args.gene_groups_file}",
              file=sys.stderr)
        print("Run download_references.py to download FlyBase data", file=sys.stderr)
        return 1

    # Parse gene groups
    try:
        groups, group_info = parse_gene_groups(args.gene_groups_file)
    except Exception as e:
        print(f"Error parsing gene groups: {e}", file=sys.stderr)
        return 1

    # Filter by search term if specified
    if args.search:
        search_term = args.search.lower()
        groups = {
            name: genes for name, genes in groups.items()
            if search_term in name.lower()
        }

    # Filter by minimum genes
    if args.min_genes > 0:
        groups = {
            name: genes for name, genes in groups.items()
            if len(genes) >= args.min_genes
        }

    # Sort groups
    if args.sort_by == 'name':
        sorted_groups = sorted(groups.items(), key=lambda x: x[0].lower())
    else:  # sort by size
        sorted_groups = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)

    # Print header
    print("FlyBase Gene Groups")
    print("=" * 80)

    if args.search:
        print(f"Search: {args.search}")

    if args.min_genes > 0:
        print(f"Minimum genes: {args.min_genes}")

    print(f"Total groups: {len(sorted_groups)}")
    print()

    # Print groups
    if args.verbose:
        # Verbose format with more details
        for group_name, gene_ids in sorted_groups:
            print(f"{group_name}")
            print(f"  Genes: {len(gene_ids)}")

            # Try to get additional info
            if group_name in group_info:
                info = group_info[group_name]
                print(f"  ID: {info['id']}")
                if info.get('parent_symbol'):
                    print(f"  Parent: {info['parent_symbol']}")

            print()
    else:
        # Compact format
        # Calculate column width
        max_name_len = max((len(name) for name in groups.keys()), default=0)
        col_width = min(max_name_len + 2, 60)

        for group_name, gene_ids in sorted_groups:
            # Truncate long names
            display_name = group_name
            if len(display_name) > col_width - 2:
                display_name = display_name[:col_width - 5] + "..."

            print(f"{len(gene_ids):5d} genes  {display_name}")

    # Summary
    if sorted_groups:
        total_genes = sum(len(genes) for genes in groups.values())
        avg_genes = total_genes / len(sorted_groups) if sorted_groups else 0

        print()
        print("=" * 80)
        print(f"Summary:")
        print(f"  Groups: {len(sorted_groups)}")
        print(f"  Average genes per group: {avg_genes:.1f}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
