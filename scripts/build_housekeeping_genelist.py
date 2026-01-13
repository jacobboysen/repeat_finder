#!/usr/bin/env python3
"""
Build housekeeping gene list for control comparisons.

Uses well-characterized Drosophila housekeeping genes commonly used
for normalization in expression studies:
Act5C, RpL32 (rp49), Gapdh1, alphaTub84B, Ef1alpha48D, RpS17, etc.

These genes are ubiquitously expressed and not specifically enriched
in germ plasm, serving as negative controls for TE fossil detection.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


# Canonical housekeeping genes with their FlyBase IDs
# Sources: FlyBase, commonly used normalization genes in Drosophila studies
HOUSEKEEPING_GENES = {
    'Act5C': {
        'flybase_id': 'FBgn0000042',
        'full_name': 'Actin 5C',
        'description': 'Cytoskeletal actin, ubiquitously expressed'
    },
    'RpL32': {
        'flybase_id': 'FBgn0002626',
        'full_name': 'Ribosomal protein L32 (rp49)',
        'description': 'Ribosomal protein, standard normalization gene'
    },
    'Gapdh1': {
        'flybase_id': 'FBgn0001091',
        'full_name': 'Glyceraldehyde 3 phosphate dehydrogenase 1',
        'description': 'Glycolytic enzyme, ubiquitously expressed'
    },
    'alphaTub84B': {
        'flybase_id': 'FBgn0003884',
        'full_name': 'alpha-Tubulin at 84B',
        'description': 'Cytoskeletal tubulin, constitutively expressed'
    },
    'Ef1alpha48D': {
        'flybase_id': 'FBgn0000556',
        'full_name': 'Elongation factor 1alpha 48D',
        'description': 'Translation factor, highly expressed'
    },
    'RpS17': {
        'flybase_id': 'FBgn0005533',
        'full_name': 'Ribosomal protein S17',
        'description': 'Ribosomal protein, constitutive expression'
    },
    'Tbp': {
        'flybase_id': 'FBgn0010851',
        'full_name': 'TATA binding protein',
        'description': 'General transcription factor'
    },
    'eIF4A': {
        'flybase_id': 'FBgn0001942',
        'full_name': 'Eukaryotic initiation factor 4A',
        'description': 'Translation initiation factor'
    },
    'SdhA': {
        'flybase_id': 'FBgn0261439',
        'full_name': 'Succinate dehydrogenase A',
        'description': 'Mitochondrial enzyme, metabolic housekeeping'
    },
    'Rpl4': {
        'flybase_id': 'FBgn0003279',
        'full_name': 'Ribosomal protein L4',
        'description': 'Ribosomal protein, constitutive expression'
    }
}


def write_consolidated_tsv(output_file, genes):
    """Write consolidated gene list as TSV."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        # Header (tier=0 for housekeeping/control)
        f.write('gene_symbol\tflybase_id\tfull_name\ttier\tdescription\tsource\n')

        for symbol, info in genes.items():
            f.write(f"{symbol}\t{info['flybase_id']}\t{info['full_name']}\t")
            f.write(f"0\t{info['description']}\thousekeeping\n")

    return len(genes)


def write_fbgn_list(output_file, genes):
    """Write FBgn ID list for filtering."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        for info in genes.values():
            f.write(f"{info['flybase_id']}\n")

    return len(genes)


def write_status_json(output_file, genes_count):
    """Write status JSON with metadata."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    status = {
        'success': True,
        'timestamp': datetime.now().isoformat(),
        'source': 'housekeeping',
        'genes_found': genes_count,
        'notes': 'Canonical Drosophila housekeeping genes for control comparison'
    }

    with open(output_file, 'w') as f:
        json.dump(status, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/gene_lists'),
        help='Output directory (default: data/gene_lists)'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List housekeeping genes and exit'
    )

    args = parser.parse_args()

    # List genes if requested
    if args.list:
        print("\nCanonical Housekeeping Genes:")
        print("=" * 70)
        print(f"{'Symbol':<15} {'FBgn ID':<15} {'Name'}")
        print("-" * 70)

        for symbol, info in sorted(HOUSEKEEPING_GENES.items()):
            print(f"{symbol:<15} {info['flybase_id']:<15} {info['full_name']}")

        return 0

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Building Housekeeping Gene List")
    print("=" * 60)
    print(f"Source: Canonical housekeeping genes ({len(HOUSEKEEPING_GENES)} genes)")
    print(f"Output directory: {args.output_dir}\n")

    # Write consolidated TSV
    consolidated_file = args.output_dir / 'housekeeping_genes_consolidated.tsv'
    genes_count = write_consolidated_tsv(consolidated_file, HOUSEKEEPING_GENES)
    print(f"Wrote {genes_count} genes to: {consolidated_file.name}")

    # Write FBgn ID list for filtering
    fbgn_file = args.output_dir / 'housekeeping_fbgn_ids.txt'
    write_fbgn_list(fbgn_file, HOUSEKEEPING_GENES)
    print(f"Wrote FBgn IDs to: {fbgn_file.name}")

    # Write status JSON
    status_file = args.output_dir / 'housekeeping_status.json'
    write_status_json(status_file, genes_count)
    print(f"Wrote status to: {status_file.name}")

    print("\n" + "=" * 60)
    print("Housekeeping gene list complete!")
    print(f"  Total genes: {genes_count}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
