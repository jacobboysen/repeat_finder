#!/usr/bin/env python3
"""
Build consolidated germ plasm gene list using canonical genes.

Uses well-characterized germ plasm genes from the literature:
nanos, oskar, pgc, gcl, vasa, aubergine, tudor, piwi, AGO3, deadhead, CycB, Hsp83

These genes have documented roles in germ plasm assembly, pole cell formation,
or localization to the posterior pole.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


# Canonical germ plasm genes with their FlyBase IDs
# Sources: FlyBase (https://flybase.org)
CANONICAL_GENES = {
    'nos': {
        'flybase_id': 'FBgn0002962',
        'full_name': 'nanos',
        'tier': 1,
        'description': 'Maternal RNA-binding protein, posterior determinant'
    },
    'osk': {
        'flybase_id': 'FBgn0003015',
        'full_name': 'oskar',
        'tier': 1,
        'description': 'Germ plasm organizer, directs nanos localization'
    },
    'pgc': {
        'flybase_id': 'FBgn0016053',
        'full_name': 'polar granule component',
        'tier': 1,
        'description': 'Polar granule structural component'
    },
    'gcl': {
        'flybase_id': 'FBgn0005695',
        'full_name': 'germ cell-less',
        'tier': 1,
        'description': 'BTB protein, transcriptional repression in germline'
    },
    'vas': {
        'flybase_id': 'FBgn0283442',
        'full_name': 'vasa',
        'tier': 1,
        'description': 'DEAD-box RNA helicase, polar granule component'
    },
    'aub': {
        'flybase_id': 'FBgn0000146',
        'full_name': 'aubergine',
        'tier': 1,
        'description': 'Piwi subfamily, piRNA pathway, germ cell mRNA localization'
    },
    'CycB': {
        'flybase_id': 'FBgn0000405',
        'full_name': 'Cyclin B',
        'tier': 1,
        'description': 'Cell cycle regulator, transcript localized to posterior pole'
    },
    'tud': {
        'flybase_id': 'FBgn0003891',
        'full_name': 'tudor',
        'tier': 2,
        'description': 'Tudor domain protein, polar granule assembly scaffold'
    },
    'piwi': {
        'flybase_id': 'FBgn0004872',
        'full_name': 'piwi',
        'tier': 2,
        'description': 'piRNA pathway, transposon silencing in germline'
    },
    'AGO3': {
        'flybase_id': 'FBgn0250816',
        'full_name': 'Argonaute 3',
        'tier': 2,
        'description': 'piRNA pathway, transposon silencing'
    },
    'dhd': {
        'flybase_id': 'FBgn0011761',
        'full_name': 'deadhead',
        'tier': 2,
        'description': 'Thioredoxin, required for female meiosis'
    },
    'Hsp83': {
        'flybase_id': 'FBgn0001233',
        'full_name': 'Heat shock protein 83',
        'tier': 2,
        'description': 'Chaperone, required for piRNA loading into PIWI proteins'
    }
}

# Tier 1 genes for initial testing (well-characterized with known TE content)
TIER1_GENES = ['nos', 'osk', 'pgc', 'gcl', 'vas', 'aub', 'CycB']


def write_consolidated_tsv(output_file, genes):
    """Write consolidated gene list as TSV."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        # Header
        f.write('gene_symbol\tflybase_id\tfull_name\ttier\tdescription\tsource\n')

        for symbol, info in genes.items():
            f.write(f"{symbol}\t{info['flybase_id']}\t{info['full_name']}\t")
            f.write(f"{info['tier']}\t{info['description']}\tcanonical\n")

    return len(genes)


def write_tier1_list(output_file, tier1_genes):
    """Write tier 1 gene list (one gene symbol per line)."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        for gene in tier1_genes:
            f.write(f"{gene}\n")

    return len(tier1_genes)


def write_status_json(output_file, genes_count, tier1_count):
    """Write status JSON with metadata."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    status = {
        'success': True,
        'timestamp': datetime.now().isoformat(),
        'source': 'canonical',
        'genes_found': genes_count,
        'tier1_genes': tier1_count,
        'flyfish_used': False,
        'bdgp_used': False,
        'notes': 'Using canonical germ plasm gene list. Fly-FISH/BDGP scraping skipped.'
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
        help='List canonical genes and exit'
    )

    args = parser.parse_args()

    # List genes if requested
    if args.list:
        print("\nCanonical Germ Plasm Genes:")
        print("=" * 70)
        print(f"{'Symbol':<10} {'FBgn ID':<15} {'Tier':<6} {'Name'}")
        print("-" * 70)

        for symbol, info in sorted(CANONICAL_GENES.items(), key=lambda x: (x[1]['tier'], x[0])):
            print(f"{symbol:<10} {info['flybase_id']:<15} {info['tier']:<6} {info['full_name']}")

        print("\nTier 1 (initial testing):", ', '.join(TIER1_GENES))
        return 0

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Building Germ Plasm Gene List")
    print("=" * 60)
    print(f"Source: Canonical gene list ({len(CANONICAL_GENES)} genes)")
    print(f"Output directory: {args.output_dir}\n")

    # Write consolidated TSV
    consolidated_file = args.output_dir / 'germ_plasm_genes_consolidated.tsv'
    genes_count = write_consolidated_tsv(consolidated_file, CANONICAL_GENES)
    print(f"Wrote {genes_count} genes to: {consolidated_file.name}")

    # Write tier 1 list
    tier1_file = args.output_dir / 'germ_plasm_genes_tier1.txt'
    tier1_count = write_tier1_list(tier1_file, TIER1_GENES)
    print(f"Wrote {tier1_count} tier1 genes to: {tier1_file.name}")

    # Write FBgn ID list for filtering
    fbgn_file = args.output_dir / 'germ_plasm_fbgn_ids.txt'
    with open(fbgn_file, 'w') as f:
        for info in CANONICAL_GENES.values():
            f.write(f"{info['flybase_id']}\n")
    print(f"Wrote FBgn IDs to: {fbgn_file.name}")

    # Write status JSON
    status_file = args.output_dir / 'gene_list_status.json'
    write_status_json(status_file, genes_count, tier1_count)
    print(f"Wrote status to: {status_file.name}")

    print("\n" + "=" * 60)
    print("Gene list build complete!")
    print(f"  Total genes: {genes_count}")
    print(f"  Tier 1 genes: {tier1_count}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
