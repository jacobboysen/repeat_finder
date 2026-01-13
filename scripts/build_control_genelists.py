#!/usr/bin/env python3
"""
Build control gene lists for TE fossil analysis.

Creates multiple control groups with different expression patterns:
1. Somatic-specific: muscle, neural, cuticle genes (never in germline)
2. Maternally-cleared: mRNAs degraded in early embryo/pole cells
3. Adult-specific: genes not expressed in embryo

These provide better negative controls than always-on housekeeping genes.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


# Somatic-specific genes - never expressed in germline
# Muscle, neural, eye, cuticle genes
SOMATIC_GENES = {
    # Muscle-specific
    'Mhc': {
        'flybase_id': 'FBgn0086783',
        'full_name': 'Myosin heavy chain',
        'category': 'muscle',
        'description': 'Muscle myosin, flight/jump muscles'
    },
    'Act88F': {
        'flybase_id': 'FBgn0000047',
        'full_name': 'Actin 88F',
        'category': 'muscle',
        'description': 'Flight muscle-specific actin'
    },
    'Tm1': {
        'flybase_id': 'FBgn0003721',
        'full_name': 'Tropomyosin 1',
        'category': 'muscle',
        'description': 'Muscle tropomyosin'
    },
    'Mf': {
        'flybase_id': 'FBgn0002741',
        'full_name': 'Muscle-specific protein 300',
        'category': 'muscle',
        'description': 'Muscle structural protein'
    },
    # Neural-specific
    'elav': {
        'flybase_id': 'FBgn0260400',
        'full_name': 'embryonic lethal abnormal vision',
        'category': 'neural',
        'description': 'Pan-neuronal RNA-binding protein'
    },
    'nSyb': {
        'flybase_id': 'FBgn0013342',
        'full_name': 'neuronal Synaptobrevin',
        'category': 'neural',
        'description': 'Synaptic vesicle protein, neurons only'
    },
    'Syt1': {
        'flybase_id': 'FBgn0004242',
        'full_name': 'Synaptotagmin 1',
        'category': 'neural',
        'description': 'Synaptic calcium sensor'
    },
    'ChAT': {
        'flybase_id': 'FBgn0000303',
        'full_name': 'Choline acetyltransferase',
        'category': 'neural',
        'description': 'Cholinergic neuron marker'
    },
    # Eye-specific
    'Rh1': {
        'flybase_id': 'FBgn0002940',
        'full_name': 'Rhodopsin 1',
        'category': 'eye',
        'description': 'Major rhodopsin, R1-6 photoreceptors'
    },
    'ninaE': {
        'flybase_id': 'FBgn0002940',  # Same as Rh1
        'full_name': 'neither inactivation nor afterpotential E',
        'category': 'eye',
        'description': 'Rhodopsin 1 (ninaE = Rh1)'
    },
    # Cuticle genes
    'Cpr65Av': {
        'flybase_id': 'FBgn0035689',
        'full_name': 'Cuticular protein 65Av',
        'category': 'cuticle',
        'description': 'Larval/pupal cuticle protein'
    },
}

# Maternally-cleared genes - mRNAs degraded during MZT or cleared from pole cells
# Many are Smaug targets or cleared by piRNA pathway
CLEARED_GENES = {
    'Hsp70Aa': {
        'flybase_id': 'FBgn0013275',
        'full_name': 'Heat shock protein 70Aa',
        'category': 'heat_shock',
        'description': 'Stress response, cleared from pole cells'
    },
    'Hsp70Ab': {
        'flybase_id': 'FBgn0013276',
        'full_name': 'Heat shock protein 70Ab',
        'category': 'heat_shock',
        'description': 'Stress response, cleared from pole cells'
    },
    'stg': {
        'flybase_id': 'FBgn0003525',
        'full_name': 'string (cdc25)',
        'category': 'cell_cycle',
        'description': 'Cell cycle, degraded in early embryo'
    },
    'twine': {
        'flybase_id': 'FBgn0003900',
        'full_name': 'twine',
        'category': 'cell_cycle',
        'description': 'Meiotic cdc25, maternal-specific'
    },
    'smaug': {
        'flybase_id': 'FBgn0016070',
        'full_name': 'smaug',
        'category': 'rna_decay',
        'description': 'RNA-binding, drives MZT decay'
    },
    'hb': {
        'flybase_id': 'FBgn0001180',
        'full_name': 'hunchback',
        'category': 'patterning',
        'description': 'Gap gene, maternal cleared by MZT'
    },
    'kni': {
        'flybase_id': 'FBgn0001320',
        'full_name': 'knirps',
        'category': 'patterning',
        'description': 'Gap gene, zygotic only'
    },
    'Kr': {
        'flybase_id': 'FBgn0001325',
        'full_name': 'Kruppel',
        'category': 'patterning',
        'description': 'Gap gene, zygotic expression'
    },
    'run': {
        'flybase_id': 'FBgn0003300',
        'full_name': 'runt',
        'category': 'patterning',
        'description': 'Pair-rule gene'
    },
    'ftz': {
        'flybase_id': 'FBgn0001077',
        'full_name': 'fushi tarazu',
        'category': 'patterning',
        'description': 'Pair-rule, maternal cleared'
    },
}

# Adult-specific genes - not expressed in embryo
# Yolk proteins, accessory gland proteins, adult structures
ADULT_GENES = {
    'Yp1': {
        'flybase_id': 'FBgn0004045',
        'full_name': 'Yolk protein 1',
        'category': 'yolk',
        'description': 'Female fat body, adult only'
    },
    'Yp2': {
        'flybase_id': 'FBgn0004047',
        'full_name': 'Yolk protein 2',
        'category': 'yolk',
        'description': 'Female fat body, adult only'
    },
    'Yp3': {
        'flybase_id': 'FBgn0004048',
        'full_name': 'Yolk protein 3',
        'category': 'yolk',
        'description': 'Female fat body, adult only'
    },
    'Acp26Aa': {
        'flybase_id': 'FBgn0002855',
        'full_name': 'Accessory gland protein 26Aa',
        'category': 'accessory_gland',
        'description': 'Male accessory gland, adult only'
    },
    'Acp36DE': {
        'flybase_id': 'FBgn0011559',
        'full_name': 'Accessory gland protein 36DE',
        'category': 'accessory_gland',
        'description': 'Male accessory gland, adult only'
    },
    'Acp70A': {
        'flybase_id': 'FBgn0013753',
        'full_name': 'Accessory gland protein 70A (Sex peptide)',
        'category': 'accessory_gland',
        'description': 'Sex peptide, male adult only'
    },
    'Lsp1alpha': {
        'flybase_id': 'FBgn0002562',
        'full_name': 'Larval serum protein 1 alpha',
        'category': 'larval_serum',
        'description': 'Larval fat body, pupae storage'
    },
    'Lsp1beta': {
        'flybase_id': 'FBgn0002563',
        'full_name': 'Larval serum protein 1 beta',
        'category': 'larval_serum',
        'description': 'Larval fat body, pupae storage'
    },
    'Lsp2': {
        'flybase_id': 'FBgn0002565',
        'full_name': 'Larval serum protein 2',
        'category': 'larval_serum',
        'description': 'Larval fat body, metamorphosis'
    },
    'Obp99b': {
        'flybase_id': 'FBgn0039685',
        'full_name': 'Odorant-binding protein 99b',
        'category': 'sensory',
        'description': 'Adult antenna, olfaction'
    },
}


def write_gene_list(output_file, genes, list_type):
    """Write gene list as TSV."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        f.write('gene_symbol\tflybase_id\tfull_name\ttier\tdescription\tsource\n')

        for symbol, info in genes.items():
            # Use tier=0 for all controls
            f.write(f"{symbol}\t{info['flybase_id']}\t{info['full_name']}\t")
            f.write(f"0\t{info['description']}\t{list_type}\n")

    return len(genes)


def write_fbgn_list(output_file, genes):
    """Write FBgn ID list."""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Use set to avoid duplicates (ninaE = Rh1)
    fbgn_ids = set(info['flybase_id'] for info in genes.values())

    with open(output_file, 'w') as f:
        for fbgn in sorted(fbgn_ids):
            f.write(f"{fbgn}\n")

    return len(fbgn_ids)


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
        choices=['somatic', 'cleared', 'adult', 'all'],
        help='List specific gene set and exit'
    )

    args = parser.parse_args()

    # List genes if requested
    if args.list:
        gene_sets = {
            'somatic': ('Somatic-specific (never in germline)', SOMATIC_GENES),
            'cleared': ('Maternally-cleared (MZT degraded)', CLEARED_GENES),
            'adult': ('Adult-specific (not in embryo)', ADULT_GENES)
        }

        if args.list == 'all':
            sets_to_show = ['somatic', 'cleared', 'adult']
        else:
            sets_to_show = [args.list]

        for set_name in sets_to_show:
            title, genes = gene_sets[set_name]
            print(f"\n{title}:")
            print("=" * 70)
            print(f"{'Symbol':<15} {'FBgn ID':<15} {'Category':<12} {'Name'}")
            print("-" * 70)

            for symbol, info in sorted(genes.items()):
                cat = info.get('category', 'unknown')
                print(f"{symbol:<15} {info['flybase_id']:<15} {cat:<12} {info['full_name']}")

            print(f"\nTotal: {len(genes)} genes")

        return 0

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Building Control Gene Lists")
    print("=" * 60)

    all_stats = {}

    # Process each gene set
    gene_sets = [
        ('somatic', 'Somatic-specific', SOMATIC_GENES),
        ('cleared', 'Maternally-cleared', CLEARED_GENES),
        ('adult', 'Adult-specific', ADULT_GENES),
    ]

    for set_id, set_name, genes in gene_sets:
        print(f"\n{set_name} genes ({len(genes)} genes):")

        # Write consolidated TSV
        tsv_file = args.output_dir / f'{set_id}_genes_consolidated.tsv'
        count = write_gene_list(tsv_file, genes, set_id)
        print(f"  Wrote: {tsv_file.name}")

        # Write FBgn list
        fbgn_file = args.output_dir / f'{set_id}_fbgn_ids.txt'
        unique_count = write_fbgn_list(fbgn_file, genes)
        print(f"  Wrote: {fbgn_file.name} ({unique_count} unique IDs)")

        all_stats[set_id] = {'total': count, 'unique_ids': unique_count}

    # Write combined status
    status = {
        'timestamp': datetime.now().isoformat(),
        'gene_sets': all_stats,
        'notes': 'Control gene lists for TE fossil analysis'
    }

    status_file = args.output_dir / 'control_lists_status.json'
    with open(status_file, 'w') as f:
        json.dump(status, f, indent=2)
    print(f"\nWrote status to: {status_file.name}")

    print("\n" + "=" * 60)
    print("Control gene lists complete!")
    for set_id, stats in all_stats.items():
        print(f"  {set_id}: {stats['total']} genes ({stats['unique_ids']} unique IDs)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
