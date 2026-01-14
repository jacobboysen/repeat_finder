#!/usr/bin/env python3
"""
Compare 5'UTR and 3'UTR TE content across all gene sets.
"""

from pathlib import Path
from collections import defaultdict


def load_gene_list(path):
    """Load gene list from TSV file."""
    genes = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                symbol = parts[0]
                fbgn = parts[1]
                genes[fbgn] = symbol
    return genes


def load_strand_bias(path):
    """Load strand bias data, aggregated by gene."""
    gene_data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                fbgn = parts[1]
                utr_len = int(parts[2])
                total_hits = int(parts[3])
                sense_hits = int(parts[4])
                anti_hits = int(parts[5])

                if fbgn not in gene_data:
                    gene_data[fbgn] = {
                        'total_hits': 0,
                        'sense_hits': 0,
                        'anti_hits': 0,
                        'total_utr_len': 0
                    }
                gene_data[fbgn]['total_hits'] += total_hits
                gene_data[fbgn]['sense_hits'] += sense_hits
                gene_data[fbgn]['anti_hits'] += anti_hits
                gene_data[fbgn]['total_utr_len'] += utr_len

    # Calculate percentages
    for fbgn, data in gene_data.items():
        total = data['total_hits']
        if total > 0:
            data['sense_pct'] = 100 * data['sense_hits'] / total
        else:
            data['sense_pct'] = 50

    return gene_data


def print_gene_set_table(set_name, genes, utr5_data, utr3_data):
    """Print table for a gene set."""
    print(f"\n{'='*90}")
    print(f"  {set_name.upper()} ({len(genes)} genes)")
    print(f"{'='*90}")

    print(f"\n{'Gene':<12} {'5UTR Hits':>10} {'5UTR %S':>8} {'5UTR Len':>9} "
          f"{'3UTR Hits':>10} {'3UTR %S':>8} {'3UTR Len':>9}")
    print("-" * 90)

    total_5utr_hits = 0
    total_5utr_sense = 0
    total_3utr_hits = 0
    total_3utr_sense = 0
    genes_with_data = 0

    for fbgn, symbol in sorted(genes.items(), key=lambda x: x[1]):
        d5 = utr5_data.get(fbgn, {'total_hits': 0, 'sense_pct': 0, 'total_utr_len': 0, 'sense_hits': 0})
        d3 = utr3_data.get(fbgn, {'total_hits': 0, 'sense_pct': 0, 'total_utr_len': 0, 'sense_hits': 0})

        h5 = d5['total_hits']
        s5 = d5['sense_pct'] if h5 > 0 else 0
        l5 = d5['total_utr_len']
        h3 = d3['total_hits']
        s3 = d3['sense_pct'] if h3 > 0 else 0
        l3 = d3['total_utr_len']

        # Format sense percentage with marker for strong bias
        s5_str = f"{s5:.1f}%" if h5 > 0 else "-"
        s3_str = f"{s3:.1f}%" if h3 > 0 else "-"

        print(f"{symbol:<12} {h5:>10,} {s5_str:>8} {l5:>9,} "
              f"{h3:>10,} {s3_str:>8} {l3:>9,}")

        total_5utr_hits += h5
        total_5utr_sense += d5.get('sense_hits', 0)
        total_3utr_hits += h3
        total_3utr_sense += d3.get('sense_hits', 0)
        if h5 > 0 or h3 > 0:
            genes_with_data += 1

    # Summary row
    print("-" * 90)
    avg_5utr_sense = 100 * total_5utr_sense / total_5utr_hits if total_5utr_hits > 0 else 0
    avg_3utr_sense = 100 * total_3utr_sense / total_3utr_hits if total_3utr_hits > 0 else 0

    print(f"{'TOTAL':<12} {total_5utr_hits:>10,} {avg_5utr_sense:>7.1f}% {'-':>9} "
          f"{total_3utr_hits:>10,} {avg_3utr_sense:>7.1f}% {'-':>9}")

    return {
        'name': set_name,
        'n_genes': len(genes),
        'genes_with_hits': genes_with_data,
        'utr5_hits': total_5utr_hits,
        'utr5_sense_pct': avg_5utr_sense,
        'utr3_hits': total_3utr_hits,
        'utr3_sense_pct': avg_3utr_sense
    }


def main():
    base_dir = Path('/Users/jacobboysen/git_repos/repeat_finder')

    # Load strand bias data
    print("Loading TE content data...")
    utr5_data = load_strand_bias(base_dir / 'results/5utr_analysis/strand_bias_by_utr.tsv')
    utr3_data = load_strand_bias(base_dir / 'results/strand_bias_by_utr.tsv')

    # Load all gene sets
    gene_sets = {}
    gene_list_dir = base_dir / 'data/gene_lists'

    set_files = [
        ('germ_plasm', 'germ_plasm_genes_consolidated.tsv'),
        ('housekeeping', 'housekeeping_genes_consolidated.tsv'),
        ('somatic', 'somatic_genes_consolidated.tsv'),
        ('cleared', 'cleared_genes_consolidated.tsv'),
        ('adult', 'adult_genes_consolidated.tsv'),
    ]

    for set_name, filename in set_files:
        filepath = gene_list_dir / filename
        if filepath.exists():
            gene_sets[set_name] = load_gene_list(filepath)
            print(f"  Loaded {set_name}: {len(gene_sets[set_name])} genes")

    print("\n" + "#" * 90)
    print("  5'UTR vs 3'UTR TE CONTENT - ALL GENE SETS")
    print("#" * 90)

    # Print each gene set
    summaries = []
    for set_name, genes in gene_sets.items():
        summary = print_gene_set_table(set_name, genes, utr5_data, utr3_data)
        summaries.append(summary)

    # Print summary comparison
    print("\n" + "=" * 90)
    print("  SUMMARY COMPARISON ACROSS GENE SETS")
    print("=" * 90)

    print(f"\n{'Gene Set':<15} {'N Genes':>8} {'5UTR Hits':>12} {'5UTR %Sense':>12} "
          f"{'3UTR Hits':>12} {'3UTR %Sense':>12} {'5/3 Ratio':>10}")
    print("-" * 90)

    for s in summaries:
        ratio = s['utr5_hits'] / s['utr3_hits'] if s['utr3_hits'] > 0 else 0
        print(f"{s['name']:<15} {s['n_genes']:>8} {s['utr5_hits']:>12,} {s['utr5_sense_pct']:>11.1f}% "
              f"{s['utr3_hits']:>12,} {s['utr3_sense_pct']:>11.1f}% {ratio:>10.2f}x")

    # Key observations
    print("\n" + "=" * 90)
    print("  KEY OBSERVATIONS")
    print("=" * 90)

    # Find highest/lowest sense bias
    max_5utr_sense = max(summaries, key=lambda x: x['utr5_sense_pct'])
    min_5utr_sense = min(summaries, key=lambda x: x['utr5_sense_pct'])
    max_3utr_sense = max(summaries, key=lambda x: x['utr3_sense_pct'])
    min_3utr_sense = min(summaries, key=lambda x: x['utr3_sense_pct'])

    print(f"""
1. 5'UTR Sense Bias Range:
   - Highest: {max_5utr_sense['name']} ({max_5utr_sense['utr5_sense_pct']:.1f}%)
   - Lowest:  {min_5utr_sense['name']} ({min_5utr_sense['utr5_sense_pct']:.1f}%)

2. 3'UTR Sense Bias Range:
   - Highest: {max_3utr_sense['name']} ({max_3utr_sense['utr3_sense_pct']:.1f}%)
   - Lowest:  {min_3utr_sense['name']} ({min_3utr_sense['utr3_sense_pct']:.1f}%)

3. 5'/3' Hit Ratio by Gene Set:""")

    for s in sorted(summaries, key=lambda x: x['utr5_hits']/x['utr3_hits'] if x['utr3_hits'] > 0 else 0, reverse=True):
        ratio = s['utr5_hits'] / s['utr3_hits'] if s['utr3_hits'] > 0 else 0
        print(f"   - {s['name']}: {ratio:.2f}x (5'UTR has {'more' if ratio > 1 else 'fewer'} hits)")


if __name__ == '__main__':
    main()
