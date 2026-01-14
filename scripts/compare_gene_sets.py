#!/usr/bin/env python3
"""
Compare original gene sets against genome-wide TE metrics.
"""

from pathlib import Path
from collections import defaultdict
import statistics


def load_gene_list(path):
    """Load FBgn IDs from a file."""
    genes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and line.startswith('FBgn'):
                genes.append(line)
    return genes


def load_gene_te_data(path):
    """Load gene-level TE density data."""
    gene_data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                rank = int(parts[0])
                fbgn = parts[1]
                density = float(parts[2])
                hits = int(parts[3])
                hit_bp = int(parts[4])
                utr_len = int(parts[5])
                gene_data[fbgn] = {
                    'rank': rank,
                    'density': density,
                    'hits': hits,
                    'hit_bp': hit_bp,
                    'utr_len': utr_len
                }
    return gene_data


def load_strand_bias_by_utr(path):
    """Load UTR-level strand bias data."""
    utr_data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                fbtr = parts[0]
                fbgn = parts[1]
                total_hits = int(parts[3])
                sense_hits = int(parts[4])
                anti_hits = int(parts[5])
                sense_pct = float(parts[6])

                if fbgn not in utr_data:
                    utr_data[fbgn] = {
                        'total_hits': 0,
                        'sense_hits': 0,
                        'anti_hits': 0
                    }
                utr_data[fbgn]['total_hits'] += total_hits
                utr_data[fbgn]['sense_hits'] += sense_hits
                utr_data[fbgn]['anti_hits'] += anti_hits

    # Calculate percentages
    for fbgn, data in utr_data.items():
        total = data['total_hits']
        if total > 0:
            data['sense_pct'] = 100 * data['sense_hits'] / total
            data['anti_pct'] = 100 * data['anti_hits'] / total
        else:
            data['sense_pct'] = 50
            data['anti_pct'] = 50

    return utr_data


def load_gene_symbols(consolidated_path):
    """Load gene symbols from consolidated TSV."""
    symbols = {}
    if not Path(consolidated_path).exists():
        return symbols
    with open(consolidated_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Format varies, try to find FBgn and symbol
                for i, part in enumerate(parts):
                    if part.startswith('FBgn'):
                        # Symbol is usually first or second column
                        if i > 0:
                            symbols[part] = parts[0]
                        elif len(parts) > 1:
                            symbols[part] = parts[1]
                        break
    return symbols


def main():
    base_dir = Path('/Users/jacobboysen/git_repos/repeat_finder')

    # Load all gene-level TE data (combine top and bottom 100 files with full dataset)
    # First, let's create a complete ranking from the strand bias file
    strand_data = load_strand_bias_by_utr(base_dir / 'results/strand_bias_by_utr.tsv')

    # Load top/bottom 100 for density ranking
    top_100 = load_gene_te_data(base_dir / 'results/top_100_te_genes_FIXED.tsv')
    bottom_100 = load_gene_te_data(base_dir / 'results/bottom_100_te_genes_FIXED.tsv')

    # Merge
    all_gene_data = {**top_100, **bottom_100}

    # Gene sets to analyze
    gene_sets = {
        'germ_plasm': base_dir / 'data/gene_lists/germ_plasm_fbgn_ids.txt',
        'housekeeping': base_dir / 'data/gene_lists/housekeeping_fbgn_ids.txt',
        'somatic': base_dir / 'data/gene_lists/somatic_fbgn_ids.txt',
        'cleared': base_dir / 'data/gene_lists/cleared_fbgn_ids.txt',
        'adult': base_dir / 'data/gene_lists/adult_fbgn_ids.txt',
    }

    # Symbol files
    symbol_files = {
        'germ_plasm': base_dir / 'data/gene_lists/germ_plasm_genes_consolidated.tsv',
        'housekeeping': base_dir / 'data/gene_lists/housekeeping_genes_consolidated.tsv',
        'somatic': base_dir / 'data/gene_lists/somatic_genes_consolidated.tsv',
        'cleared': base_dir / 'data/gene_lists/cleared_genes_consolidated.tsv',
        'adult': base_dir / 'data/gene_lists/adult_genes_consolidated.tsv',
    }

    # Total genes for percentile calculation
    total_genes = 13298

    print("=" * 80)
    print("GENE SET COMPARISON: TE DENSITY & STRAND BIAS")
    print("=" * 80)

    # Summary statistics by group
    group_stats = {}

    for group_name, gene_list_path in gene_sets.items():
        genes = load_gene_list(gene_list_path)
        symbols = load_gene_symbols(symbol_files.get(group_name, ''))

        print(f"\n{'='*60}")
        print(f"## {group_name.upper()} ({len(genes)} genes)")
        print('='*60)

        densities = []
        ranks = []
        sense_pcts = []

        print(f"\n{'Gene':<15} {'Symbol':<12} {'Rank':>8} {'%ile':>8} {'Density':>12} {'Sense%':>8}")
        print("-" * 70)

        for fbgn in genes:
            symbol = symbols.get(fbgn, fbgn[-6:])

            # Get density rank
            if fbgn in all_gene_data:
                rank = all_gene_data[fbgn]['rank']
                density = all_gene_data[fbgn]['density']
                percentile = 100 * (total_genes - rank + 1) / total_genes
            else:
                # Not in top/bottom 100, estimate from strand data
                rank = '—'
                density = None
                percentile = None

            # Get strand bias
            if fbgn in strand_data:
                sense_pct = strand_data[fbgn]['sense_pct']
                total_hits = strand_data[fbgn]['total_hits']
            else:
                sense_pct = None
                total_hits = 0

            # Format output
            rank_str = f"{rank}" if isinstance(rank, int) else rank
            pct_str = f"{percentile:.1f}%" if percentile else "—"
            dens_str = f"{density:.0f}" if density else "—"
            sense_str = f"{sense_pct:.1f}%" if sense_pct else "—"

            print(f"{fbgn:<15} {symbol:<12} {rank_str:>8} {pct_str:>8} {dens_str:>12} {sense_str:>8}")

            if density:
                densities.append(density)
                ranks.append(rank)
            if sense_pct:
                sense_pcts.append(sense_pct)

        # Group summary
        if densities:
            avg_density = statistics.mean(densities)
            avg_rank = statistics.mean(ranks)
            avg_percentile = 100 * (total_genes - avg_rank + 1) / total_genes
        else:
            avg_density = None
            avg_percentile = None

        if sense_pcts:
            avg_sense = statistics.mean(sense_pcts)
        else:
            avg_sense = None

        group_stats[group_name] = {
            'n_genes': len(genes),
            'avg_density': avg_density,
            'avg_percentile': avg_percentile,
            'avg_sense_pct': avg_sense,
            'densities': densities,
            'sense_pcts': sense_pcts
        }

        print("-" * 70)
        if avg_density:
            print(f"{'AVERAGE':<15} {'':<12} {avg_rank:>8.0f} {avg_percentile:>7.1f}% {avg_density:>12.0f} {avg_sense:>7.1f}%")

    # Cross-group comparison
    print("\n" + "=" * 80)
    print("CROSS-GROUP SUMMARY")
    print("=" * 80)
    print(f"\n{'Group':<15} {'N':>5} {'Avg %ile':>10} {'Avg Density':>14} {'Avg Sense%':>12} {'Strand Bias':>12}")
    print("-" * 70)

    for group_name, stats in group_stats.items():
        n = stats['n_genes']
        pct = f"{stats['avg_percentile']:.1f}%" if stats['avg_percentile'] else "—"
        dens = f"{stats['avg_density']:.0f}" if stats['avg_density'] else "—"
        sense = f"{stats['avg_sense_pct']:.1f}%" if stats['avg_sense_pct'] else "—"

        # Interpret strand bias
        if stats['avg_sense_pct']:
            if stats['avg_sense_pct'] > 65:
                bias = "Strong sense"
            elif stats['avg_sense_pct'] > 55:
                bias = "Sense"
            elif stats['avg_sense_pct'] < 45:
                bias = "Antisense"
            elif stats['avg_sense_pct'] < 35:
                bias = "Strong anti"
            else:
                bias = "Balanced"
        else:
            bias = "—"

        print(f"{group_name:<15} {n:>5} {pct:>10} {dens:>14} {sense:>12} {bias:>12}")

    # Genome-wide baseline
    print("-" * 70)
    print(f"{'GENOME-WIDE':<15} {total_genes:>5} {'50.0%':>10} {'—':>14} {'60.2%':>12} {'Sense':>12}")

    print("\n" + "=" * 80)
    print("KEY OBSERVATIONS")
    print("=" * 80)


if __name__ == '__main__':
    main()
