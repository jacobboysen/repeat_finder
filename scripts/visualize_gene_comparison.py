#!/usr/bin/env python3
"""
Visualize TE signal comparison across germ plasm genes.

Generates publication-quality visualizations including:
1. Summary table of TE signal by gene
2. Ranking bar chart of TE signal
3. TE family breakdown by gene
4. Length-normalized comparison
"""

import argparse
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_density_data(density_file):
    """Load density data from pickle file."""
    with open(density_file, 'rb') as f:
        return pickle.load(f)


def extract_gene_name(qseqid):
    """Extract gene symbol from query sequence ID."""
    # Format: gene_FBtr...
    parts = qseqid.split('_')
    if len(parts) >= 2:
        return parts[0]
    return qseqid


def aggregate_by_gene(density_dict):
    """Aggregate statistics by gene (combining isoforms)."""
    gene_stats = defaultdict(lambda: {
        'isoforms': [],
        'total_length': 0,
        'total_hits': 0,
        'max_signal': 0,
        'total_signal': 0,
        'te_families': defaultdict(int)
    })

    for qseqid, data in density_dict.items():
        gene = extract_gene_name(qseqid)

        gene_stats[gene]['isoforms'].append(qseqid)
        gene_stats[gene]['total_length'] += data['length']
        gene_stats[gene]['total_hits'] += data['total_hits']
        gene_stats[gene]['max_signal'] = max(
            gene_stats[gene]['max_signal'],
            np.max(data['combined']) if len(data['combined']) > 0 else 0
        )
        gene_stats[gene]['total_signal'] += np.sum(data['combined'])

        # Count TE families from hits
        for hit in data.get('hits', []):
            te_id = hit['sseqid']
            # Extract family
            if ':' in te_id:
                te_family = te_id.split(':')[0]
            elif '_' in te_id:
                te_family = te_id.split('_')[0]
            else:
                te_family = te_id
            gene_stats[gene]['te_families'][te_family] += 1

    return dict(gene_stats)


def create_summary_table(gene_stats):
    """Create summary DataFrame."""
    rows = []
    for gene, stats in gene_stats.items():
        # Get top TE family
        te_fams = stats['te_families']
        top_te = max(te_fams.keys(), key=lambda x: te_fams[x]) if te_fams else None

        rows.append({
            'gene': gene,
            'num_isoforms': len(stats['isoforms']),
            'total_length': stats['total_length'],
            'total_hits': stats['total_hits'],
            'max_signal': stats['max_signal'],
            'total_signal': stats['total_signal'],
            'signal_per_kb': stats['total_signal'] / (stats['total_length'] / 1000) if stats['total_length'] > 0 else 0,
            'hits_per_kb': stats['total_hits'] / (stats['total_length'] / 1000) if stats['total_length'] > 0 else 0,
            'top_te_family': top_te,
            'num_te_families': len(te_fams)
        })

    return pd.DataFrame(rows)


def plot_signal_ranking(df, output_dir):
    """Bar chart ranking genes by TE signal."""
    df_sorted = df.sort_values('max_signal', ascending=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Absolute signal
    ax1 = axes[0]
    colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(df_sorted)))
    bars = ax1.barh(df_sorted['gene'], df_sorted['max_signal'], color=colors)
    ax1.set_xlabel('Maximum TE Signal (Combined Score)', fontsize=11)
    ax1.set_ylabel('Gene', fontsize=11)
    ax1.set_title('TE Signal by Gene (Absolute)', fontsize=12)

    # Use log scale if needed
    if df_sorted['max_signal'].max() / (df_sorted['max_signal'].min() + 1e-10) > 100:
        ax1.set_xscale('log')

    # Plot 2: Length-normalized
    ax2 = axes[1]
    df_sorted_norm = df.sort_values('signal_per_kb', ascending=True)
    colors = plt.cm.RdYlBu_r(np.linspace(0.2, 0.8, len(df_sorted_norm)))
    ax2.barh(df_sorted_norm['gene'], df_sorted_norm['signal_per_kb'], color=colors)
    ax2.set_xlabel('TE Signal per kb', fontsize=11)
    ax2.set_ylabel('Gene', fontsize=11)
    ax2.set_title('TE Signal by Gene (Length-Normalized)', fontsize=12)

    if df_sorted_norm['signal_per_kb'].max() / (df_sorted_norm['signal_per_kb'].min() + 1e-10) > 100:
        ax2.set_xscale('log')

    plt.tight_layout()

    out_file = output_dir / 'gene_ranking.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_te_family_breakdown(gene_stats, output_dir, top_n_families=8):
    """Stacked bar chart of TE families by gene."""
    # Collect all TE families
    all_families = defaultdict(int)
    for gene, stats in gene_stats.items():
        for family, count in stats['te_families'].items():
            all_families[family] += count

    # Get top families
    top_families = sorted(all_families.keys(), key=lambda x: -all_families[x])[:top_n_families]

    # Build data for stacked bar
    genes = list(gene_stats.keys())
    data = {family: [] for family in top_families}
    data['Other'] = []

    for gene in genes:
        fams = gene_stats[gene]['te_families']
        for family in top_families:
            data[family].append(fams.get(family, 0))
        # Sum others
        other = sum(count for fam, count in fams.items() if fam not in top_families)
        data['Other'].append(other)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(12, 6))

    bottom = np.zeros(len(genes))
    colors = plt.cm.tab20(np.linspace(0, 1, len(top_families) + 1))

    for i, family in enumerate(top_families + ['Other']):
        ax.bar(genes, data[family], bottom=bottom, label=family[:20], color=colors[i])
        bottom += np.array(data[family])

    ax.set_xlabel('Gene', fontsize=11)
    ax.set_ylabel('Number of BLAST Hits', fontsize=11)
    ax.set_title('TE Family Breakdown by Gene', fontsize=12)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    out_file = output_dir / 'te_family_breakdown.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_length_vs_signal(df, output_dir):
    """Scatter plot of 3'UTR length vs TE signal."""
    fig, ax = plt.subplots(figsize=(10, 6))

    scatter = ax.scatter(
        df['total_length'],
        df['total_signal'],
        s=df['total_hits'] / 10 + 50,  # Size by hits
        c=df['num_te_families'],
        cmap='viridis',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )

    # Add labels
    for _, row in df.iterrows():
        ax.annotate(row['gene'], (row['total_length'], row['total_signal']),
                   fontsize=9, alpha=0.8)

    plt.colorbar(scatter, label='Number of TE Families')
    ax.set_xlabel('Total 3\'UTR Length (bp)', fontsize=11)
    ax.set_ylabel('Total TE Signal', fontsize=11)
    ax.set_title('3\'UTR Length vs TE Signal\n(bubble size = hit count)', fontsize=12)

    if df['total_signal'].max() / (df['total_signal'].min() + 1e-10) > 100:
        ax.set_yscale('log')

    plt.tight_layout()

    out_file = output_dir / 'length_vs_signal.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'density_dir',
        type=Path,
        nargs='?',
        help='Directory containing density_data.pkl'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('figures/gene_comparison'),
        help='Output directory for figures'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find density data
    if args.density_dir is None:
        sweep_dir = Path('results/parameter_sweep')
        if sweep_dir.exists():
            subdirs = sorted([d for d in sweep_dir.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                combo_dirs = sorted([d for d in subdirs[0].iterdir()
                                    if d.is_dir() and d.name.startswith('combo_')])
                if combo_dirs:
                    density_dir = combo_dirs[0] / 'density'
                    if density_dir.exists():
                        args.density_dir = density_dir

        if args.density_dir is None:
            print("Error: No density directory found", file=sys.stderr)
            return 1

    density_file = args.density_dir / 'density_data.pkl'
    if not density_file.exists():
        print(f"Error: Density file not found: {density_file}", file=sys.stderr)
        return 1

    print("Gene Comparison Analysis")
    print("=" * 60)
    print(f"Density data: {args.density_dir}")
    print(f"Output: {args.output_dir}")
    print()

    # Load data
    print("Loading density data...")
    density_dict = load_density_data(density_file)
    print(f"  Loaded data for {len(density_dict)} sequences")

    # Aggregate by gene
    print("\nAggregating by gene...")
    gene_stats = aggregate_by_gene(density_dict)
    print(f"  Found {len(gene_stats)} unique genes")

    # Create summary table
    df = create_summary_table(gene_stats)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Save summary table
    summary_file = args.output_dir / 'gene_summary.tsv'
    df.to_csv(summary_file, sep='\t', index=False)
    print(f"\nSaved summary table: {summary_file.name}")

    # Generate plots
    print("\nGenerating plots...")

    out_file = plot_signal_ranking(df, args.output_dir)
    print(f"  Saved: {out_file.name}")

    out_file = plot_te_family_breakdown(gene_stats, args.output_dir)
    print(f"  Saved: {out_file.name}")

    out_file = plot_length_vs_signal(df, args.output_dir)
    print(f"  Saved: {out_file.name}")

    # Print summary
    print("\n" + "=" * 60)
    print("Summary Table:")
    print("-" * 60)
    print(df.sort_values('max_signal', ascending=False).to_string(index=False))

    print("\n" + "=" * 60)
    print("Gene comparison complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
