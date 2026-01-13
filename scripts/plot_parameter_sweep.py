#!/usr/bin/env python3
"""
Visualize BLAST parameter sweep results.

Creates plots showing:
1. Heatmaps of hit counts across parameter combinations
2. Sensitivity vs specificity comparison (sense vs antisense)
3. E-value distribution comparisons
4. Parameter ranking by hit count
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_sweep_results(sweep_dir):
    """Load sweep summary from directory."""
    summary_file = sweep_dir / 'sweep_summary.tsv'
    if not summary_file.exists():
        raise FileNotFoundError(f"Summary not found: {summary_file}")
    return pd.read_csv(summary_file, sep='\t')


def plot_parameter_heatmap(df, output_dir):
    """Create heatmap of hits by word_size and gap parameters."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Filter successful runs
    df_ok = df[df['success'] == True].copy()

    if len(df_ok) == 0:
        print("No successful runs to plot")
        return None

    # Heatmap 1: word_size vs gapopen (aggregate by task)
    ax1 = axes[0]
    pivot1 = df_ok.pivot_table(
        values='total_hits',
        index='gapopen',
        columns='word_size',
        aggfunc='mean'
    )
    sns.heatmap(pivot1, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax1)
    ax1.set_title('Mean Hits by Word Size and Gap Open', fontsize=12)
    ax1.set_xlabel('Word Size')
    ax1.set_ylabel('Gap Open Penalty')

    # Heatmap 2: gapopen vs gapextend
    ax2 = axes[1]
    pivot2 = df_ok.pivot_table(
        values='total_hits',
        index='gapextend',
        columns='gapopen',
        aggfunc='mean'
    )
    sns.heatmap(pivot2, annot=True, fmt='.0f', cmap='YlOrRd', ax=ax2)
    ax2.set_title('Mean Hits by Gap Open and Gap Extend', fontsize=12)
    ax2.set_xlabel('Gap Open Penalty')
    ax2.set_ylabel('Gap Extend Penalty')

    plt.tight_layout()

    out_file = output_dir / 'parameter_heatmap.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_task_comparison(df, output_dir):
    """Compare blastn vs blastn-short."""
    df_ok = df[df['success'] == True].copy()

    if len(df_ok) == 0:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Box plot of hits by task
    ax1 = axes[0]
    sns.boxplot(data=df_ok, x='task', y='total_hits', ax=ax1)
    ax1.set_title('Total Hits by BLAST Task', fontsize=12)
    ax1.set_ylabel('Total Hits')
    ax1.set_xlabel('Task')

    # Box plot of hits by word_size, colored by task
    ax2 = axes[1]
    sns.boxplot(data=df_ok, x='word_size', y='total_hits', hue='task', ax=ax2)
    ax2.set_title('Hits by Word Size and Task', fontsize=12)
    ax2.set_ylabel('Total Hits')
    ax2.set_xlabel('Word Size')
    ax2.legend(title='Task')

    plt.tight_layout()

    out_file = output_dir / 'task_comparison.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_reward_penalty(df, output_dir):
    """Plot effect of reward/penalty combinations."""
    df_ok = df[df['success'] == True].copy()

    if len(df_ok) == 0:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    # Create reward_penalty string for grouping
    df_ok['reward_penalty'] = df_ok['reward'].astype(str) + '/' + df_ok['penalty'].astype(str)

    sns.boxplot(data=df_ok, x='reward_penalty', y='total_hits', ax=ax)
    ax.set_title('Total Hits by Reward/Penalty Combination', fontsize=12)
    ax.set_ylabel('Total Hits')
    ax.set_xlabel('Reward/Penalty')

    plt.tight_layout()

    out_file = output_dir / 'reward_penalty_comparison.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_top_combinations(df, output_dir, top_n=20):
    """Bar chart of top parameter combinations by hits."""
    df_ok = df[df['success'] == True].copy()

    if len(df_ok) == 0:
        return None

    # Sort by hits
    df_sorted = df_ok.nlargest(top_n, 'total_hits')

    # Create label
    df_sorted['label'] = (
        'ws=' + df_sorted['word_size'].astype(str) +
        ' go=' + df_sorted['gapopen'].astype(str) +
        ' r/p=' + df_sorted['reward'].astype(str) + '/' + df_sorted['penalty'].astype(str)
    )

    fig, ax = plt.subplots(figsize=(12, 8))

    colors = ['#377eb8' if t == 'blastn' else '#ff7f00' for t in df_sorted['task']]
    bars = ax.barh(df_sorted['label'], df_sorted['total_hits'], color=colors)

    ax.set_xlabel('Total Hits', fontsize=12)
    ax.set_ylabel('Parameter Combination', fontsize=12)
    ax.set_title(f'Top {top_n} Parameter Combinations by Hit Count', fontsize=12)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#377eb8', label='blastn'),
                       Patch(facecolor='#ff7f00', label='blastn-short')]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()

    out_file = output_dir / 'top_combinations.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def plot_evalue_distribution(df, output_dir):
    """Plot e-value statistics across parameter combinations."""
    df_ok = df[df['success'] == True].copy()

    if len(df_ok) == 0 or 'hits_lt_0.01' not in df_ok.columns:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Distribution of strong hits (e-value < 0.01)
    ax1 = axes[0]
    sns.boxplot(data=df_ok, x='word_size', y='hits_lt_0.01', hue='task', ax=ax1)
    ax1.set_title('Strong Hits (E-value < 0.01)', fontsize=12)
    ax1.set_ylabel('Number of Hits')
    ax1.set_xlabel('Word Size')

    # Plot 2: Scatter of total hits vs strong hits
    ax2 = axes[1]
    scatter = ax2.scatter(df_ok['total_hits'], df_ok['hits_lt_0.01'],
                         c=df_ok['word_size'], cmap='viridis', alpha=0.6)
    plt.colorbar(scatter, ax=ax2, label='Word Size')
    ax2.set_xlabel('Total Hits', fontsize=12)
    ax2.set_ylabel('Strong Hits (E < 0.01)', fontsize=12)
    ax2.set_title('Total vs Strong Hits', fontsize=12)

    plt.tight_layout()

    out_file = output_dir / 'evalue_distribution.png'
    plt.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close()

    return out_file


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'sweep_dir',
        type=Path,
        nargs='?',
        help='Parameter sweep results directory'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('figures/parameter_sweep'),
        help='Output directory for figures'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find sweep directory
    if args.sweep_dir is None:
        sweep_base = Path('results/parameter_sweep')
        if sweep_base.exists():
            subdirs = sorted([d for d in sweep_base.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                args.sweep_dir = subdirs[0]

        if args.sweep_dir is None:
            print("Error: No sweep directory found", file=sys.stderr)
            return 1

    print("Parameter Sweep Visualization")
    print("=" * 60)
    print(f"Sweep directory: {args.sweep_dir}")
    print(f"Output: {args.output_dir}")
    print()

    # Load results
    print("Loading sweep results...")
    try:
        df = load_sweep_results(args.sweep_dir)
        print(f"  Loaded {len(df)} parameter combinations")
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    print("\nGenerating plots...")

    out_file = plot_parameter_heatmap(df, args.output_dir)
    if out_file:
        print(f"  Saved: {out_file.name}")

    out_file = plot_task_comparison(df, args.output_dir)
    if out_file:
        print(f"  Saved: {out_file.name}")

    out_file = plot_reward_penalty(df, args.output_dir)
    if out_file:
        print(f"  Saved: {out_file.name}")

    out_file = plot_top_combinations(df, args.output_dir)
    if out_file:
        print(f"  Saved: {out_file.name}")

    out_file = plot_evalue_distribution(df, args.output_dir)
    if out_file:
        print(f"  Saved: {out_file.name}")

    print("\n" + "=" * 60)
    print("Plotting complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
