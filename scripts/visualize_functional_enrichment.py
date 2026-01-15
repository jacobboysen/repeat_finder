#!/usr/bin/env python3
"""
Visualize TE enrichment results from functional gene set analysis.

Creates:
- Heatmap of TE density by functional category
- Forest plot of odds ratios with confidence intervals
- Bar chart of top/bottom enriched categories
- Volcano plot (effect size vs significance)

Output: figures/enrichment/
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir


def setup_style():
    """Set up matplotlib style for publication-quality figures."""
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    })


def plot_density_heatmap(df: pd.DataFrame, output_path: Path):
    """
    Create heatmap of TE density by gene set category.
    """
    # Pivot data for heatmap
    # Select key metrics
    metrics = ['mean_density', 'pct_with_hits', 'sense_pct']

    # Filter to reasonable number of gene sets per category
    plot_data = []
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat].nlargest(10, 'mean_density')
        plot_data.append(cat_df)

    if not plot_data:
        print("  Warning: No data for heatmap")
        return

    subset = pd.concat(plot_data)

    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, max(6, len(subset) * 0.25)))

    # Mean density heatmap
    ax = axes[0]
    pivot = subset.pivot_table(
        values='mean_density',
        index='gene_set',
        columns='category',
        aggfunc='first'
    )
    sns.heatmap(
        pivot.fillna(0),
        ax=ax,
        cmap='YlOrRd',
        cbar_kws={'label': 'Mean density (hits/kb)'}
    )
    ax.set_title('TE Density by Gene Set')
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Percent with hits
    ax = axes[1]
    pivot = subset.pivot_table(
        values='pct_with_hits',
        index='gene_set',
        columns='category',
        aggfunc='first'
    )
    sns.heatmap(
        pivot.fillna(0),
        ax=ax,
        cmap='Blues',
        cbar_kws={'label': '% genes with hits'}
    )
    ax.set_title('TE Hit Prevalence')
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Strand bias
    ax = axes[2]
    pivot = subset.pivot_table(
        values='sense_pct',
        index='gene_set',
        columns='category',
        aggfunc='first'
    )
    sns.heatmap(
        pivot.fillna(50),
        ax=ax,
        cmap='coolwarm',
        center=50,
        vmin=20,
        vmax=80,
        cbar_kws={'label': '% sense strand'}
    )
    ax.set_title('Strand Bias')
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_forest(df: pd.DataFrame, output_path: Path, sig_threshold: float = 0.05):
    """
    Create forest plot of odds ratios with significance indicators.
    """
    # Sort by odds ratio - handle inf values by capping
    plot_df = df.copy()
    # Replace inf with a large but finite value for visualization
    plot_df['fisher_or_capped'] = plot_df['fisher_or'].replace([np.inf, -np.inf], np.nan)
    plot_df.loc[plot_df['fisher_or'] == np.inf, 'fisher_or_capped'] = 100  # Cap at 100x
    plot_df['log2_or'] = np.log2(plot_df['fisher_or_capped'].replace(0, 0.01))  # Cap at 0.01x
    plot_df = plot_df.dropna(subset=['log2_or'])
    plot_df = plot_df.sort_values('log2_or')

    # Limit to most extreme values
    n_show = min(40, len(plot_df))
    top_enriched = plot_df.nlargest(n_show // 2, 'log2_or')
    top_depleted = plot_df.nsmallest(n_show // 2, 'log2_or')
    plot_df = pd.concat([top_depleted, top_enriched])

    if len(plot_df) == 0:
        print("  Warning: No data for forest plot")
        return

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.3)))

    # Colors based on significance
    colors = []
    for _, row in plot_df.iterrows():
        if row['fisher_q'] < sig_threshold:
            if row['log2_or'] > 0:
                colors.append('#d62728')  # Red for enriched
            else:
                colors.append('#1f77b4')  # Blue for depleted
        else:
            colors.append('#7f7f7f')  # Gray for non-significant

    y_pos = range(len(plot_df))

    # Plot horizontal bars
    ax.barh(y_pos, plot_df['log2_or'], color=colors, height=0.7, alpha=0.8)

    # Add significance markers
    for i, (_, row) in enumerate(plot_df.iterrows()):
        if row['fisher_q'] < 0.001:
            marker = '***'
        elif row['fisher_q'] < 0.01:
            marker = '**'
        elif row['fisher_q'] < sig_threshold:
            marker = '*'
        else:
            marker = ''

        x_pos = row['log2_or']
        offset = 0.1 if x_pos >= 0 else -0.1
        ax.text(x_pos + offset, i, marker, va='center', fontsize=8)

    # Reference line at 0
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['gene_set'])
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_title('TE Enrichment by Gene Set (Fisher\'s Exact Test)')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='Enriched (q < 0.05)'),
        Patch(facecolor='#1f77b4', label='Depleted (q < 0.05)'),
        Patch(facecolor='#7f7f7f', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_volcano(df: pd.DataFrame, output_path: Path, sig_threshold: float = 0.05):
    """
    Create volcano plot (effect size vs significance).
    """
    plot_df = df.copy()
    plot_df['log2_or'] = np.log2(plot_df['fisher_or'].replace(0, np.nan).replace(np.inf, 100))
    plot_df['-log10_q'] = -np.log10(plot_df['fisher_q'].replace(0, 1e-200).clip(lower=1e-200))
    plot_df = plot_df.dropna(subset=['log2_or', '-log10_q'])

    # Cap extreme values for better visualization
    plot_df['log2_or'] = plot_df['log2_or'].clip(-10, 10)
    plot_df['-log10_q'] = plot_df['-log10_q'].clip(0, 200)

    if len(plot_df) == 0:
        print("  Warning: No data for volcano plot")
        return

    fig, ax = plt.subplots(figsize=(12, 10))

    # Color by significance and direction
    colors = []
    sizes = []
    for _, row in plot_df.iterrows():
        if row['fisher_q'] < sig_threshold:
            if row['log2_or'] > 0:
                colors.append('#d62728')
            else:
                colors.append('#1f77b4')
            sizes.append(80)
        else:
            colors.append('#cccccc')
            sizes.append(40)

    ax.scatter(
        plot_df['log2_or'],
        plot_df['-log10_q'],
        c=colors,
        s=sizes,
        alpha=0.7,
        edgecolors='white',
        linewidth=0.5
    )

    # Significance threshold line
    sig_line = -np.log10(sig_threshold)
    ax.axhline(y=sig_line, color='gray', linestyle='--', alpha=0.5, label=f'q={sig_threshold}')
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3)

    # Label only top significant points (avoid overlap)
    sig_df = plot_df[plot_df['fisher_q'] < sig_threshold].copy()
    sig_df['importance'] = sig_df['-log10_q'] + abs(sig_df['log2_or'])
    top_labels = sig_df.nlargest(20, 'importance')

    # Try to use adjustText for better label placement
    try:
        from adjustText import adjust_text
        use_adjust = True
    except ImportError:
        use_adjust = False

    texts = []
    for i, (_, row) in enumerate(top_labels.iterrows()):
        # Clean up label
        label = row['gene_set'].replace('expr_', '').replace('go_', 'GO:').replace('flyfish_', 'FF:')
        label = label.replace('_high', ' (high)').replace('_enriched', ' (enr)')
        label = label.replace('_specific', ' (spec)').replace('_', ' ').title()

        # Add jitter based on position to reduce overlap
        x_jitter = 0.1 * (i % 3 - 1)
        y_jitter = 2 * (i % 4 - 1.5)

        txt = ax.text(
            row['log2_or'] + x_jitter,
            row['-log10_q'] + y_jitter,
            label,
            fontsize=8,
            fontweight='bold' if row['-log10_q'] > 100 else 'normal',
            alpha=0.9,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none')
        )
        texts.append(txt)

    # Use adjustText if available
    if use_adjust and texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5))

    ax.set_xlabel('log2(Odds Ratio)', fontsize=12)
    ax.set_ylabel('-log10(q-value)', fontsize=12)
    ax.set_title('TE Enrichment: Effect Size vs Significance', fontsize=14, fontweight='bold')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='Enriched (more TE hits)'),
        Patch(facecolor='#1f77b4', label='Depleted (fewer TE hits)'),
        Patch(facecolor='#cccccc', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    # Add quadrant labels
    ax.text(0.95, 0.95, 'TE Enriched\n(highly significant)',
            transform=ax.transAxes, ha='right', va='top', fontsize=9, color='#d62728', alpha=0.7)
    ax.text(0.05, 0.95, 'TE Depleted\n(highly significant)',
            transform=ax.transAxes, ha='left', va='top', fontsize=9, color='#1f77b4', alpha=0.7)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_category_summary(df: pd.DataFrame, output_path: Path):
    """
    Create bar chart summarizing results by category.
    """
    # Compute category statistics
    summary = []
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        summary.append({
            'category': cat,
            'n_sets': len(cat_df),
            'mean_density': cat_df['mean_density'].mean(),
            'sig_enriched': ((cat_df['fisher_q'] < 0.05) & (cat_df['fisher_or'] > 1)).sum(),
            'sig_depleted': ((cat_df['fisher_q'] < 0.05) & (cat_df['fisher_or'] < 1)).sum(),
        })

    summary_df = pd.DataFrame(summary)

    fig, axes = plt.subplots(1, 3, figsize=(12, 5))

    # Number of gene sets per category
    ax = axes[0]
    ax.bar(summary_df['category'], summary_df['n_sets'], color='steelblue')
    ax.set_xlabel('Category')
    ax.set_ylabel('Number of gene sets')
    ax.set_title('Gene Sets by Category')
    ax.tick_params(axis='x', rotation=45)

    # Mean TE density per category
    ax = axes[1]
    ax.bar(summary_df['category'], summary_df['mean_density'], color='coral')
    ax.set_xlabel('Category')
    ax.set_ylabel('Mean TE density (hits/kb)')
    ax.set_title('TE Density by Category')
    ax.tick_params(axis='x', rotation=45)

    # Significant enrichment/depletion
    ax = axes[2]
    x = range(len(summary_df))
    width = 0.35
    ax.bar([i - width/2 for i in x], summary_df['sig_enriched'],
           width, label='Enriched', color='#d62728')
    ax.bar([i + width/2 for i in x], summary_df['sig_depleted'],
           width, label='Depleted', color='#1f77b4')
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df['category'])
    ax.set_xlabel('Category')
    ax.set_ylabel('Number of significant sets')
    ax.set_title('Significant Enrichment/Depletion')
    ax.legend()
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_top_sets(df: pd.DataFrame, output_path: Path, n_top: int = 15):
    """
    Create bar chart of top enriched/depleted gene sets.
    """
    # Filter significant results
    sig_df = df[df['fisher_q'] < 0.05].copy()

    if len(sig_df) == 0:
        print("  Warning: No significant results for top sets plot")
        sig_df = df.copy()

    # Handle inf values by capping
    sig_df['fisher_or_capped'] = sig_df['fisher_or'].replace([np.inf, -np.inf], np.nan)
    sig_df.loc[sig_df['fisher_or'] == np.inf, 'fisher_or_capped'] = 100
    sig_df['log2_or'] = np.log2(sig_df['fisher_or_capped'].replace(0, 0.01))
    sig_df = sig_df.dropna(subset=['log2_or'])

    # Get top enriched and depleted
    top_enriched = sig_df.nlargest(n_top, 'log2_or')
    top_depleted = sig_df.nsmallest(n_top, 'log2_or')

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    def clean_label(name):
        """Clean gene set name for display."""
        label = name.replace('expr_', '').replace('go_', 'GO: ').replace('flyfish_', 'FlyFISH: ')
        label = label.replace('_high', ' (high expr)')
        label = label.replace('_enriched', ' (tissue-enriched)')
        label = label.replace('_specific', ' (tissue-specific)')
        label = label.replace('_', ' ').title()
        return label

    # Top enriched
    ax = axes[0]
    if len(top_enriched) > 0:
        colors = [plt.cm.Reds(0.4 + 0.4 * i / len(top_enriched))
                  for i in range(len(top_enriched))]
        labels = [clean_label(x) for x in top_enriched['gene_set']]
        bars = ax.barh(range(len(top_enriched)), top_enriched['log2_or'], color=colors)
        ax.set_yticks(range(len(top_enriched)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel('log2(Odds Ratio)', fontsize=11)
        ax.set_title(f'Top {n_top} TE-Enriched Gene Sets', fontsize=12, fontweight='bold')
        ax.invert_yaxis()
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, top_enriched['log2_or'])):
            ax.text(val + 0.1, i, f'{2**val:.1f}x', va='center', fontsize=8)

    # Top depleted
    ax = axes[1]
    if len(top_depleted) > 0:
        colors = [plt.cm.Blues(0.4 + 0.4 * i / len(top_depleted))
                  for i in range(len(top_depleted))]
        labels = [clean_label(x) for x in top_depleted['gene_set']]
        bars = ax.barh(range(len(top_depleted)), top_depleted['log2_or'], color=colors)
        ax.set_yticks(range(len(top_depleted)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel('log2(Odds Ratio)', fontsize=11)
        ax.set_title(f'Top {n_top} TE-Depleted Gene Sets', fontsize=12, fontweight='bold')
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, top_depleted['log2_or'])):
            ax.text(val - 0.3, i, f'{2**val:.2f}x', va='center', fontsize=8, ha='right')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--results',
        type=Path,
        default=None,
        help='Enrichment results file (default: results/functional_te_enrichment.tsv)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Output directory for figures (default: figures/enrichment)'
    )
    parser.add_argument(
        '--sig-threshold',
        type=float,
        default=0.05,
        help='Significance threshold for q-values (default: 0.05)'
    )
    parser.add_argument(
        '--format',
        choices=['png', 'pdf', 'svg'],
        default='png',
        help='Output format (default: png)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()

    # Set paths
    if args.results:
        results_path = args.results
    else:
        results_path = results_dir / "functional_te_enrichment.tsv"

    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = project_root / "figures" / "enrichment"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Check input
    if not results_path.exists():
        print(f"Error: Results file not found: {results_path}")
        print("Run analyze_functional_enrichment.py first.")
        return 1

    # Load results
    print(f"Loading results from: {results_path}")
    df = pd.read_csv(results_path, sep='\t')
    print(f"  {len(df)} gene sets")

    # Set up plotting style
    setup_style()

    # Create figures
    print()
    print("Generating figures...")

    plot_density_heatmap(
        df,
        output_dir / f"te_density_heatmap.{args.format}"
    )

    plot_forest(
        df,
        output_dir / f"te_enrichment_forest.{args.format}",
        sig_threshold=args.sig_threshold
    )

    plot_volcano(
        df,
        output_dir / f"te_enrichment_volcano.{args.format}",
        sig_threshold=args.sig_threshold
    )

    plot_category_summary(
        df,
        output_dir / f"te_enrichment_by_category.{args.format}"
    )

    plot_top_sets(
        df,
        output_dir / f"te_enrichment_top_sets.{args.format}"
    )

    print()
    print(f"All figures saved to: {output_dir}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
