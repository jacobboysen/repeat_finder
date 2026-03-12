#!/usr/bin/env python3
"""
Visualization suite for TE motif analysis results.

Creates:
1. Volcano plot: enrichment vs -log10(q-value)
2. Top 20 enriched/depleted motifs bar chart
3. Positional heatmap: motifs x deciles
4. Gene set enrichment heatmap
5. Isoform concordance histogram

Output: figures/motif_analysis/
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
from utils.paths import get_project_root, get_results_dir, get_figures_dir


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
        'savefig.facecolor': 'white',
    })


def plot_volcano(df: pd.DataFrame, output_path: Path, sig_threshold: float = 0.05):
    """
    Create volcano plot: enrichment vs significance.
    """
    plot_df = df.copy()

    # Compute transformed values
    plot_df['log2_enrich'] = np.log2(plot_df['enrichment'].clip(lower=0.01, upper=100))
    plot_df['-log10_q'] = -np.log10(plot_df['q_value'].clip(lower=1e-300))

    # Cap extreme values for visualization
    plot_df['-log10_q'] = plot_df['-log10_q'].clip(upper=300)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Color by significance and direction
    colors = []
    sizes = []
    for _, row in plot_df.iterrows():
        if row['q_value'] < sig_threshold:
            if row['log2_enrich'] > 1:  # >2x enriched
                colors.append('#d62728')  # Red
                sizes.append(60)
            elif row['log2_enrich'] < -1:  # <0.5x depleted
                colors.append('#1f77b4')  # Blue
                sizes.append(60)
            else:
                colors.append('#ff7f0e')  # Orange (significant but small effect)
                sizes.append(40)
        else:
            colors.append('#cccccc')  # Gray
            sizes.append(20)

    ax.scatter(
        plot_df['log2_enrich'],
        plot_df['-log10_q'],
        c=colors,
        s=sizes,
        alpha=0.6,
        edgecolors='white',
        linewidth=0.3
    )

    # Significance threshold line
    sig_line = -np.log10(sig_threshold)
    ax.axhline(y=sig_line, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=0, color='gray', linestyle='-', alpha=0.3, linewidth=1)
    ax.axvline(x=1, color='gray', linestyle=':', alpha=0.3, linewidth=1)  # 2x enriched
    ax.axvline(x=-1, color='gray', linestyle=':', alpha=0.3, linewidth=1)  # 0.5x

    # Label top significant motifs
    sig_df = plot_df[(plot_df['q_value'] < sig_threshold) &
                     ((plot_df['log2_enrich'] > 1) | (plot_df['log2_enrich'] < -1))]
    sig_df = sig_df.copy()
    sig_df['importance'] = sig_df['-log10_q'] * abs(sig_df['log2_enrich'])

    top_labels = sig_df.nlargest(25, 'importance')
    for _, row in top_labels.iterrows():
        ax.annotate(
            row['motif'],
            xy=(row['log2_enrich'], row['-log10_q']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=7,
            alpha=0.9,
            fontweight='bold' if row['-log10_q'] > 50 else 'normal'
        )

    ax.set_xlabel('log2(Enrichment)', fontsize=12)
    ax.set_ylabel('-log10(q-value)', fontsize=12)
    ax.set_title('TE Motif Enrichment: Real vs Shuffled Controls', fontsize=14, fontweight='bold')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', label='Enriched >2x (q < 0.05)'),
        Patch(facecolor='#1f77b4', label='Depleted <0.5x (q < 0.05)'),
        Patch(facecolor='#ff7f0e', label='Significant (small effect)'),
        Patch(facecolor='#cccccc', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_top_motifs(df: pd.DataFrame, output_path: Path, n_top: int = 20):
    """
    Create bar chart of top enriched and depleted motifs.
    """
    # Filter significant
    sig_df = df[df['q_value'] < 0.05].copy()

    if len(sig_df) == 0:
        print("  Warning: No significant motifs for bar chart")
        sig_df = df.copy()

    sig_df['log2_enrich'] = np.log2(sig_df['enrichment'].clip(lower=0.01, upper=100))

    # Get top enriched and depleted
    top_enriched = sig_df.nlargest(n_top, 'enrichment')
    top_depleted = sig_df.nsmallest(n_top, 'enrichment')

    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    # Enriched motifs
    ax = axes[0]
    if len(top_enriched) > 0:
        colors = plt.cm.Reds(np.linspace(0.4, 0.9, len(top_enriched)))
        bars = ax.barh(range(len(top_enriched)), top_enriched['log2_enrich'],
                       color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(top_enriched)))
        ax.set_yticklabels(top_enriched['motif'], fontfamily='monospace', fontsize=10)
        ax.set_xlabel('log2(Enrichment)', fontsize=11)
        ax.set_title(f'Top {n_top} Enriched Motifs', fontsize=12, fontweight='bold')
        ax.invert_yaxis()

        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, top_enriched['enrichment'])):
            ax.text(bar.get_width() + 0.1, i, f'{val:.1f}x', va='center', fontsize=8)

    # Depleted motifs
    ax = axes[1]
    if len(top_depleted) > 0:
        colors = plt.cm.Blues(np.linspace(0.4, 0.9, len(top_depleted)))
        bars = ax.barh(range(len(top_depleted)), top_depleted['log2_enrich'],
                       color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(top_depleted)))
        ax.set_yticklabels(top_depleted['motif'], fontfamily='monospace', fontsize=10)
        ax.set_xlabel('log2(Enrichment)', fontsize=11)
        ax.set_title(f'Top {n_top} Depleted Motifs', fontsize=12, fontweight='bold')

        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, top_depleted['enrichment'])):
            ax.text(bar.get_width() - 0.1, i, f'{val:.2f}x', va='center', ha='right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_position_heatmap(df: pd.DataFrame, output_path: Path, n_motifs: int = 30):
    """
    Create positional heatmap: motifs x position bins.
    """
    # Select top motifs by total count
    motif_totals = df.groupby('motif')['real_count'].sum().sort_values(ascending=False)
    top_motifs = motif_totals.head(n_motifs).index.tolist()

    # Filter to top motifs
    plot_df = df[df['motif'].isin(top_motifs)]

    # Pivot for heatmap
    pivot = plot_df.pivot(index='motif', columns='bin', values='real_density')

    # Reorder by 5' vs 3' bias
    five_prime = pivot[[1, 2, 3]].sum(axis=1)
    three_prime = pivot[[8, 9, 10]].sum(axis=1)
    bias = five_prime - three_prime
    pivot = pivot.loc[bias.sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(12, max(8, len(top_motifs) * 0.3)))

    sns.heatmap(
        pivot,
        ax=ax,
        cmap='YlOrRd',
        cbar_kws={'label': 'Density (fraction of motif occurrences)'},
        linewidths=0.5,
        linecolor='white'
    )

    ax.set_xlabel('Position Decile (1=5\' end, 10=3\' end)', fontsize=11)
    ax.set_ylabel('Motif', fontsize=11)
    ax.set_title('Motif Positional Distribution within 3\'UTRs', fontsize=12, fontweight='bold')

    # Y-axis labels
    ax.set_yticklabels(ax.get_yticklabels(), fontfamily='monospace', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_geneset_heatmap(df: pd.DataFrame, output_path: Path,
                         n_motifs: int = 20, n_genesets: int = 20):
    """
    Create gene set enrichment heatmap.
    """
    # Filter significant results
    sig_df = df[df['fisher_q'] < 0.1].copy()

    if len(sig_df) < 5:
        print("  Warning: Too few significant gene set associations for heatmap")
        sig_df = df.nsmallest(100, 'fisher_p').copy()

    # Get top motifs and gene sets by significance
    top_motifs = sig_df.groupby('motif')['fisher_p'].min().nsmallest(n_motifs).index.tolist()
    top_genesets = sig_df.groupby('gene_set')['fisher_p'].min().nsmallest(n_genesets).index.tolist()

    # Filter
    plot_df = sig_df[sig_df['motif'].isin(top_motifs) & sig_df['gene_set'].isin(top_genesets)]

    if len(plot_df) == 0:
        print("  Warning: No data for gene set heatmap")
        return

    # Pivot for heatmap - use log2 odds ratio
    plot_df['log2_or'] = np.log2(plot_df['fisher_or'].clip(lower=0.1, upper=10))
    pivot = plot_df.pivot(index='motif', columns='gene_set', values='log2_or')

    # Clean gene set names for display
    pivot.columns = [c.replace('expr_', '').replace('go_', 'GO:').replace('flyfish_', 'FF:')[:25]
                     for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(14, max(8, len(top_motifs) * 0.35)))

    sns.heatmap(
        pivot.fillna(0),
        ax=ax,
        cmap='coolwarm',
        center=0,
        vmin=-2,
        vmax=2,
        cbar_kws={'label': 'log2(Odds Ratio)'},
        linewidths=0.5,
        linecolor='white'
    )

    ax.set_xlabel('Gene Set', fontsize=11)
    ax.set_ylabel('Motif', fontsize=11)
    ax.set_title('Motif-Gene Set Associations', fontsize=12, fontweight='bold')

    ax.set_yticklabels(ax.get_yticklabels(), fontfamily='monospace', fontsize=9)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_isoform_concordance(df: pd.DataFrame, output_path: Path):
    """
    Create isoform concordance histogram.
    """
    if len(df) == 0:
        print("  Warning: No isoform data for histogram")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel 1: Distribution of isoform specificity
    ax = axes[0]
    specificity_counts = df['specificity'].value_counts()
    colors = {'unique_to_one': '#d62728', 'minority': '#ff7f0e', 'majority': '#2ca02c'}
    ax.bar(specificity_counts.index, specificity_counts.values,
           color=[colors.get(x, 'gray') for x in specificity_counts.index],
           edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Specificity Type', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('A. Isoform Specificity Distribution', fontsize=12, fontweight='bold')

    # Panel 2: Number of specific motifs per gene
    ax = axes[1]
    gene_counts = df.groupby('fbgn')['motif'].nunique()
    ax.hist(gene_counts, bins=20, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Number of Isoform-Specific Motifs', fontsize=11)
    ax.set_ylabel('Number of Genes', fontsize=11)
    ax.set_title('B. Motifs per Gene', fontsize=12, fontweight='bold')
    ax.axvline(gene_counts.median(), color='red', linestyle='--', label=f'Median: {gene_counts.median():.1f}')
    ax.legend()

    # Panel 3: Number of genes per motif
    ax = axes[2]
    motif_counts = df.groupby('motif')['fbgn'].nunique().sort_values(ascending=False)
    top_motifs = motif_counts.head(15)
    ax.barh(range(len(top_motifs)), top_motifs.values, color='coral', edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(top_motifs)))
    ax.set_yticklabels(top_motifs.index, fontfamily='monospace', fontsize=9)
    ax.set_xlabel('Number of Genes', fontsize=11)
    ax.set_title('C. Top Motifs by Gene Count', fontsize=12, fontweight='bold')
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_summary_dashboard(
    motif_df: pd.DataFrame,
    position_df: pd.DataFrame,
    geneset_df: pd.DataFrame,
    isoform_df: pd.DataFrame,
    output_path: Path
):
    """
    Create summary dashboard combining all analyses.
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Panel A: Mini volcano plot
    ax = axes[0, 0]
    plot_df = motif_df.copy()
    plot_df['log2_enrich'] = np.log2(plot_df['enrichment'].clip(lower=0.01, upper=100))
    plot_df['-log10_q'] = -np.log10(plot_df['q_value'].clip(lower=1e-300)).clip(upper=200)

    colors = ['#d62728' if (q < 0.05 and e > 1) else '#1f77b4' if (q < 0.05 and e < -1) else '#cccccc'
              for q, e in zip(plot_df['q_value'], plot_df['log2_enrich'])]
    ax.scatter(plot_df['log2_enrich'], plot_df['-log10_q'], c=colors, alpha=0.5, s=10)
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
    ax.set_xlabel('log2(Enrichment)')
    ax.set_ylabel('-log10(q-value)')
    ax.set_title('A. Motif Enrichment Volcano')

    # Panel B: Top motifs bar
    ax = axes[0, 1]
    sig_motifs = motif_df[motif_df['q_value'] < 0.05].nlargest(10, 'enrichment')
    if len(sig_motifs) > 0:
        ax.barh(range(len(sig_motifs)), np.log2(sig_motifs['enrichment']),
                color='#d62728', edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(sig_motifs)))
        ax.set_yticklabels(sig_motifs['motif'], fontfamily='monospace', fontsize=9)
        ax.set_xlabel('log2(Enrichment)')
        ax.invert_yaxis()
    ax.set_title('B. Top Enriched Motifs')

    # Panel C: Position bias summary
    ax = axes[0, 2]
    if len(position_df) > 0:
        # Average position by motif enrichment
        motif_positions = position_df.groupby('motif').apply(
            lambda x: (x['bin'] * x['real_density']).sum() / x['real_density'].sum()
            if x['real_density'].sum() > 0 else 5.5
        )
        motif_enrich = motif_df.set_index('motif')['enrichment']
        common = motif_positions.index.intersection(motif_enrich.index)

        if len(common) > 0:
            ax.scatter(np.log2(motif_enrich[common].clip(0.1, 100)),
                       motif_positions[common], alpha=0.5, s=20, c='steelblue')
            ax.axhline(5.5, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel('log2(Enrichment)')
            ax.set_ylabel('Mean Position (1=5\', 10=3\')')
    ax.set_title('C. Enrichment vs Position')

    # Panel D: Gene set summary
    ax = axes[1, 0]
    if len(geneset_df) > 0:
        sig_assoc = geneset_df[geneset_df['fisher_q'] < 0.05]
        motif_assoc = sig_assoc.groupby('motif').size().sort_values(ascending=False).head(10)
        if len(motif_assoc) > 0:
            ax.barh(range(len(motif_assoc)), motif_assoc.values,
                    color='coral', edgecolor='black', linewidth=0.5)
            ax.set_yticks(range(len(motif_assoc)))
            ax.set_yticklabels(motif_assoc.index, fontfamily='monospace', fontsize=9)
            ax.set_xlabel('N Gene Set Associations')
            ax.invert_yaxis()
    ax.set_title('D. Motifs with Gene Set Associations')

    # Panel E: Isoform specificity
    ax = axes[1, 1]
    if len(isoform_df) > 0:
        specificity_counts = isoform_df['specificity'].value_counts()
        colors = {'unique_to_one': '#d62728', 'minority': '#ff7f0e', 'majority': '#2ca02c'}
        ax.bar(specificity_counts.index, specificity_counts.values,
               color=[colors.get(x, 'gray') for x in specificity_counts.index],
               edgecolor='black', linewidth=0.5)
        ax.set_xlabel('Specificity Type')
        ax.set_ylabel('Count')
    ax.set_title('E. Isoform Specificity')

    # Panel F: Summary statistics
    ax = axes[1, 2]
    ax.axis('off')

    n_enriched = len(motif_df[(motif_df['q_value'] < 0.05) & (motif_df['enrichment'] > 2)])
    n_depleted = len(motif_df[(motif_df['q_value'] < 0.05) & (motif_df['enrichment'] < 0.5)])
    n_geneset_assoc = len(geneset_df[geneset_df['fisher_q'] < 0.05]) if len(geneset_df) > 0 else 0
    n_isoform_specific = len(isoform_df[isoform_df['specificity'] == 'unique_to_one']) if len(isoform_df) > 0 else 0

    summary_text = f"""
MOTIF ANALYSIS SUMMARY

K-mer Analysis:
  Total motifs analyzed: {len(motif_df):,}
  Significantly enriched (>2x): {n_enriched}
  Significantly depleted (<0.5x): {n_depleted}

Gene Set Associations:
  Total significant (q < 0.05): {n_geneset_assoc}
  Unique motif-geneset pairs

Isoform Specificity:
  Unique-to-one-isoform: {n_isoform_specific}
  Total genes with specific motifs: {isoform_df['fbgn'].nunique() if len(isoform_df) > 0 else 0}
"""

    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
            fontsize=11, fontfamily='monospace', verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    ax.set_title('F. Summary Statistics')

    plt.suptitle('TE Motif Analysis Dashboard', fontsize=16, fontweight='bold', y=1.02)
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
        '--motif-dir',
        type=Path,
        default=None,
        help='Motif analysis results directory (default: results/motif_analysis)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Output directory for figures (default: figures/motif_analysis)'
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
    figures_dir = get_figures_dir()

    # Set paths
    motif_dir = args.motif_dir or results_dir / "motif_analysis"
    output_dir = args.output_dir or figures_dir / "motif_analysis"

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("MOTIF ANALYSIS VISUALIZATION")
    print("=" * 70)

    # Set up style
    setup_style()

    # Load data files
    print(f"\nLoading data from: {motif_dir}")

    # Motif enrichment results
    motif_path = motif_dir / "motif_enrichment_6mer.tsv"
    if motif_path.exists():
        motif_df = pd.read_csv(motif_path, sep='\t')
        print(f"  Motif enrichment: {len(motif_df)} motifs")
    else:
        print(f"  Warning: {motif_path} not found")
        motif_df = pd.DataFrame()

    # Position results
    position_path = motif_dir / "motif_position_density.tsv"
    if position_path.exists():
        position_df = pd.read_csv(position_path, sep='\t')
        print(f"  Position data: {len(position_df)} rows")
    else:
        print(f"  Warning: {position_path} not found")
        position_df = pd.DataFrame()

    # Gene set results
    geneset_path = motif_dir / "motif_gene_set_enrichment.tsv"
    if geneset_path.exists():
        geneset_df = pd.read_csv(geneset_path, sep='\t')
        print(f"  Gene set associations: {len(geneset_df)} rows")
    else:
        print(f"  Warning: {geneset_path} not found")
        geneset_df = pd.DataFrame()

    # Isoform results
    isoform_path = motif_dir / "isoform_specific_te_motifs.tsv"
    if isoform_path.exists():
        isoform_df = pd.read_csv(isoform_path, sep='\t')
        print(f"  Isoform specificity: {len(isoform_df)} rows")
    else:
        print(f"  Warning: {isoform_path} not found")
        isoform_df = pd.DataFrame()

    # Generate figures
    print(f"\nGenerating figures...")

    # 1. Volcano plot
    if len(motif_df) > 0:
        plot_volcano(motif_df, output_dir / f"motif_volcano.{args.format}")

    # 2. Top motifs bar chart
    if len(motif_df) > 0:
        plot_top_motifs(motif_df, output_dir / f"top_motifs_barplot.{args.format}")

    # 3. Position heatmap
    if len(position_df) > 0:
        plot_position_heatmap(position_df, output_dir / f"motif_position_heatmap.{args.format}")

    # 4. Gene set heatmap
    if len(geneset_df) > 0:
        plot_geneset_heatmap(geneset_df, output_dir / f"motif_geneset_heatmap.{args.format}")

    # 5. Isoform concordance
    if len(isoform_df) > 0:
        plot_isoform_concordance(isoform_df, output_dir / f"isoform_concordance_hist.{args.format}")

    # 6. Summary dashboard
    if len(motif_df) > 0:
        plot_summary_dashboard(
            motif_df, position_df, geneset_df, isoform_df,
            output_dir / f"motif_analysis_dashboard.{args.format}"
        )

    print(f"\nAll figures saved to: {output_dir}")

    print("\n" + "=" * 70)
    print("VISUALIZATION COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
