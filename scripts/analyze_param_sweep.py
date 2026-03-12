#!/usr/bin/env python3
"""
Comprehensive analysis of parameter sweep results.

Analyzes:
1. Positional distribution along UTRs
2. Match quality metrics (pident, length, evalue, gaps)
3. TE family enrichment

Usage:
    python scripts/analyze_param_sweep.py --sweep-dir results/param_sweep_full
"""

import argparse
import sys
import json
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.blast_io import BLAST_COLUMNS

# Extended columns for this analysis (17 columns with stitle)
SWEEP_COLUMNS = BLAST_COLUMNS[:16] + ['stitle']

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# UTR length bins
UTR_LEN_BINS = [
    (0, 200, '<200bp'),
    (200, 500, '200-500bp'),
    (500, 1000, '500-1kb'),
    (1000, 2000, '1-2kb'),
    (2000, 100000, '>2kb'),
]

# Quality bins
PIDENT_BINS = [(60, 70), (70, 80), (80, 90), (90, 100)]
LENGTH_BINS = [(15, 30), (30, 50), (50, 100), (100, 200), (200, 1000)]


def load_blast_results(filepath: Path) -> pd.DataFrame:
    """Load BLAST results with derived columns."""
    if not filepath.exists() or filepath.stat().st_size == 0:
        return pd.DataFrame()

    df = pd.read_csv(filepath, sep='\t', header=None, names=SWEEP_COLUMNS)

    # Derived columns
    df['rel_start'] = df[['qstart', 'qend']].min(axis=1) / df['qlen']
    df['rel_end'] = df[['qstart', 'qend']].max(axis=1) / df['qlen']
    df['rel_mid'] = (df['rel_start'] + df['rel_end']) / 2

    # UTR length bin
    def get_len_bin(qlen):
        for lo, hi, label in UTR_LEN_BINS:
            if lo <= qlen < hi:
                return label
        return '>2kb'
    df['utr_len_bin'] = df['qlen'].apply(get_len_bin)

    # TE family from sseqid (e.g., "FBti0019293" or "roo{}" -> extract family)
    def get_te_family(sseqid):
        # Handle different TE ID formats
        if '{}' in str(sseqid):
            return sseqid.split('{}')[0]
        elif '_' in str(sseqid):
            return sseqid.split('_')[0]
        else:
            return sseqid
    df['te_family'] = df['sseqid'].apply(get_te_family)

    # Has gaps
    df['has_gaps'] = df['gapopen'] > 0

    return df


def analyze_positional_distribution(df: pd.DataFrame, n_bins: int = 10) -> Dict:
    """Analyze positional distribution along UTRs."""
    results = {
        'overall': {},
        'by_utr_len': {},
    }

    if len(df) == 0:
        return results

    # Overall distribution
    counts, bin_edges = np.histogram(df['rel_mid'], bins=n_bins, range=(0, 1))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    density = counts / counts.sum() if counts.sum() > 0 else counts

    results['overall'] = {
        'bin_centers': bin_centers.tolist(),
        'counts': counts.tolist(),
        'density': density.tolist(),
        'mean': float(df['rel_mid'].mean()),
        'std': float(df['rel_mid'].std()),
        'five_prime_frac': float((df['rel_mid'] < 0.25).mean()),
        'three_prime_frac': float((df['rel_mid'] >= 0.75).mean()),
    }

    # By UTR length
    for lo, hi, label in UTR_LEN_BINS:
        subset = df[(df['qlen'] >= lo) & (df['qlen'] < hi)]
        if len(subset) < 10:
            continue

        counts, _ = np.histogram(subset['rel_mid'], bins=n_bins, range=(0, 1))
        density = counts / counts.sum() if counts.sum() > 0 else counts

        results['by_utr_len'][label] = {
            'n_hits': len(subset),
            'n_seqs': subset['qseqid'].nunique(),
            'counts': counts.tolist(),
            'density': density.tolist(),
            'mean': float(subset['rel_mid'].mean()),
        }

    return results


def analyze_quality_metrics(df: pd.DataFrame) -> Dict:
    """Analyze match quality metrics."""
    results = {}

    if len(df) == 0:
        return results

    # Basic stats
    for col in ['pident', 'length', 'evalue', 'gapopen', 'mismatch']:
        if col in df.columns:
            results[col] = {
                'mean': float(df[col].mean()),
                'median': float(df[col].median()),
                'std': float(df[col].std()),
                'min': float(df[col].min()),
                'max': float(df[col].max()),
                'q25': float(df[col].quantile(0.25)),
                'q75': float(df[col].quantile(0.75)),
            }

    # Percent with gaps
    results['pct_with_gaps'] = float((df['gapopen'] > 0).mean() * 100)

    # E-value thresholds
    results['pct_e_lt_1'] = float((df['evalue'] < 1).mean() * 100)
    results['pct_e_lt_01'] = float((df['evalue'] < 0.1).mean() * 100)
    results['pct_e_lt_001'] = float((df['evalue'] < 0.01).mean() * 100)

    # Length thresholds
    results['pct_len_ge_50'] = float((df['length'] >= 50).mean() * 100)
    results['pct_len_ge_100'] = float((df['length'] >= 100).mean() * 100)
    results['pct_len_ge_200'] = float((df['length'] >= 200).mean() * 100)

    # Percent identity distribution
    pident_dist = []
    for lo, hi in PIDENT_BINS:
        pct = float(((df['pident'] >= lo) & (df['pident'] < hi)).mean() * 100)
        pident_dist.append({'range': f'{lo}-{hi}%', 'pct': pct})
    results['pident_distribution'] = pident_dist

    # Length distribution
    length_dist = []
    for lo, hi in LENGTH_BINS:
        pct = float(((df['length'] >= lo) & (df['length'] < hi)).mean() * 100)
        length_dist.append({'range': f'{lo}-{hi}bp', 'pct': pct})
    results['length_distribution'] = length_dist

    return results


def analyze_te_families(df: pd.DataFrame, top_n: int = 20) -> Dict:
    """Analyze TE family enrichment."""
    results = {}

    if len(df) == 0:
        return results

    # Count hits by TE family
    family_counts = df['te_family'].value_counts()

    results['n_families'] = len(family_counts)
    results['total_hits'] = len(df)

    # Top families
    top_families = family_counts.head(top_n)
    results['top_families'] = [
        {'family': fam, 'hits': int(count), 'pct': float(count / len(df) * 100)}
        for fam, count in top_families.items()
    ]

    # Family diversity (Shannon entropy)
    probs = family_counts / family_counts.sum()
    entropy = -np.sum(probs * np.log2(probs + 1e-10))
    results['shannon_entropy'] = float(entropy)

    # Concentration (top 5 families)
    top5_pct = float(family_counts.head(5).sum() / len(df) * 100)
    results['top5_concentration'] = top5_pct

    return results


def load_all_sweep_results(sweep_dir: Path) -> pd.DataFrame:
    """Load summary of all sweep results."""
    summary_file = sweep_dir / 'sweep_summary.tsv'
    if summary_file.exists():
        return pd.read_csv(summary_file, sep='\t')
    return pd.DataFrame()


def run_full_analysis(sweep_dir: Path, output_dir: Path):
    """Run comprehensive analysis on all sweep results."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PARAMETER SWEEP ANALYSIS")
    print(f"Input: {sweep_dir}")
    print(f"Output: {output_dir}")
    print("=" * 70)

    # Load sweep summary
    summary_df = load_all_sweep_results(sweep_dir)
    if len(summary_df) == 0:
        print("Error: No sweep summary found")
        return 1

    print(f"\nFound {len(summary_df)} result files")

    # Get unique parameter combinations
    param_combos = summary_df[summary_df['seq_type'] == 'real']['combo'].unique()
    print(f"Parameter combinations: {len(param_combos)}")

    # Analyze each combination
    all_results = []

    for combo in param_combos:
        print(f"\nAnalyzing: {combo}")

        # Load real and shuffled data
        real_file = sweep_dir / f"{combo}_real.tsv"
        shuf_file = sweep_dir / f"{combo}_shuffled.tsv"

        real_df = load_blast_results(real_file)
        shuf_df = load_blast_results(shuf_file)

        if len(real_df) == 0:
            print(f"  -> No real data, skipping")
            continue

        # Extract parameters from combo name
        # Format: ws{ws}_go{go}_ge{ge}_p{p}_dust{dust}
        parts = combo.split('_')
        params = {
            'combo': combo,
            'word_size': int(parts[0].replace('ws', '')),
            'gapopen': int(parts[1].replace('go', '')),
            'gapextend': int(parts[2].replace('ge', '')),
            'penalty': int(parts[3].replace('p', '')),
            'dust': parts[4].replace('dust', ''),
        }

        # Basic counts
        result = {
            **params,
            'n_hits_real': len(real_df),
            'n_hits_shuf': len(shuf_df),
            'n_seqs_real': real_df['qseqid'].nunique() if len(real_df) > 0 else 0,
            'n_seqs_shuf': shuf_df['qseqid'].nunique() if len(shuf_df) > 0 else 0,
            'real_shuf_ratio': len(real_df) / len(shuf_df) if len(shuf_df) > 0 else np.nan,
        }

        # Positional analysis
        pos_real = analyze_positional_distribution(real_df)
        pos_shuf = analyze_positional_distribution(shuf_df)

        if pos_real.get('overall'):
            result['pos_mean_real'] = pos_real['overall']['mean']
            result['pos_5prime_real'] = pos_real['overall']['five_prime_frac']
            result['pos_3prime_real'] = pos_real['overall']['three_prime_frac']

        if pos_shuf.get('overall'):
            result['pos_mean_shuf'] = pos_shuf['overall']['mean']
            result['pos_5prime_shuf'] = pos_shuf['overall']['five_prime_frac']
            result['pos_3prime_shuf'] = pos_shuf['overall']['three_prime_frac']

        # Positional enrichment
        if pos_real.get('overall') and pos_shuf.get('overall'):
            r5 = pos_real['overall']['five_prime_frac']
            s5 = pos_shuf['overall']['five_prime_frac']
            r3 = pos_real['overall']['three_prime_frac']
            s3 = pos_shuf['overall']['three_prime_frac']
            result['pos_5prime_enrich'] = r5 / s5 if s5 > 0 else np.nan
            result['pos_3prime_enrich'] = r3 / s3 if s3 > 0 else np.nan

        # Quality analysis
        qual_real = analyze_quality_metrics(real_df)
        qual_shuf = analyze_quality_metrics(shuf_df)

        for metric in ['pident', 'length']:
            if metric in qual_real:
                result[f'{metric}_mean_real'] = qual_real[metric]['mean']
                result[f'{metric}_median_real'] = qual_real[metric]['median']
            if metric in qual_shuf:
                result[f'{metric}_mean_shuf'] = qual_shuf[metric]['mean']
                result[f'{metric}_median_shuf'] = qual_shuf[metric]['median']

        if 'pct_with_gaps' in qual_real:
            result['pct_gaps_real'] = qual_real['pct_with_gaps']
        if 'pct_with_gaps' in qual_shuf:
            result['pct_gaps_shuf'] = qual_shuf['pct_with_gaps']

        # TE family analysis
        te_real = analyze_te_families(real_df)
        te_shuf = analyze_te_families(shuf_df)

        if te_real:
            result['n_te_families_real'] = te_real['n_families']
            result['te_entropy_real'] = te_real['shannon_entropy']
            result['te_top5_conc_real'] = te_real['top5_concentration']
            if te_real.get('top_families'):
                result['top_te_family_real'] = te_real['top_families'][0]['family']
                result['top_te_pct_real'] = te_real['top_families'][0]['pct']

        if te_shuf:
            result['n_te_families_shuf'] = te_shuf['n_families']
            result['te_entropy_shuf'] = te_shuf['shannon_entropy']

        all_results.append(result)
        print(f"  -> {len(real_df):,} real, {len(shuf_df):,} shuffled, ratio={result['real_shuf_ratio']:.2f}x")

    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(output_dir / 'analysis_summary.tsv', sep='\t', index=False, float_format='%.4f')
    print(f"\nSaved analysis summary: {output_dir / 'analysis_summary.tsv'}")

    # Generate visualizations
    print("\n" + "=" * 70)
    print("GENERATING VISUALIZATIONS")
    print("=" * 70)

    generate_visualizations(results_df, sweep_dir, output_dir)

    return 0


def generate_visualizations(results_df: pd.DataFrame, sweep_dir: Path, output_dir: Path):
    """Generate comprehensive visualizations."""

    # =========================================================================
    # FIGURE 1: Parameter effect summary
    # =========================================================================
    print("\nGenerating Figure 1: Parameter effects...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Panel A: Hit counts by word_size and dust
    ax = axes[0, 0]
    for ws in [7, 11]:
        for dust in ['yes', 'no']:
            subset = results_df[(results_df['word_size'] == ws) & (results_df['dust'] == dust)]
            label = f'ws={ws}, dust={dust}'
            color = 'steelblue' if dust == 'yes' else 'coral'
            marker = 'o' if ws == 7 else 's'
            ax.scatter([f"go{r['gapopen']}_ge{r['gapextend']}" for _, r in subset.iterrows()],
                      subset['n_hits_real'], label=label, color=color, marker=marker, s=80, alpha=0.7)
    ax.set_xlabel('Gap parameters')
    ax.set_ylabel('Total hits (real)')
    ax.set_title('A. Hit counts by parameters')
    ax.legend(fontsize=8)
    ax.tick_params(axis='x', rotation=45)

    # Panel B: Real/Shuffled ratio
    ax = axes[0, 1]
    for dust in ['yes', 'no']:
        subset = results_df[results_df['dust'] == dust]
        color = 'steelblue' if dust == 'yes' else 'coral'
        ws7 = subset[subset['word_size'] == 7]['real_shuf_ratio'].values
        ws11 = subset[subset['word_size'] == 11]['real_shuf_ratio'].values
        x = np.arange(len(ws7))
        width = 0.35
        if dust == 'yes':
            ax.bar(x - width/2, ws7, width, label=f'dust={dust}, ws=7', color='steelblue', alpha=0.7)
            ax.bar(x + width/2, ws11, width, label=f'dust={dust}, ws=11', color='steelblue', alpha=0.4)
        else:
            ax.bar(x - width/2 + 0.02, ws7, width, label=f'dust={dust}, ws=7', color='coral', alpha=0.7)
            ax.bar(x + width/2 + 0.02, ws11, width, label=f'dust={dust}, ws=11', color='coral', alpha=0.4)
    ax.axhline(1, color='gray', linestyle='--', linewidth=2)
    ax.set_xlabel('Parameter combination')
    ax.set_ylabel('Real/Shuffled ratio')
    ax.set_title('B. Signal vs noise ratio')
    ax.legend(fontsize=7)

    # Panel C: Positional enrichment (3' end)
    ax = axes[0, 2]
    for ws in [7, 11]:
        for dust in ['yes', 'no']:
            subset = results_df[(results_df['word_size'] == ws) & (results_df['dust'] == dust)]
            label = f'ws={ws}, dust={dust}'
            color = 'green' if dust == 'yes' else 'purple'
            marker = 'o' if ws == 7 else 's'
            ax.scatter(subset['pos_5prime_enrich'], subset['pos_3prime_enrich'],
                      label=label, color=color, marker=marker, s=80, alpha=0.7)
    ax.axhline(1, color='gray', linestyle='--', linewidth=1)
    ax.axvline(1, color='gray', linestyle='--', linewidth=1)
    ax.set_xlabel("5' enrichment (real/shuf)")
    ax.set_ylabel("3' enrichment (real/shuf)")
    ax.set_title("C. Positional enrichment")
    ax.legend(fontsize=8)

    # Panel D: Mean alignment length
    ax = axes[1, 0]
    dust_yes = results_df[results_df['dust'] == 'yes']
    dust_no = results_df[results_df['dust'] == 'no']
    x = np.arange(len(dust_yes))
    width = 0.35
    ax.bar(x - width/2, dust_yes['length_mean_real'], width, label='DUST=yes', color='steelblue')
    ax.bar(x + width/2, dust_no['length_mean_real'], width, label='DUST=no', color='coral')
    ax.set_xlabel('Parameter combination')
    ax.set_ylabel('Mean alignment length (bp)')
    ax.set_title('D. Alignment length')
    ax.legend()

    # Panel E: Percent identity
    ax = axes[1, 1]
    ax.bar(x - width/2, dust_yes['pident_mean_real'], width, label='DUST=yes', color='steelblue')
    ax.bar(x + width/2, dust_no['pident_mean_real'], width, label='DUST=no', color='coral')
    ax.set_xlabel('Parameter combination')
    ax.set_ylabel('Mean percent identity')
    ax.set_title('E. Alignment quality')
    ax.legend()

    # Panel F: TE family diversity
    ax = axes[1, 2]
    ax.bar(x - width/2, dust_yes['te_entropy_real'], width, label='DUST=yes', color='steelblue')
    ax.bar(x + width/2, dust_no['te_entropy_real'], width, label='DUST=no', color='coral')
    ax.set_xlabel('Parameter combination')
    ax.set_ylabel('Shannon entropy (TE families)')
    ax.set_title('F. TE family diversity')
    ax.legend()

    fig.suptitle('Parameter Sweep: Effect on BLAST Results', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'param_effects_summary.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: param_effects_summary.png")

    # =========================================================================
    # FIGURE 2: Positional distributions (selected params)
    # =========================================================================
    print("\nGenerating Figure 2: Positional distributions...")

    # Select key parameter combinations to compare
    key_combos = [
        ('ws7_go2_ge1_p-1_dustyes', 'ws=7, optimal, DUST=yes'),
        ('ws7_go2_ge1_p-1_dustno', 'ws=7, optimal, DUST=no'),
        ('ws11_go2_ge1_p-1_dustyes', 'ws=11, optimal, DUST=yes'),
        ('ws11_go2_ge1_p-1_dustno', 'ws=11, optimal, DUST=no'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    n_bins = 20

    for idx, (combo, label) in enumerate(key_combos):
        ax = axes[idx]

        real_file = sweep_dir / f"{combo}_real.tsv"
        shuf_file = sweep_dir / f"{combo}_shuffled.tsv"

        if not real_file.exists():
            ax.set_title(f'{label}\n(no data)')
            continue

        real_df = load_blast_results(real_file)
        shuf_df = load_blast_results(shuf_file)

        if len(real_df) == 0:
            continue

        # Compute distributions
        real_counts, bin_edges = np.histogram(real_df['rel_mid'], bins=n_bins, range=(0, 1))
        shuf_counts, _ = np.histogram(shuf_df['rel_mid'], bins=n_bins, range=(0, 1))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        real_density = real_counts / real_counts.sum()
        shuf_density = shuf_counts / shuf_counts.sum() if shuf_counts.sum() > 0 else shuf_counts

        # Plot
        ax.fill_between(bin_centers, real_density, alpha=0.3, color='steelblue', label='Real')
        ax.plot(bin_centers, real_density, color='steelblue', linewidth=2)
        ax.fill_between(bin_centers, shuf_density, alpha=0.3, color='coral', label='Shuffled')
        ax.plot(bin_centers, shuf_density, color='coral', linewidth=2)
        ax.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)

        ax.set_xlabel("Relative position (0=5' end, 1=3' end)")
        ax.set_ylabel('Density')
        ax.set_title(f'{label}\n(n={len(real_df):,} real, {len(shuf_df):,} shuf)')
        ax.legend(fontsize=9)
        ax.set_xlim(0, 1)

    fig.suptitle('Positional Distribution: Key Parameter Comparisons', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'position_distributions.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_distributions.png")

    # =========================================================================
    # FIGURE 3: Word size comparison
    # =========================================================================
    print("\nGenerating Figure 3: Word size comparison...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: Hits by word size
    ax = axes[0, 0]
    ws7 = results_df[results_df['word_size'] == 7]
    ws11 = results_df[results_df['word_size'] == 11]
    ax.scatter(ws7['n_hits_real'], ws11['n_hits_real'].values,
              c=['steelblue' if d == 'yes' else 'coral' for d in ws7['dust']], s=100, alpha=0.7)
    max_val = max(ws7['n_hits_real'].max(), ws11['n_hits_real'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1)
    ax.set_xlabel('Hits with word_size=7')
    ax.set_ylabel('Hits with word_size=11')
    ax.set_title('A. Total hits: ws=7 vs ws=11')

    # Panel B: Real/Shuf ratio by word size
    ax = axes[0, 1]
    ax.scatter(ws7['real_shuf_ratio'], ws11['real_shuf_ratio'].values,
              c=['steelblue' if d == 'yes' else 'coral' for d in ws7['dust']], s=100, alpha=0.7)
    max_val = max(ws7['real_shuf_ratio'].max(), ws11['real_shuf_ratio'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1)
    ax.set_xlabel('Real/Shuf ratio (ws=7)')
    ax.set_ylabel('Real/Shuf ratio (ws=11)')
    ax.set_title('B. Signal ratio: ws=7 vs ws=11')

    # Panel C: Length by word size
    ax = axes[1, 0]
    ax.scatter(ws7['length_mean_real'], ws11['length_mean_real'].values,
              c=['steelblue' if d == 'yes' else 'coral' for d in ws7['dust']], s=100, alpha=0.7)
    max_val = max(ws7['length_mean_real'].max(), ws11['length_mean_real'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1)
    ax.set_xlabel('Mean length (ws=7)')
    ax.set_ylabel('Mean length (ws=11)')
    ax.set_title('C. Alignment length: ws=7 vs ws=11')

    # Panel D: Summary bar chart
    ax = axes[1, 1]
    metrics = ['n_hits_real', 'real_shuf_ratio', 'length_mean_real', 'pident_mean_real']
    metric_labels = ['Hits', 'R/S Ratio', 'Length', 'Identity']
    x = np.arange(len(metrics))
    width = 0.35

    # Normalize to ws=7 values for comparison
    ws7_means = [ws7[m].mean() for m in metrics]
    ws11_means = [ws11[m].mean() for m in metrics]
    ws7_norm = [1.0] * len(metrics)
    ws11_norm = [ws11_means[i] / ws7_means[i] if ws7_means[i] > 0 else 0 for i in range(len(metrics))]

    ax.bar(x - width/2, ws7_norm, width, label='ws=7', color='steelblue')
    ax.bar(x + width/2, ws11_norm, width, label='ws=11', color='coral')
    ax.axhline(1, color='gray', linestyle='--', linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels)
    ax.set_ylabel('Relative to ws=7')
    ax.set_title('D. Summary: ws=11 relative to ws=7')
    ax.legend()

    fig.suptitle('Word Size Comparison (7 vs 11)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'wordsize_comparison.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: wordsize_comparison.png")

    # =========================================================================
    # FIGURE 4: TE Family Analysis
    # =========================================================================
    print("\nGenerating Figure 4: TE family analysis...")

    # Load data for optimal params with dust=yes
    optimal_combo = 'ws7_go2_ge1_p-1_dustyes'
    real_file = sweep_dir / f"{optimal_combo}_real.tsv"
    shuf_file = sweep_dir / f"{optimal_combo}_shuffled.tsv"

    if real_file.exists():
        real_df = load_blast_results(real_file)
        shuf_df = load_blast_results(shuf_file)

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Panel A: Top TE families (real)
        ax = axes[0, 0]
        real_counts = real_df['te_family'].value_counts().head(15)
        ax.barh(range(len(real_counts)), real_counts.values, color='steelblue')
        ax.set_yticks(range(len(real_counts)))
        ax.set_yticklabels(real_counts.index)
        ax.set_xlabel('Number of hits')
        ax.set_title('A. Top 15 TE families (Real UTRs)')
        ax.invert_yaxis()

        # Panel B: Real vs Shuffled enrichment
        ax = axes[0, 1]
        real_fam = real_df['te_family'].value_counts()
        shuf_fam = shuf_df['te_family'].value_counts()

        # Calculate enrichment for top families
        top_fams = real_fam.head(15).index
        enrichment = []
        for fam in top_fams:
            r = real_fam.get(fam, 0)
            s = shuf_fam.get(fam, 1)
            # Normalize by total hits
            r_norm = r / len(real_df)
            s_norm = s / len(shuf_df) if len(shuf_df) > 0 else 0.001
            enrichment.append(r_norm / s_norm if s_norm > 0 else r_norm / 0.001)

        colors = ['steelblue' if e > 1 else 'coral' for e in enrichment]
        ax.barh(range(len(top_fams)), enrichment, color=colors)
        ax.axvline(1, color='gray', linestyle='--', linewidth=2)
        ax.set_yticks(range(len(top_fams)))
        ax.set_yticklabels(top_fams)
        ax.set_xlabel('Enrichment (Real/Shuffled)')
        ax.set_title('B. TE Family Enrichment')
        ax.invert_yaxis()

        # Panel C: TE diversity by parameter
        ax = axes[1, 0]
        for dust in ['yes', 'no']:
            subset = results_df[results_df['dust'] == dust]
            ax.scatter(subset['n_te_families_real'], subset['te_entropy_real'],
                      c='steelblue' if dust == 'yes' else 'coral',
                      label=f'DUST={dust}', s=100, alpha=0.7)
        ax.set_xlabel('Number of TE families')
        ax.set_ylabel('Shannon entropy')
        ax.set_title('C. TE Family Diversity by Parameters')
        ax.legend()

        # Panel D: Top family concentration
        ax = axes[1, 1]
        for ws in [7, 11]:
            subset = results_df[results_df['word_size'] == ws]
            color = 'steelblue' if ws == 7 else 'coral'
            ax.scatter(subset['te_top5_conc_real'], subset['real_shuf_ratio'],
                      c=color, label=f'ws={ws}', s=100, alpha=0.7,
                      marker='o' if ws == 7 else 's')
        ax.set_xlabel('Top 5 TE families (% of hits)')
        ax.set_ylabel('Real/Shuffled ratio')
        ax.set_title('D. TE Concentration vs Signal')
        ax.legend()

        fig.suptitle('TE Family Analysis', fontsize=14, fontweight='bold')
        plt.tight_layout()
        fig.savefig(output_dir / 'te_family_analysis.png', dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        print(f"  Saved: te_family_analysis.png")

    print("\nVisualization complete!")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sweep-dir', type=str, required=True,
                        help='Directory containing sweep results')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: figures/param_sweep_analysis)')

    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir)
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path('figures/param_sweep_analysis')

    return run_full_analysis(sweep_dir, output_dir)


if __name__ == '__main__':
    sys.exit(main())
