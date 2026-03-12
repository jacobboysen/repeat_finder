#!/usr/bin/env python3
"""
Analyze positional distribution of TE hits along UTRs across BLAST parameter sets.

Compares real vs shuffled sequences, with and without DUST filtering,
stratified by UTR length.

Usage:
    python scripts/analyze_param_position_distribution.py
"""

import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

# Add project root for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir
from utils.blast_io import BLAST_COLUMNS

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Color palette
COLORS = {
    'real': '#2171b5',       # Blue
    'shuffled': '#cb181d',   # Red
    'dustyes': '#238b45',    # Green
    'dustno': '#6a51a3',     # Purple
}

# UTR length bins
UTR_LEN_BINS = [
    (0, 200, '<200bp'),
    (200, 500, '200-500bp'),
    (500, 1000, '500-1kb'),
    (1000, 2000, '1-2kb'),
    (2000, 100000, '>2kb'),
]

# Parameters to analyze
PARAM_COMBOS = [
    ('go2_ge1_r1p-1', 'go=2, ge=1, p=-1 (optimal)'),
    ('go2_ge1_r1p-2', 'go=2, ge=1, p=-2'),
    ('go2_ge2_r1p-1', 'go=2, ge=2, p=-1'),
    ('go5_ge2_r1p-1', 'go=5, ge=2, p=-1'),
    ('go10_ge2_r1p-1', 'go=10, ge=2, p=-1'),
]


def load_blast_positions(filepath: Path) -> pd.DataFrame:
    """Load BLAST results and extract positional information."""
    if not filepath.exists():
        print(f"  Warning: File not found: {filepath}")
        return pd.DataFrame()

    # Read BLAST output (no header)
    df = pd.read_csv(filepath, sep='\t', header=None, names=BLAST_COLUMNS[:16])

    # Calculate relative positions
    df['rel_start'] = df[['qstart', 'qend']].min(axis=1) / df['qlen']
    df['rel_end'] = df[['qstart', 'qend']].max(axis=1) / df['qlen']
    df['rel_mid'] = (df['rel_start'] + df['rel_end']) / 2
    df['hit_len'] = df['length']

    return df[['qseqid', 'qlen', 'rel_start', 'rel_end', 'rel_mid', 'hit_len', 'pident', 'evalue']]


def compute_position_distribution(df: pd.DataFrame, n_bins: int = 20,
                                   utr_len_range: Optional[Tuple[int, int]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute position distribution (density histogram).

    Args:
        df: DataFrame with rel_mid column
        n_bins: Number of bins for distribution
        utr_len_range: Optional (min, max) UTR length filter

    Returns:
        bin_centers, density (normalized to sum to 1)
    """
    if len(df) == 0:
        return np.linspace(0.025, 0.975, n_bins), np.zeros(n_bins)

    data = df.copy()

    # Filter by UTR length if specified
    if utr_len_range:
        lo, hi = utr_len_range
        data = data[(data['qlen'] >= lo) & (data['qlen'] < hi)]

    if len(data) == 0:
        return np.linspace(0.025, 0.975, n_bins), np.zeros(n_bins)

    # Compute histogram
    counts, bin_edges = np.histogram(data['rel_mid'], bins=n_bins, range=(0, 1))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Normalize to density
    density = counts / counts.sum() if counts.sum() > 0 else counts

    return bin_centers, density


def analyze_all_files(results_dir: Path) -> Dict:
    """
    Load and analyze all parameter comparison files.

    Returns:
        Dict with structure: {param_combo: {dust: {seq_type: DataFrame}}}
    """
    data = {}

    for param_code, param_name in PARAM_COMBOS:
        data[param_code] = {'name': param_name}

        for dust in ['yes', 'no']:
            data[param_code][f'dust{dust}'] = {}

            for seq_type in ['real', 'shuffled']:
                filename = f"{param_code}_dust{dust}_{seq_type}.tsv"
                filepath = results_dir / filename

                print(f"Loading {filename}...", end=" ")
                df = load_blast_positions(filepath)
                data[param_code][f'dust{dust}'][seq_type] = df
                print(f"{len(df):,} hits")

    return data


def plot_position_distributions(data: Dict, output_dir: Path):
    """Create comprehensive position distribution visualizations."""
    output_dir.mkdir(parents=True, exist_ok=True)

    n_bins = 20

    # =========================================================================
    # FIGURE 1: Real vs Shuffled comparison, DUST=yes only, all params
    # =========================================================================
    print("\nGenerating Figure 1: Real vs Shuffled (DUST=yes)...")

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    for idx, (param_code, param_name) in enumerate(PARAM_COMBOS):
        ax = axes[idx]

        real_df = data[param_code]['dustyes']['real']
        shuf_df = data[param_code]['dustyes']['shuffled']

        # Compute distributions
        bins, real_dist = compute_position_distribution(real_df, n_bins)
        _, shuf_dist = compute_position_distribution(shuf_df, n_bins)

        # Plot
        ax.fill_between(bins, real_dist, alpha=0.3, color=COLORS['real'], label='Real')
        ax.plot(bins, real_dist, color=COLORS['real'], linewidth=2)

        ax.fill_between(bins, shuf_dist, alpha=0.3, color=COLORS['shuffled'], label='Shuffled')
        ax.plot(bins, shuf_dist, color=COLORS['shuffled'], linewidth=2)

        ax.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7, label='Uniform')

        ax.set_xlabel("Relative position (0=5' end, 1=3' end)")
        ax.set_ylabel('Density')
        ax.set_title(f'{data[param_code]["name"]}\n(n={len(real_df):,} real, {len(shuf_df):,} shuf)')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, max(real_dist.max(), shuf_dist.max()) * 1.2)

        if idx == 0:
            ax.legend(loc='upper right', fontsize=9)

    # Hide last subplot if odd number
    if len(PARAM_COMBOS) < 6:
        axes[-1].axis('off')

    fig.suptitle('TE Hit Position Distribution: Real vs Shuffled (DUST=yes)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'position_real_vs_shuffled_dustyes.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_real_vs_shuffled_dustyes.png")

    # =========================================================================
    # FIGURE 2: DUST effect comparison (Real sequences only)
    # =========================================================================
    print("\nGenerating Figure 2: DUST effect (Real only)...")

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    for idx, (param_code, param_name) in enumerate(PARAM_COMBOS):
        ax = axes[idx]

        dustyes_df = data[param_code]['dustyes']['real']
        dustno_df = data[param_code]['dustno']['real']

        bins, dustyes_dist = compute_position_distribution(dustyes_df, n_bins)
        _, dustno_dist = compute_position_distribution(dustno_df, n_bins)

        ax.fill_between(bins, dustyes_dist, alpha=0.3, color=COLORS['dustyes'], label='DUST=yes')
        ax.plot(bins, dustyes_dist, color=COLORS['dustyes'], linewidth=2)

        ax.fill_between(bins, dustno_dist, alpha=0.3, color=COLORS['dustno'], label='DUST=no')
        ax.plot(bins, dustno_dist, color=COLORS['dustno'], linewidth=2)

        ax.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)

        ax.set_xlabel("Relative position (0=5' end, 1=3' end)")
        ax.set_ylabel('Density')
        ax.set_title(f'{data[param_code]["name"]}\n(DUST=yes: {len(dustyes_df):,}, DUST=no: {len(dustno_df):,})')
        ax.set_xlim(0, 1)

        if idx == 0:
            ax.legend(loc='upper right', fontsize=9)

    if len(PARAM_COMBOS) < 6:
        axes[-1].axis('off')

    fig.suptitle('DUST Filtering Effect on Position Distribution (Real UTRs)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'position_dust_effect_real.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_dust_effect_real.png")

    # =========================================================================
    # FIGURE 3: Position distribution by UTR length (optimal params, DUST=yes)
    # =========================================================================
    print("\nGenerating Figure 3: By UTR length (optimal params)...")

    param_code = 'go2_ge1_r1p-1'  # Optimal params

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.flatten()

    for idx, (lo, hi, label) in enumerate(UTR_LEN_BINS):
        ax = axes[idx]

        real_df = data[param_code]['dustyes']['real']
        shuf_df = data[param_code]['dustyes']['shuffled']

        bins, real_dist = compute_position_distribution(real_df, n_bins, (lo, hi))
        _, shuf_dist = compute_position_distribution(shuf_df, n_bins, (lo, hi))

        # Count hits in this bin
        n_real = len(real_df[(real_df['qlen'] >= lo) & (real_df['qlen'] < hi)])
        n_shuf = len(shuf_df[(shuf_df['qlen'] >= lo) & (shuf_df['qlen'] < hi)])

        ax.fill_between(bins, real_dist, alpha=0.3, color=COLORS['real'], label='Real')
        ax.plot(bins, real_dist, color=COLORS['real'], linewidth=2)

        ax.fill_between(bins, shuf_dist, alpha=0.3, color=COLORS['shuffled'], label='Shuffled')
        ax.plot(bins, shuf_dist, color=COLORS['shuffled'], linewidth=2)

        ax.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)

        ax.set_xlabel("Relative position (0=5' end, 1=3' end)")
        ax.set_ylabel('Density')
        ax.set_title(f'UTR length: {label}\n(n={n_real:,} real, {n_shuf:,} shuf)')
        ax.set_xlim(0, 1)

        if idx == 0:
            ax.legend(loc='upper right', fontsize=9)

    axes[-1].axis('off')

    fig.suptitle(f'Position Distribution by UTR Length\n(BLAST params: {data[param_code]["name"]}, DUST=yes)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'position_by_utr_length.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_by_utr_length.png")

    # =========================================================================
    # FIGURE 4: Comprehensive dashboard - single figure with all key comparisons
    # =========================================================================
    print("\nGenerating Figure 4: Comprehensive dashboard...")

    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)

    param_code = 'go2_ge1_r1p-1'  # Focus on optimal params

    # Row 1: Real vs Shuffled for different DUST settings
    # Panel A: DUST=yes, Real vs Shuffled
    ax1 = fig.add_subplot(gs[0, 0])
    real_df = data[param_code]['dustyes']['real']
    shuf_df = data[param_code]['dustyes']['shuffled']
    bins, real_dist = compute_position_distribution(real_df, n_bins)
    _, shuf_dist = compute_position_distribution(shuf_df, n_bins)

    ax1.fill_between(bins, real_dist, alpha=0.3, color=COLORS['real'])
    ax1.plot(bins, real_dist, color=COLORS['real'], linewidth=2, label='Real')
    ax1.fill_between(bins, shuf_dist, alpha=0.3, color=COLORS['shuffled'])
    ax1.plot(bins, shuf_dist, color=COLORS['shuffled'], linewidth=2, label='Shuffled')
    ax1.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax1.set_xlabel("Position (0=5', 1=3')")
    ax1.set_ylabel('Density')
    ax1.set_title('A. DUST=yes: Real vs Shuffled')
    ax1.legend(fontsize=8)
    ax1.set_xlim(0, 1)

    # Panel B: DUST=no, Real vs Shuffled
    ax2 = fig.add_subplot(gs[0, 1])
    real_df_no = data[param_code]['dustno']['real']
    shuf_df_no = data[param_code]['dustno']['shuffled']
    bins, real_dist_no = compute_position_distribution(real_df_no, n_bins)
    _, shuf_dist_no = compute_position_distribution(shuf_df_no, n_bins)

    ax2.fill_between(bins, real_dist_no, alpha=0.3, color=COLORS['real'])
    ax2.plot(bins, real_dist_no, color=COLORS['real'], linewidth=2, label='Real')
    ax2.fill_between(bins, shuf_dist_no, alpha=0.3, color=COLORS['shuffled'])
    ax2.plot(bins, shuf_dist_no, color=COLORS['shuffled'], linewidth=2, label='Shuffled')
    ax2.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax2.set_xlabel("Position (0=5', 1=3')")
    ax2.set_ylabel('Density')
    ax2.set_title('B. DUST=no: Real vs Shuffled')
    ax2.legend(fontsize=8)
    ax2.set_xlim(0, 1)

    # Panel C: DUST effect on Real
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.fill_between(bins, real_dist, alpha=0.3, color=COLORS['dustyes'])
    ax3.plot(bins, real_dist, color=COLORS['dustyes'], linewidth=2, label='DUST=yes')
    ax3.fill_between(bins, real_dist_no, alpha=0.3, color=COLORS['dustno'])
    ax3.plot(bins, real_dist_no, color=COLORS['dustno'], linewidth=2, label='DUST=no')
    ax3.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax3.set_xlabel("Position (0=5', 1=3')")
    ax3.set_ylabel('Density')
    ax3.set_title('C. Real UTRs: DUST Effect')
    ax3.legend(fontsize=8)
    ax3.set_xlim(0, 1)

    # Panel D: Real/Shuffled enrichment ratio
    ax4 = fig.add_subplot(gs[0, 3])
    ratio_yes = real_dist / (shuf_dist + 1e-10)
    ratio_no = real_dist_no / (shuf_dist_no + 1e-10)

    ax4.plot(bins, ratio_yes, color=COLORS['dustyes'], linewidth=2, marker='o', markersize=4, label='DUST=yes')
    ax4.plot(bins, ratio_no, color=COLORS['dustno'], linewidth=2, marker='s', markersize=4, label='DUST=no')
    ax4.axhline(1, color='gray', linestyle='--', linewidth=2)
    ax4.set_xlabel("Position (0=5', 1=3')")
    ax4.set_ylabel('Real/Shuffled Ratio')
    ax4.set_title('D. Enrichment Ratio by Position')
    ax4.legend(fontsize=8)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0.5, 2.0)

    # Row 2: By UTR length (DUST=yes)
    for idx, (lo, hi, label) in enumerate(UTR_LEN_BINS[:4]):
        ax = fig.add_subplot(gs[1, idx])

        bins, real_dist = compute_position_distribution(real_df, n_bins, (lo, hi))
        _, shuf_dist = compute_position_distribution(shuf_df, n_bins, (lo, hi))

        n_real = len(real_df[(real_df['qlen'] >= lo) & (real_df['qlen'] < hi)])
        n_shuf = len(shuf_df[(shuf_df['qlen'] >= lo) & (shuf_df['qlen'] < hi)])

        ax.fill_between(bins, real_dist, alpha=0.3, color=COLORS['real'])
        ax.plot(bins, real_dist, color=COLORS['real'], linewidth=2)
        ax.fill_between(bins, shuf_dist, alpha=0.3, color=COLORS['shuffled'])
        ax.plot(bins, shuf_dist, color=COLORS['shuffled'], linewidth=2)
        ax.axhline(1/n_bins, color='gray', linestyle='--', linewidth=1, alpha=0.7)

        panel_letter = chr(ord('E') + idx)
        ax.set_xlabel("Position (0=5', 1=3')")
        ax.set_ylabel('Density')
        ax.set_title(f'{panel_letter}. UTR {label} (n={n_real:,}/{n_shuf:,})')
        ax.set_xlim(0, 1)

    # Row 3: Heatmaps and summary
    # Panel I: Real position heatmap by UTR length
    ax_heat1 = fig.add_subplot(gs[2, 0:2])
    heatmap_data_real = []
    for lo, hi, label in UTR_LEN_BINS:
        _, dist = compute_position_distribution(real_df, n_bins, (lo, hi))
        heatmap_data_real.append(dist)

    heatmap_real = np.array(heatmap_data_real)
    im1 = ax_heat1.imshow(heatmap_real, cmap='YlOrRd', aspect='auto',
                          extent=[0, 1, len(UTR_LEN_BINS)-0.5, -0.5])
    ax_heat1.set_xlabel("Relative Position (0=5' end, 1=3' end)")
    ax_heat1.set_yticks(range(len(UTR_LEN_BINS)))
    ax_heat1.set_yticklabels([b[2] for b in UTR_LEN_BINS])
    ax_heat1.set_ylabel('UTR Length')
    ax_heat1.set_title('I. Real Hits: Position Density by UTR Length')
    plt.colorbar(im1, ax=ax_heat1, label='Density', shrink=0.8)

    # Panel J: Shuffled position heatmap by UTR length
    ax_heat2 = fig.add_subplot(gs[2, 2:4])
    heatmap_data_shuf = []
    for lo, hi, label in UTR_LEN_BINS:
        _, dist = compute_position_distribution(shuf_df, n_bins, (lo, hi))
        heatmap_data_shuf.append(dist)

    heatmap_shuf = np.array(heatmap_data_shuf)
    im2 = ax_heat2.imshow(heatmap_shuf, cmap='YlOrRd', aspect='auto',
                          extent=[0, 1, len(UTR_LEN_BINS)-0.5, -0.5])
    ax_heat2.set_xlabel("Relative Position (0=5' end, 1=3' end)")
    ax_heat2.set_yticks(range(len(UTR_LEN_BINS)))
    ax_heat2.set_yticklabels([b[2] for b in UTR_LEN_BINS])
    ax_heat2.set_ylabel('UTR Length')
    ax_heat2.set_title('J. Shuffled Hits: Position Density by UTR Length')
    plt.colorbar(im2, ax=ax_heat2, label='Density', shrink=0.8)

    fig.suptitle(f'TE Hit Position Analysis: Optimal Parameters ({data[param_code]["name"]})',
                 fontsize=16, fontweight='bold', y=0.98)

    fig.savefig(output_dir / 'position_dashboard.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_dashboard.png")

    # =========================================================================
    # FIGURE 5: Parameter comparison summary
    # =========================================================================
    print("\nGenerating Figure 5: Parameter comparison summary...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Collect summary stats for each parameter
    summary_data = []

    for param_code, param_name in PARAM_COMBOS:
        for dust in ['yes', 'no']:
            real_df = data[param_code][f'dust{dust}']['real']
            shuf_df = data[param_code][f'dust{dust}']['shuffled']

            bins, real_dist = compute_position_distribution(real_df, n_bins)
            _, shuf_dist = compute_position_distribution(shuf_df, n_bins)

            # 5' quarter enrichment
            five_prime_real = real_dist[:5].sum()
            five_prime_shuf = shuf_dist[:5].sum()
            five_prime_enrich = five_prime_real / (five_prime_shuf + 1e-10)

            # 3' quarter enrichment
            three_prime_real = real_dist[-5:].sum()
            three_prime_shuf = shuf_dist[-5:].sum()
            three_prime_enrich = three_prime_real / (three_prime_shuf + 1e-10)

            # Middle enrichment
            mid_real = real_dist[5:15].sum()
            mid_shuf = shuf_dist[5:15].sum()
            mid_enrich = mid_real / (mid_shuf + 1e-10)

            # Distribution uniformity (KL-divergence from uniform)
            uniform = np.ones(n_bins) / n_bins
            kl_real = np.sum(real_dist * np.log((real_dist + 1e-10) / uniform))
            kl_shuf = np.sum(shuf_dist * np.log((shuf_dist + 1e-10) / uniform))

            summary_data.append({
                'param': param_code,
                'param_name': param_name,
                'dust': dust,
                'n_real': len(real_df),
                'n_shuf': len(shuf_df),
                'real_shuf_ratio': len(real_df) / (len(shuf_df) + 1),
                'five_prime_enrich': five_prime_enrich,
                'three_prime_enrich': three_prime_enrich,
                'mid_enrich': mid_enrich,
                'kl_real': kl_real,
                'kl_shuf': kl_shuf,
            })

    summary_df = pd.DataFrame(summary_data)

    # Panel A: Total hits by parameter and DUST
    ax = axes[0, 0]
    x = np.arange(len(PARAM_COMBOS))
    width = 0.35

    dustyes_hits = summary_df[summary_df['dust'] == 'yes'].groupby('param')['n_real'].first()
    dustno_hits = summary_df[summary_df['dust'] == 'no'].groupby('param')['n_real'].first()

    # Reorder to match PARAM_COMBOS order
    dustyes_vals = [dustyes_hits.get(p[0], 0) for p in PARAM_COMBOS]
    dustno_vals = [dustno_hits.get(p[0], 0) for p in PARAM_COMBOS]

    ax.bar(x - width/2, dustyes_vals, width, label='DUST=yes', color=COLORS['dustyes'], edgecolor='black')
    ax.bar(x + width/2, dustno_vals, width, label='DUST=no', color=COLORS['dustno'], edgecolor='black')
    ax.set_xlabel('Parameter Set')
    ax.set_ylabel('Total Hits (Real)')
    ax.set_title('A. Total Hits by Parameter & DUST')
    ax.set_xticks(x)
    ax.set_xticklabels([p[0].replace('_', '\n') for p in PARAM_COMBOS], fontsize=8)
    ax.legend()

    # Panel B: Real/Shuffled ratio
    ax = axes[0, 1]
    dustyes_ratio = summary_df[summary_df['dust'] == 'yes'].groupby('param')['real_shuf_ratio'].first()
    dustno_ratio = summary_df[summary_df['dust'] == 'no'].groupby('param')['real_shuf_ratio'].first()

    dustyes_r = [dustyes_ratio.get(p[0], 0) for p in PARAM_COMBOS]
    dustno_r = [dustno_ratio.get(p[0], 0) for p in PARAM_COMBOS]

    ax.bar(x - width/2, dustyes_r, width, label='DUST=yes', color=COLORS['dustyes'], edgecolor='black')
    ax.bar(x + width/2, dustno_r, width, label='DUST=no', color=COLORS['dustno'], edgecolor='black')
    ax.axhline(1, color='gray', linestyle='--', linewidth=2)
    ax.set_xlabel('Parameter Set')
    ax.set_ylabel('Real/Shuffled Ratio')
    ax.set_title('B. Real/Shuffled Hit Ratio')
    ax.set_xticks(x)
    ax.set_xticklabels([p[0].replace('_', '\n') for p in PARAM_COMBOS], fontsize=8)
    ax.legend()

    # Panel C: 5' vs 3' enrichment (DUST=yes only)
    ax = axes[1, 0]
    dustyes_df = summary_df[summary_df['dust'] == 'yes']

    five_prime = [dustyes_df[dustyes_df['param'] == p[0]]['five_prime_enrich'].values[0] for p in PARAM_COMBOS]
    three_prime = [dustyes_df[dustyes_df['param'] == p[0]]['three_prime_enrich'].values[0] for p in PARAM_COMBOS]

    ax.bar(x - width/2, five_prime, width, label="5' end (0-25%)", color='#74c476', edgecolor='black')
    ax.bar(x + width/2, three_prime, width, label="3' end (75-100%)", color='#fd8d3c', edgecolor='black')
    ax.axhline(1, color='gray', linestyle='--', linewidth=2)
    ax.set_xlabel('Parameter Set')
    ax.set_ylabel('Real/Shuffled Enrichment')
    ax.set_title("C. Positional Enrichment (DUST=yes)")
    ax.set_xticks(x)
    ax.set_xticklabels([p[0].replace('_', '\n') for p in PARAM_COMBOS], fontsize=8)
    ax.legend()

    # Panel D: Distribution uniformity (KL divergence)
    ax = axes[1, 1]
    kl_real_yes = [dustyes_df[dustyes_df['param'] == p[0]]['kl_real'].values[0] for p in PARAM_COMBOS]
    kl_shuf_yes = [dustyes_df[dustyes_df['param'] == p[0]]['kl_shuf'].values[0] for p in PARAM_COMBOS]

    ax.bar(x - width/2, kl_real_yes, width, label='Real', color=COLORS['real'], edgecolor='black')
    ax.bar(x + width/2, kl_shuf_yes, width, label='Shuffled', color=COLORS['shuffled'], edgecolor='black')
    ax.set_xlabel('Parameter Set')
    ax.set_ylabel('KL Divergence from Uniform')
    ax.set_title('D. Position Non-Uniformity (DUST=yes)')
    ax.set_xticks(x)
    ax.set_xticklabels([p[0].replace('_', '\n') for p in PARAM_COMBOS], fontsize=8)
    ax.legend()

    fig.suptitle('Parameter Comparison Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig(output_dir / 'position_param_summary.png', dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: position_param_summary.png")

    # Save summary table
    summary_df.to_csv(output_dir / 'position_summary_stats.tsv', sep='\t', index=False, float_format='%.4f')
    print(f"  Saved: position_summary_stats.tsv")

    return summary_df


def main():
    print("=" * 70)
    print("POSITIONAL DISTRIBUTION ANALYSIS")
    print("Real vs Shuffled UTRs across BLAST Parameter Sets")
    print("=" * 70)

    results_dir = get_results_dir() / 'param_comparison_shuffle_v2'
    output_dir = Path('figures/param_position_analysis')

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        return 1

    # Load all data
    print("\nLoading BLAST results...")
    data = analyze_all_files(results_dir)

    # Generate visualizations
    print("\n" + "=" * 70)
    print("GENERATING VISUALIZATIONS")
    print("=" * 70)

    summary_df = plot_position_distributions(data, output_dir)

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    print("\nKey metrics (DUST=yes, optimal params go2_ge1_r1p-1):")
    opt_row = summary_df[(summary_df['param'] == 'go2_ge1_r1p-1') & (summary_df['dust'] == 'yes')].iloc[0]
    print(f"  Total real hits: {opt_row['n_real']:,}")
    print(f"  Total shuffled hits: {opt_row['n_shuf']:,}")
    print(f"  Real/Shuffled ratio: {opt_row['real_shuf_ratio']:.2f}x")
    print(f"  5' enrichment: {opt_row['five_prime_enrich']:.3f}x")
    print(f"  3' enrichment: {opt_row['three_prime_enrich']:.3f}x")
    print(f"  KL divergence (real): {opt_row['kl_real']:.4f}")
    print(f"  KL divergence (shuf): {opt_row['kl_shuf']:.4f}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print(f"Figures saved to: {output_dir}")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
