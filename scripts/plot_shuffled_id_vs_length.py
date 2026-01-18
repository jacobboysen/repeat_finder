#!/usr/bin/env python3
"""
Generate identity vs length 2D histograms for real and 10 shuffled datasets.

Creates an 11-panel figure comparing the distribution of hits across
identity and length dimensions.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Directories
shuffled_dir = Path('results/shuffled_controls/deduplicated')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("IDENTITY vs LENGTH: Real + 10 Shuffled Controls")
print("=" * 70)

def load_blast_hits(filepath):
    """Load pident and length from BLAST TSV file."""
    pidents = []
    lengths = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                try:
                    pidents.append(float(parts[2]))
                    lengths.append(int(parts[3]))
                except (ValueError, IndexError):
                    continue
    return np.array(pidents), np.array(lengths)

# Load all datasets
datasets = []

# Real data first
print("\nLoading datasets...")
real_file = shuffled_dir / 'real_deduplicated_hits.tsv'
if real_file.exists():
    pidents, lengths = load_blast_hits(real_file)
    datasets.append(('Real (Deduplicated)', pidents, lengths))
    print(f"  Real: {len(pidents):,} hits")
else:
    print(f"  Error: {real_file} not found")
    sys.exit(1)

# 10 shuffled replicates
for i in range(1, 11):
    shuf_file = shuffled_dir / f'shuffled_rep{i}_deduplicated_hits.tsv'
    if shuf_file.exists():
        pidents, lengths = load_blast_hits(shuf_file)
        datasets.append((f'Shuffled Rep {i}', pidents, lengths))
        print(f"  Shuffled {i}: {len(pidents):,} hits")
    else:
        print(f"  Warning: {shuf_file} not found")

print(f"\nTotal datasets: {len(datasets)}")

# Create 11-panel figure (3 rows x 4 cols, with one empty)
fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.flatten()

# Common parameters for 2D histogram
bins_pident = np.linspace(60, 100, 41)
bins_length = np.linspace(0, 300, 31)

# Find global max for consistent color scale
all_hists = []
for name, pidents, lengths in datasets:
    h, _, _ = np.histogram2d(pidents, np.minimum(lengths, 300), bins=[bins_pident, bins_length])
    all_hists.append(h.max())
vmax = max(all_hists)

print(f"\nGenerating figure (vmax={vmax:.0f})...")

for idx, (name, pidents, lengths) in enumerate(datasets):
    ax = axes[idx]

    # Cap lengths at 300 for visualization
    lengths_capped = np.minimum(lengths, 300)

    # 2D histogram
    h = ax.hist2d(pidents, lengths_capped, bins=[bins_pident, bins_length],
                  cmap='YlOrRd', vmin=1, vmax=vmax)

    # Reference lines
    ax.axhline(50, color='black', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(80, color='black', linestyle=':', linewidth=1, alpha=0.7)

    # Labels
    ax.set_xlabel('Identity (%)', fontsize=9)
    ax.set_ylabel('Length (bp)', fontsize=9)

    # Title with hit count
    if 'Real' in name:
        ax.set_title(f'{name}\n(n={len(pidents):,})', fontsize=10, fontweight='bold', color='darkblue')
    else:
        ax.set_title(f'{name}\n(n={len(pidents):,})', fontsize=10)

# Hide the last subplot (12th position) since we only have 11 datasets
axes[11].axis('off')

# Add a single colorbar
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
sm = plt.cm.ScalarMappable(cmap='YlOrRd', norm=plt.Normalize(vmin=1, vmax=vmax))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Number of Hits', fontsize=11)

# Main title
fig.suptitle('Identity vs Length Distribution: Real vs Shuffled Controls\n(Deduplicated Data)',
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 0.91, 0.95])
plt.savefig(figures_dir / '30_real_vs_shuffled_id_length.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"\nSaved: {figures_dir}/30_real_vs_shuffled_id_length.png")
plt.close()

# Also create summary statistics table
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)
print(f"\n{'Dataset':<20} {'Hits':>10} {'Med ID':>8} {'Med Len':>8} {'HQ Hits':>10} {'HQ %':>8}")
print("-" * 70)

for name, pidents, lengths in datasets:
    hq_count = np.sum((pidents >= 80) & (lengths >= 50))
    hq_pct = 100 * hq_count / len(pidents) if len(pidents) > 0 else 0
    print(f"{name:<20} {len(pidents):>10,} {np.median(pidents):>8.1f} {np.median(lengths):>8.0f} {hq_count:>10,} {hq_pct:>7.1f}%")

print("\n" + "=" * 70)
print("FIGURE COMPLETE")
print("=" * 70)
