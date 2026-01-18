#!/usr/bin/env python3
"""
Generate ancient TE candidate figures (19-20) using v2 candidate data.

Creates publication-quality figures visualizing ancient TE fossil candidates.
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import re

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Colorblind-friendly palette
COLORS = {
    'primary': '#0077BB',     # Blue
    'secondary': '#EE7733',   # Orange
    'tertiary': '#009988',    # Teal
    'accent': '#CC3311',      # Red
    'neutral': '#BBBBBB',     # Gray
    'highlight': '#33BBEE',   # Cyan
}

# Directories
candidates_file = Path('results/repeatmasker_analysis_v2/ancient_te_candidates_clean.tsv')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("ANCIENT TE CANDIDATE FIGURES (19-20)")
print("=" * 70)

# Load ancient candidates
print("\nLoading ancient TE candidates...")

if not candidates_file.exists():
    print(f"  Error: {candidates_file} not found")
    sys.exit(1)

candidates = []
with open(candidates_file) as f:
    header = next(f).strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 11:
            try:
                candidates.append({
                    'rank': int(parts[0]),
                    'transcript': parts[1],
                    'chrom': parts[2],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'te_family': parts[5],
                    'te_id': parts[6],
                    'phyloP': float(parts[7]) if parts[7] not in ['NA', ''] else 0,
                    'pident': float(parts[8]) if parts[8] not in ['NA', ''] else 0,
                    'length': int(float(parts[9])) if parts[9] not in ['NA', ''] else 0,
                    'syntenic_species': int(parts[10]) if parts[10] not in ['NA', ''] else 0,
                })
            except (ValueError, IndexError):
                continue

print(f"  Loaded {len(candidates):,} ancient TE candidates")

# Extract metrics
phyloPs = np.array([c['phyloP'] for c in candidates])
pidents = np.array([c['pident'] for c in candidates])
lengths = np.array([c['length'] for c in candidates])
syntenic = np.array([c['syntenic_species'] for c in candidates])

print(f"  Conservation: median={np.median(phyloPs):.2f}, mean={np.mean(phyloPs):.2f}")
print(f"  Identity: median={np.median(pidents):.1f}%, mean={np.mean(pidents):.1f}%")
print(f"  Length: median={np.median(lengths):.0f}bp, mean={np.mean(lengths):.0f}bp")

# ============================================================
# FIGURE 19: Ancient Candidates Overview (4-panel)
# ============================================================
print("\nGenerating Figure 19: Ancient Candidates Overview...")

fig = plt.figure(figsize=(16, 14))

# Panel 1: Conservation vs Identity scatter
ax1 = fig.add_subplot(2, 2, 1)

# Sample for plotting if too many
if len(candidates) > 5000:
    idx = np.random.choice(len(candidates), 5000, replace=False)
    plot_phyloP = phyloPs[idx]
    plot_pident = pidents[idx]
else:
    plot_phyloP = phyloPs
    plot_pident = pidents

ax1.scatter(plot_pident, plot_phyloP, alpha=0.3, s=10, c=COLORS['primary'])
ax1.axhline(1, color=COLORS['accent'], linestyle='--', linewidth=1.5, label='High conservation')
ax1.axhline(2, color=COLORS['tertiary'], linestyle=':', linewidth=1.5, label='Very high conservation')
ax1.set_xlabel('Percent Identity (%)', fontsize=11)
ax1.set_ylabel('phyloP Conservation Score', fontsize=11)
ax1.set_title(f'Ancient TE Candidates: Conservation vs Identity\n(n={len(candidates):,})', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add marginal histograms
ax1_histx = ax1.inset_axes([0, 1.02, 1, 0.12])
ax1_histy = ax1.inset_axes([1.02, 0, 0.12, 1])
ax1_histx.hist(pidents, bins=30, color=COLORS['primary'], alpha=0.7, edgecolor='white')
ax1_histx.set_xlim(ax1.get_xlim())
ax1_histx.axis('off')
ax1_histy.hist(phyloPs, bins=30, orientation='horizontal', color=COLORS['primary'], alpha=0.7, edgecolor='white')
ax1_histy.set_ylim(ax1.get_ylim())
ax1_histy.axis('off')

# Panel 2: Identity distribution comparison
ax2 = fig.add_subplot(2, 2, 2)

# Compare ancient candidates to typical range
bins = np.linspace(60, 100, 30)
ax2.hist(pidents, bins=bins, alpha=0.8, color=COLORS['tertiary'],
         edgecolor='white', label=f'Ancient candidates (n={len(candidates):,})')
ax2.axvline(np.median(pidents), color=COLORS['accent'], linestyle='--',
            linewidth=2, label=f'Median: {np.median(pidents):.1f}%')
ax2.set_xlabel('Percent Identity (%)', fontsize=11)
ax2.set_ylabel('Number of Candidates', fontsize=11)
ax2.set_title('Identity Distribution of Ancient Candidates', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Panel 3: TE Family distribution (top 15)
ax3 = fig.add_subplot(2, 2, 3)

te_family_counts = defaultdict(int)
for c in candidates:
    te_family_counts[c['te_family']] += 1

top_families = sorted(te_family_counts.items(), key=lambda x: -x[1])[:15]
families = [f[0] for f in top_families]
counts = [f[1] for f in top_families]

colors_fam = plt.cm.tab20(np.linspace(0, 1, len(families)))
bars = ax3.barh(range(len(families)), counts, color=colors_fam, edgecolor='black', linewidth=0.5)
ax3.set_yticks(range(len(families)))
ax3.set_yticklabels(families, fontsize=10)
ax3.invert_yaxis()
ax3.set_xlabel('Number of Ancient Candidates', fontsize=11)
ax3.set_title('Top 15 TE Families in Ancient Candidates', fontsize=12, fontweight='bold')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

# Add count labels
for bar, count in zip(bars, counts):
    ax3.text(count + max(counts)*0.02, bar.get_y() + bar.get_height()/2,
             f'{count:,}', va='center', fontsize=9)

# Panel 4: Synteny species distribution
ax4 = fig.add_subplot(2, 2, 4)

syn_counts = defaultdict(int)
for c in candidates:
    syn_counts[c['syntenic_species']] += 1

species_nums = sorted(syn_counts.keys())
species_counts = [syn_counts[s] for s in species_nums]

colors_syn = [COLORS['neutral'], COLORS['secondary'], COLORS['tertiary'], COLORS['primary']][:len(species_nums)]
bars = ax4.bar(species_nums, species_counts, color=colors_syn, edgecolor='black', linewidth=1)
ax4.set_xlabel('Number of Species with Synteny', fontsize=11)
ax4.set_ylabel('Number of Candidates', fontsize=11)
ax4.set_title('Synteny Conservation Across Species', fontsize=12, fontweight='bold')
ax4.set_xticks(species_nums)
ax4.set_xticklabels([f'{s} species' for s in species_nums], fontsize=10)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)

# Add percentage labels
total = sum(species_counts)
for bar, count in zip(bars, species_counts):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height + total*0.01,
             f'{count:,}\n({100*count/total:.1f}%)',
             ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / '19_ancient_candidates_overview.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/19_ancient_candidates_overview.png")
plt.close()

# ============================================================
# FIGURE 20: Ancient Candidates Details (4-panel)
# ============================================================
print("\nGenerating Figure 20: Ancient Candidates Details...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: Conservation histogram
ax = axes[0, 0]
ax.hist(phyloPs, bins=50, color=COLORS['tertiary'], edgecolor='white', alpha=0.8)
ax.axvline(np.median(phyloPs), color=COLORS['accent'], linestyle='--', linewidth=2,
           label=f'Median: {np.median(phyloPs):.2f}')
ax.axvline(np.mean(phyloPs), color=COLORS['primary'], linestyle=':', linewidth=2,
           label=f'Mean: {np.mean(phyloPs):.2f}')
ax.axvline(1, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('phyloP Conservation Score', fontsize=11)
ax.set_ylabel('Number of Candidates', fontsize=11)
ax.set_title('A. Conservation Distribution', fontsize=12, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Length histogram
ax = axes[0, 1]
ax.hist(lengths, bins=50, color=COLORS['primary'], edgecolor='white', alpha=0.8)
ax.axvline(np.median(lengths), color=COLORS['accent'], linestyle='--', linewidth=2,
           label=f'Median: {np.median(lengths):.0f}bp')
ax.set_xlabel('Alignment Length (bp)', fontsize=11)
ax.set_ylabel('Number of Candidates', fontsize=11)
ax.set_title('B. Length Distribution', fontsize=12, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Conservation by length bins
ax = axes[1, 0]
len_bins = [(0, 30), (30, 50), (50, 100), (100, 200), (200, 500)]
bin_labels = ['<30bp', '30-50bp', '50-100bp', '100-200bp', '>200bp']
bin_phyloPs = []

for lo, hi in len_bins:
    subset = [c['phyloP'] for c in candidates if lo <= c['length'] < hi]
    bin_phyloPs.append(subset if subset else [0])

bp = ax.boxplot(bin_phyloPs, labels=bin_labels, patch_artist=True)
colors_bp = [COLORS['neutral'], COLORS['secondary'], COLORS['tertiary'], COLORS['primary'], COLORS['highlight']]
for patch, color in zip(bp['boxes'], colors_bp):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.axhline(1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel('Alignment Length', fontsize=11)
ax.set_ylabel('phyloP Conservation Score', fontsize=11)
ax.set_title('C. Conservation by Alignment Length', fontsize=12, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Identity vs Length scatter
ax = axes[1, 1]

# Color by conservation
colors = []
for c in candidates:
    if c['phyloP'] >= 2:
        colors.append(COLORS['tertiary'])
    elif c['phyloP'] >= 1:
        colors.append(COLORS['primary'])
    else:
        colors.append(COLORS['neutral'])

# Sample for plotting
if len(candidates) > 3000:
    idx = np.random.choice(len(candidates), 3000, replace=False)
    plot_pidents = pidents[idx]
    plot_lengths = lengths[idx]
    plot_colors = [colors[i] for i in idx]
else:
    plot_pidents = pidents
    plot_lengths = lengths
    plot_colors = colors

ax.scatter(plot_pidents, plot_lengths, c=plot_colors, alpha=0.5, s=20, edgecolors='black', linewidths=0.3)
ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Alignment Length (bp)', fontsize=11)
ax.set_title('D. Identity vs Length\n(color = conservation)', fontsize=12, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLORS['tertiary'], edgecolor='black', label='phyloP >= 2'),
    Patch(facecolor=COLORS['primary'], edgecolor='black', label='phyloP 1-2'),
    Patch(facecolor=COLORS['neutral'], edgecolor='black', label='phyloP < 1'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / '20_ancient_candidates_details.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/20_ancient_candidates_details.png")
plt.close()

# Summary statistics
print("\n" + "=" * 70)
print("ANCIENT TE CANDIDATE STATISTICS")
print("=" * 70)
print(f"\nTotal candidates: {len(candidates):,}")
print(f"\nConservation (phyloP):")
print(f"  Mean: {np.mean(phyloPs):.2f}")
print(f"  Median: {np.median(phyloPs):.2f}")
print(f"  Range: {np.min(phyloPs):.2f} - {np.max(phyloPs):.2f}")
print(f"\nAlignment length:")
print(f"  Mean: {np.mean(lengths):.0f}bp")
print(f"  Median: {np.median(lengths):.0f}bp")
print(f"  Range: {np.min(lengths)} - {np.max(lengths)}bp")
print(f"\nIdentity:")
print(f"  Mean: {np.mean(pidents):.1f}%")
print(f"  Median: {np.median(pidents):.1f}%")
print(f"\nSyntenic species:")
for n in sorted(syn_counts.keys()):
    print(f"  {n} species: {syn_counts[n]:,} ({100*syn_counts[n]/len(candidates):.1f}%)")

print("\n" + "=" * 70)
print("ANCIENT CANDIDATE FIGURES COMPLETE")
print("=" * 70)
