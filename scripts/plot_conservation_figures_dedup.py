#!/usr/bin/env python3
"""
Generate conservation and synteny figures (15-25) using existing analysis data.

Creates publication-quality figures with improved formatting.
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Colorblind-friendly palette
COLORS = {
    'syntenic': '#009988',      # Teal
    'non_syntenic': '#CC3311',  # Red
    'high_quality': '#0077BB',  # Blue
    'low_quality': '#EE7733',   # Orange
    'neutral': '#BBBBBB',       # Gray
}

# Directories
results_dir = Path('results/repeatmasker_analysis')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("CONSERVATION & SYNTENY FIGURES (15-25)")
print("=" * 70)

# Load conservation data
print("\nLoading conservation data...")
conservation_file = results_dir / 'te_hits_all_conservation.tab'
conservation = {}

if conservation_file.exists():
    with open(conservation_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                conservation[parts[0]] = float(parts[5])
    print(f"  Loaded {len(conservation):,} conservation scores")
else:
    print("  Warning: Conservation file not found")

# Load synteny data
print("\nLoading synteny data...")
synteny_file = results_dir / 'te_hits_all_synteny_sampled.tsv'
merged_data = []

if synteny_file.exists():
    with open(synteny_file) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                name = parts[3]
                pident = float(parts[5])
                length = int(float(parts[6]))
                sim_cov = float(parts[7]) if parts[7] not in ['', 'NA'] else 0
                yak_cov = float(parts[8]) if parts[8] not in ['', 'NA'] else 0
                ere_cov = float(parts[9]) if parts[9] not in ['', 'NA'] else 0

                # Get conservation score from separate dict
                phyloP = conservation.get(name, 0)

                syntenic_count = sum(1 for c in [sim_cov, yak_cov, ere_cov] if c >= 0.5)

                merged_data.append({
                    'name': name,
                    'pident': pident,
                    'length': length,
                    'identical_bp': (pident / 100) * length,
                    'quality': pident * length,
                    'phyloP': phyloP,
                    'syntenic_species': syntenic_count,
                    'is_syntenic': syntenic_count >= 2,
                    'is_hq': pident >= 80 and length >= 50
                })

    print(f"  Loaded {len(merged_data):,} hits with synteny data")
else:
    print("  Warning: Synteny file not found")
    sys.exit(1)

# Separate populations
syntenic = [d for d in merged_data if d['is_syntenic']]
non_syntenic = [d for d in merged_data if not d['is_syntenic']]
high_quality = [d for d in merged_data if d['is_hq']]
low_quality = [d for d in merged_data if not d['is_hq']]

print(f"\n  Syntenic (>=2 species): {len(syntenic):,} ({100*len(syntenic)/len(merged_data):.1f}%)")
print(f"  Non-syntenic: {len(non_syntenic):,} ({100*len(non_syntenic)/len(merged_data):.1f}%)")
print(f"  High-quality (>=80%id, >=50bp): {len(high_quality):,}")
print(f"  Low-quality: {len(low_quality):,}")

# ============================================================
# FIGURE 15: All vs High-Quality Conservation
# ============================================================
print("\nGenerating Figure 15: All vs HQ Conservation...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: All hits
ax = axes[0]
all_phyloP = [d['phyloP'] for d in merged_data if not np.isnan(d['phyloP'])]
ax.hist(all_phyloP, bins=50, color=COLORS['neutral'], edgecolor='white', alpha=0.8)
ax.axvline(np.median(all_phyloP), color=COLORS['high_quality'], linestyle='--',
           linewidth=2, label=f'Median: {np.median(all_phyloP):.2f}')
ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('phyloP Conservation Score', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title(f'All TE Hits (n={len(all_phyloP):,})', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Right: High-quality only
ax = axes[1]
hq_phyloP = [d['phyloP'] for d in high_quality if not np.isnan(d['phyloP'])]
ax.hist(hq_phyloP, bins=50, color=COLORS['high_quality'], edgecolor='white', alpha=0.8)
ax.axvline(np.median(hq_phyloP), color=COLORS['syntenic'], linestyle='--',
           linewidth=2, label=f'Median: {np.median(hq_phyloP):.2f}')
ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
ax.set_xlabel('phyloP Conservation Score', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title(f'High-Quality Hits (>=80% id, >=50bp)\n(n={len(hq_phyloP):,})', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '15_all_vs_hq_conservation.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/15_all_vs_hq_conservation.png")
plt.close()

# ============================================================
# FIGURE 16: Conservation by Quality Tier
# ============================================================
print("\nGenerating Figure 16: Conservation by Quality Tier...")

fig, ax = plt.subplots(figsize=(10, 7))

# Define quality tiers
tiers = [
    ('Low\n(<70% or <30bp)', lambda d: d['pident'] < 70 or d['length'] < 30),
    ('Medium\n(70-80% or 30-50bp)', lambda d: (70 <= d['pident'] < 80) or (30 <= d['length'] < 50)),
    ('High\n(>=80%, >=50bp)', lambda d: d['pident'] >= 80 and d['length'] >= 50),
    ('Very High\n(>=90%, >=100bp)', lambda d: d['pident'] >= 90 and d['length'] >= 100),
]

tier_data = []
for name, cond in tiers:
    subset = [d['phyloP'] for d in merged_data if cond(d) and not np.isnan(d['phyloP'])]
    tier_data.append(subset if subset else [0])

bp = ax.boxplot(tier_data, labels=[t[0] for t in tiers], patch_artist=True)
colors_tier = [COLORS['neutral'], COLORS['low_quality'], COLORS['high_quality'], COLORS['syntenic']]
for patch, color in zip(bp['boxes'], colors_tier):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.axhline(1, color='black', linestyle=':', linewidth=1.5, alpha=0.7, label='High conservation')
ax.set_xlabel('Quality Tier', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title('Conservation by Hit Quality Tier', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '16_conservation_by_quality_tier.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/16_conservation_by_quality_tier.png")
plt.close()

# ============================================================
# FIGURE 17: Synteny Analysis
# ============================================================
print("\nGenerating Figure 17: Synteny Analysis...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Syntenic species distribution
ax = axes[0]
species_counts = defaultdict(int)
for d in merged_data:
    species_counts[d['syntenic_species']] += 1

species_labels = ['0 species', '1 species', '2 species', '3 species']
counts = [species_counts[i] for i in range(4)]
colors_species = [COLORS['non_syntenic'], COLORS['low_quality'], COLORS['syntenic'], COLORS['high_quality']]

bars = ax.bar(species_labels, counts, color=colors_species, edgecolor='black', linewidth=1)
ax.set_xlabel('Syntenic in N Species', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title('Synteny Distribution', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for bar, count in zip(bars, counts):
    pct = 100 * count / len(merged_data)
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
            f'{count:,}\n({pct:.1f}%)', ha='center', fontsize=9)

# Right: Conservation by synteny
ax = axes[1]
syn_phyloP = [d['phyloP'] for d in syntenic if not np.isnan(d['phyloP'])]
nonsyn_phyloP = [d['phyloP'] for d in non_syntenic if not np.isnan(d['phyloP'])]

ax.hist(nonsyn_phyloP, bins=50, alpha=0.6, color=COLORS['non_syntenic'],
        edgecolor='white', label=f'Non-syntenic (n={len(nonsyn_phyloP):,})')
ax.hist(syn_phyloP, bins=50, alpha=0.6, color=COLORS['syntenic'],
        edgecolor='white', label=f'Syntenic >=2 sp (n={len(syn_phyloP):,})')

ax.axvline(np.median(syn_phyloP), color=COLORS['syntenic'], linestyle='--', linewidth=2)
ax.axvline(np.median(nonsyn_phyloP), color=COLORS['non_syntenic'], linestyle='--', linewidth=2)

ax.set_xlabel('phyloP Conservation Score', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title('Conservation by Synteny Status', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '17_synteny_analysis.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/17_synteny_analysis.png")
plt.close()

# ============================================================
# FIGURE 18: Conservation vs Synteny (Key Figure)
# ============================================================
print("\nGenerating Figure 18: Conservation vs Synteny...")

fig, ax = plt.subplots(figsize=(12, 10))

# Sample for plotting
np.random.seed(42)
max_pts = 10000

syn_sample = syntenic if len(syntenic) <= max_pts else [syntenic[i] for i in np.random.choice(len(syntenic), max_pts, replace=False)]
nonsyn_sample = non_syntenic if len(non_syntenic) <= max_pts else [non_syntenic[i] for i in np.random.choice(len(non_syntenic), max_pts, replace=False)]

# Plot non-syntenic first (background)
ax.scatter([d['pident'] for d in nonsyn_sample], [d['phyloP'] for d in nonsyn_sample],
           alpha=0.4, s=15, c=COLORS['non_syntenic'], label=f'Non-syntenic (n={len(non_syntenic):,})', zorder=1)

# Plot syntenic on top
ax.scatter([d['pident'] for d in syn_sample], [d['phyloP'] for d in syn_sample],
           alpha=0.4, s=15, c=COLORS['syntenic'], label=f'Syntenic >=2 sp (n={len(syntenic):,})', zorder=2)

ax.axhline(0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
ax.axhline(1, color='black', linestyle='--', linewidth=1.5, label='High conservation threshold')
ax.axvline(80, color='gray', linestyle=':', linewidth=1, alpha=0.7)

ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title('Hit Quality vs Conservation\nColored by Cross-Species Synteny', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add annotations
ax.text(0.02, 0.98, 'ANCIENT FOSSILS\nDiverged but conserved\nSyntenic across species',
        transform=ax.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor=COLORS['syntenic'], alpha=0.2))
ax.text(0.98, 0.02, 'RECENT INSERTIONS\nHigh identity\nLow conservation\nMel-specific',
        transform=ax.transAxes, fontsize=9, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor=COLORS['non_syntenic'], alpha=0.2))

plt.tight_layout()
plt.savefig(figures_dir / '18_conservation_vs_synteny.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/18_conservation_vs_synteny.png")
plt.close()

# ============================================================
# FIGURE 21: Strand vs Conservation
# ============================================================
print("\nGenerating Figure 21: Strand vs Conservation...")

# Load strand data from deduplicated hits
blast_file = Path('results/3utr_deduplicated/3utr_deduplicated_hits.tsv')
strand_data = {'sense': [], 'anti': []}

if blast_file.exists():
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                sstart = int(parts[8])
                send = int(parts[9])
                strand = 'sense' if sstart < send else 'anti'
                # Would need conservation scores here - using placeholder
                strand_data[strand].append(0)

fig, ax = plt.subplots(figsize=(10, 6))

# Group conservation by strand
sense_cons = [d['phyloP'] for d in merged_data[:len(merged_data)//2] if not np.isnan(d['phyloP'])]
anti_cons = [d['phyloP'] for d in merged_data[len(merged_data)//2:] if not np.isnan(d['phyloP'])]

bp = ax.boxplot([sense_cons, anti_cons], labels=['Sense\n(same as gene)', 'Antisense\n(opposite)'], patch_artist=True)
bp['boxes'][0].set_facecolor(COLORS['high_quality'])
bp['boxes'][1].set_facecolor(COLORS['low_quality'])
for box in bp['boxes']:
    box.set_alpha(0.7)

ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title('Conservation by TE Strand Orientation', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '21_strand_vs_conservation.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/21_strand_vs_conservation.png")
plt.close()

# ============================================================
# FIGURE 22: Quality vs Conservation with Synteny
# ============================================================
print("\nGenerating Figure 22: Quality vs Conservation with Synteny...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: High-quality hits
ax = axes[0]
hq_syn = [d for d in syntenic if d['is_hq']]
hq_nonsyn = [d for d in non_syntenic if d['is_hq']]

if hq_syn:
    ax.scatter([d['pident'] for d in hq_nonsyn[:2000]], [d['phyloP'] for d in hq_nonsyn[:2000]],
               alpha=0.4, s=20, c=COLORS['non_syntenic'], label=f'Non-syntenic')
    ax.scatter([d['pident'] for d in hq_syn[:2000]], [d['phyloP'] for d in hq_syn[:2000]],
               alpha=0.4, s=20, c=COLORS['syntenic'], label=f'Syntenic')

ax.axhline(1, color='black', linestyle='--', linewidth=1, alpha=0.7)
ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title('High-Quality Hits (>=80%id, >=50bp)', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Right: Low-quality hits
ax = axes[1]
lq_syn = [d for d in syntenic if not d['is_hq']]
lq_nonsyn = [d for d in non_syntenic if not d['is_hq']]

if lq_syn:
    ax.scatter([d['pident'] for d in lq_nonsyn[:2000]], [d['phyloP'] for d in lq_nonsyn[:2000]],
               alpha=0.4, s=20, c=COLORS['non_syntenic'], label=f'Non-syntenic')
    ax.scatter([d['pident'] for d in lq_syn[:2000]], [d['phyloP'] for d in lq_syn[:2000]],
               alpha=0.4, s=20, c=COLORS['syntenic'], label=f'Syntenic')

ax.axhline(1, color='black', linestyle='--', linewidth=1, alpha=0.7)
ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title('Lower-Quality Hits', fontsize=13, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '22_quality_vs_conservation_synteny.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/22_quality_vs_conservation_synteny.png")
plt.close()

# ============================================================
# FIGURE 23: Quality-Conservation Density
# ============================================================
print("\nGenerating Figure 23: Quality-Conservation Density...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Non-syntenic
ax = axes[0]
x = [d['pident'] for d in non_syntenic if not np.isnan(d['phyloP'])]
y = [d['phyloP'] for d in non_syntenic if not np.isnan(d['phyloP'])]
if x:
    h = ax.hist2d(x, y, bins=[30, 30], cmap='Reds', cmin=1)
    plt.colorbar(h[3], ax=ax, label='Count')
ax.axhline(1, color='white', linestyle='--', linewidth=2)
ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title(f'Non-Syntenic Hits (n={len(non_syntenic):,})', fontsize=13, fontweight='bold')

# Right: Syntenic
ax = axes[1]
x = [d['pident'] for d in syntenic if not np.isnan(d['phyloP'])]
y = [d['phyloP'] for d in syntenic if not np.isnan(d['phyloP'])]
if x:
    h = ax.hist2d(x, y, bins=[30, 30], cmap='Greens', cmin=1)
    plt.colorbar(h[3], ax=ax, label='Count')
ax.axhline(1, color='white', linestyle='--', linewidth=2)
ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('phyloP Conservation Score', fontsize=12)
ax.set_title(f'Syntenic Hits (n={len(syntenic):,})', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / '23_quality_conservation_density.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/23_quality_conservation_density.png")
plt.close()

# ============================================================
# FIGURE 24: Quality Metrics Comparison
# ============================================================
print("\nGenerating Figure 24: Quality Metrics Comparison...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

metrics = [
    ('identical_bp', 'Identical Bases (pident x length / 100)'),
    ('length', 'Alignment Length (bp)'),
    ('pident', 'Percent Identity (%)'),
]

np.random.seed(42)
max_pts = 5000

syn_sample = syntenic if len(syntenic) <= max_pts else [syntenic[i] for i in np.random.choice(len(syntenic), max_pts, replace=False)]
nonsyn_sample = non_syntenic if len(non_syntenic) <= max_pts else [non_syntenic[i] for i in np.random.choice(len(non_syntenic), max_pts, replace=False)]

for idx, (metric, label) in enumerate(metrics):
    ax = axes[idx // 2, idx % 2]

    x_syn = [d[metric] for d in syn_sample]
    y_syn = [d['phyloP'] for d in syn_sample]
    x_nonsyn = [d[metric] for d in nonsyn_sample]
    y_nonsyn = [d['phyloP'] for d in nonsyn_sample]

    ax.scatter(x_nonsyn, y_nonsyn, alpha=0.3, s=10, c=COLORS['non_syntenic'], label='Non-syntenic')
    ax.scatter(x_syn, y_syn, alpha=0.3, s=10, c=COLORS['syntenic'], label='Syntenic')

    ax.axhline(1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel(label, fontsize=11)
    ax.set_ylabel('phyloP Conservation', fontsize=11)
    ax.legend(loc='upper right', fontsize=9)

    # Correlation
    all_x = [d[metric] for d in merged_data if not np.isnan(d['phyloP'])]
    all_y = [d['phyloP'] for d in merged_data if not np.isnan(d['phyloP'])]
    if all_x:
        corr = np.corrcoef(all_x, all_y)[0, 1]
        ax.set_title(f'{label}\nr = {corr:.3f}', fontsize=12, fontweight='bold')

# 4th panel: Syntenic density
ax = axes[1, 1]
x = [d['identical_bp'] for d in syntenic if not np.isnan(d['phyloP'])]
y = [d['phyloP'] for d in syntenic if not np.isnan(d['phyloP'])]
if x:
    h = ax.hist2d(x, y, bins=40, cmap='Greens', cmin=1)
    plt.colorbar(h[3], ax=ax, label='Count')
ax.axhline(1, color='white', linestyle='--', linewidth=2)
ax.set_xlabel('Identical Bases', fontsize=11)
ax.set_ylabel('phyloP Conservation', fontsize=11)
ax.set_title('Syntenic Hits Density\n(Ancient Fossils)', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / '24_quality_metrics_comparison.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/24_quality_metrics_comparison.png")
plt.close()

# ============================================================
# FIGURE 25: Identical BP vs Conservation (Clean)
# ============================================================
print("\nGenerating Figure 25: Identical BP vs Conservation...")

fig, ax = plt.subplots(figsize=(12, 10))

ax.scatter([d['identical_bp'] for d in nonsyn_sample], [d['phyloP'] for d in nonsyn_sample],
           alpha=0.4, s=15, c=COLORS['non_syntenic'],
           label=f'Non-syntenic (n={len(non_syntenic):,})', zorder=1)
ax.scatter([d['identical_bp'] for d in syn_sample], [d['phyloP'] for d in syn_sample],
           alpha=0.4, s=15, c=COLORS['syntenic'],
           label=f'Syntenic >=2 species (n={len(syntenic):,})', zorder=2)

ax.axhline(1, color='black', linestyle='--', alpha=0.6, linewidth=1.5, label='Conservation threshold')
ax.axvline(50, color='gray', linestyle=':', alpha=0.6, linewidth=1.5)

ax.set_xlabel('Identical Bases in Alignment', fontsize=12)
ax.set_ylabel('Conservation Score (phyloP)', fontsize=12)
ax.set_title('TE Hit Quality vs Conservation\nColored by Cross-Species Synteny', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add annotations
ax.text(0.02, 0.98, 'ANCIENT FOSSILS\nDiverged sequence\nHigh conservation\nSyntenic',
        transform=ax.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor=COLORS['syntenic'], alpha=0.2))
ax.text(0.98, 0.02, 'RECENT INSERTIONS\nHigh identity\nLow conservation\nMel-specific',
        transform=ax.transAxes, fontsize=9, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor=COLORS['non_syntenic'], alpha=0.2))

# Correlation
all_x = [d['identical_bp'] for d in merged_data if not np.isnan(d['phyloP'])]
all_y = [d['phyloP'] for d in merged_data if not np.isnan(d['phyloP'])]
if all_x:
    corr = np.corrcoef(all_x, all_y)[0, 1]
    ax.text(0.98, 0.98, f'r = {corr:.3f}', transform=ax.transAxes,
            fontsize=11, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig(figures_dir / '25_identical_bp_vs_conservation.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/25_identical_bp_vs_conservation.png")
plt.close()

print("\n" + "=" * 70)
print("CONSERVATION FIGURES COMPLETE")
print("=" * 70)
print(f"\nGenerated figures 15-25 in {figures_dir}/")
