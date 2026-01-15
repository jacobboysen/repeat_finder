#!/usr/bin/env python3
"""
Create the killer figure: BLAST quality vs conservation, colored by synteny.

Hypothesis: Two distinct populations:
1. Recent insertions: high quality, low conservation, non-syntenic
2. Ancient fossils: low quality, high conservation, syntenic
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def load_conservation(filepath):
    """Load conservation scores keyed by hit name."""
    data = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            phyloP = float(parts[5])
            data[name] = phyloP
    return data

def load_synteny(filepath):
    """Load synteny data keyed by hit name."""
    data = {}
    with open(filepath) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            name = parts[3]
            pident = float(parts[5])
            length = float(parts[6])
            # Count syntenic species (coverage >= 0.5)
            sim_cov = float(parts[8])
            yak_cov = float(parts[9])
            ere_cov = float(parts[10])
            any_cov = float(parts[11])

            syntenic_count = sum(1 for c in [sim_cov, yak_cov, ere_cov] if c >= 0.5)

            data[name] = {
                'pident': pident,
                'length': length,
                'quality': pident * length,  # BLAST quality metric
                'syntenic_species': syntenic_count,
                'is_syntenic': syntenic_count >= 2  # syntenic in at least 2 species
            }
    return data

print("Loading data...")
conservation = load_conservation('results/repeatmasker_analysis/te_hits_all_conservation.tab')
synteny = load_synteny('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')

# Merge datasets
merged = []
for name, syn_data in synteny.items():
    if name in conservation:
        merged.append({
            'name': name,
            'quality': syn_data['quality'],
            'pident': syn_data['pident'],
            'length': syn_data['length'],
            'phyloP': conservation[name],
            'syntenic_species': syn_data['syntenic_species'],
            'is_syntenic': syn_data['is_syntenic']
        })

print(f"Merged {len(merged):,} hits with both conservation and synteny data")

# Separate by synteny status
syntenic = [d for d in merged if d['is_syntenic']]
non_syntenic = [d for d in merged if not d['is_syntenic']]

print(f"  Syntenic (≥2 species): {len(syntenic):,} ({100*len(syntenic)/len(merged):.1f}%)")
print(f"  Non-syntenic: {len(non_syntenic):,} ({100*len(non_syntenic)/len(merged):.1f}%)")

# Create figure
fig = plt.figure(figsize=(14, 12))

# Main scatter plot
ax_main = fig.add_axes([0.1, 0.1, 0.65, 0.65])

# Sample if too many points
max_points = 10000
if len(syntenic) > max_points:
    idx = np.random.choice(len(syntenic), max_points, replace=False)
    syn_sample = [syntenic[i] for i in idx]
else:
    syn_sample = syntenic

if len(non_syntenic) > max_points:
    idx = np.random.choice(len(non_syntenic), max_points, replace=False)
    nonsyn_sample = [non_syntenic[i] for i in idx]
else:
    nonsyn_sample = non_syntenic

# Plot non-syntenic first (background)
x_nonsyn = [d['quality'] for d in nonsyn_sample]
y_nonsyn = [d['phyloP'] for d in nonsyn_sample]
ax_main.scatter(x_nonsyn, y_nonsyn, alpha=0.3, s=8, c='#e74c3c', label=f'Non-syntenic (n={len(non_syntenic):,})')

# Plot syntenic on top
x_syn = [d['quality'] for d in syn_sample]
y_syn = [d['phyloP'] for d in syn_sample]
ax_main.scatter(x_syn, y_syn, alpha=0.3, s=8, c='#27ae60', label=f'Syntenic ≥2 species (n={len(syntenic):,})')

# Add quadrant lines
ax_main.axhline(1, color='black', linestyle='--', alpha=0.5, linewidth=1.5)
ax_main.axvline(5000, color='black', linestyle='--', alpha=0.5, linewidth=1.5)

# Labels
ax_main.set_xlabel('BLAST Quality (Identity % × Length)', fontsize=12)
ax_main.set_ylabel('Conservation Score (phyloP)', fontsize=12)
ax_main.set_title('TE Hit Quality vs Conservation by Synteny Status', fontsize=14, fontweight='bold')
ax_main.legend(loc='upper right', fontsize=10)

# Add quadrant labels
ax_main.text(0.02, 0.98, 'ANCIENT FOSSILS\n(low quality, high conservation)',
             transform=ax_main.transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='#27ae60', alpha=0.3))
ax_main.text(0.98, 0.02, 'RECENT INSERTIONS\n(high quality, low conservation)',
             transform=ax_main.transAxes, fontsize=10, verticalalignment='bottom',
             horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='#e74c3c', alpha=0.3))

# Top histogram (quality distribution)
ax_top = fig.add_axes([0.1, 0.77, 0.65, 0.15])
bins_qual = np.linspace(0, max(x_syn + x_nonsyn), 50)
ax_top.hist([d['quality'] for d in non_syntenic], bins=bins_qual, alpha=0.5,
            color='#e74c3c', density=True, label='Non-syntenic')
ax_top.hist([d['quality'] for d in syntenic], bins=bins_qual, alpha=0.5,
            color='#27ae60', density=True, label='Syntenic')
ax_top.set_xlim(ax_main.get_xlim())
ax_top.set_ylabel('Density')
ax_top.set_xticklabels([])
ax_top.legend(fontsize=8)

# Right histogram (conservation distribution)
ax_right = fig.add_axes([0.77, 0.1, 0.15, 0.65])
bins_cons = np.linspace(min(y_syn + y_nonsyn), max(y_syn + y_nonsyn), 50)
ax_right.hist([d['phyloP'] for d in non_syntenic], bins=bins_cons, alpha=0.5,
              color='#e74c3c', density=True, orientation='horizontal')
ax_right.hist([d['phyloP'] for d in syntenic], bins=bins_cons, alpha=0.5,
              color='#27ae60', density=True, orientation='horizontal')
ax_right.set_ylim(ax_main.get_ylim())
ax_right.set_xlabel('Density')
ax_right.set_yticklabels([])

plt.savefig('figures/repeatmasker_comparison/22_quality_vs_conservation_synteny.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/repeatmasker_comparison/22_quality_vs_conservation_synteny.png")

# Quantify the populations
print("\n" + "="*70)
print("POPULATION ANALYSIS")
print("="*70)

# Define quadrants
ancient_fossils = [d for d in merged if d['quality'] < 5000 and d['phyloP'] > 1 and d['is_syntenic']]
recent_insertions = [d for d in merged if d['quality'] >= 5000 and d['phyloP'] <= 1 and not d['is_syntenic']]
mixed_high_qual_conserved = [d for d in merged if d['quality'] >= 5000 and d['phyloP'] > 1]
mixed_low_qual_not_conserved = [d for d in merged if d['quality'] < 5000 and d['phyloP'] <= 1]

print(f"\nQuadrant Analysis (threshold: quality=5000, phyloP=1):")
print(f"  Ancient Fossils (low Q, high C, syntenic):     {len(ancient_fossils):,} ({100*len(ancient_fossils)/len(merged):.1f}%)")
print(f"  Recent Insertions (high Q, low C, non-syn):    {len(recent_insertions):,} ({100*len(recent_insertions)/len(merged):.1f}%)")
print(f"  High quality + conserved:                      {len(mixed_high_qual_conserved):,} ({100*len(mixed_high_qual_conserved)/len(merged):.1f}%)")
print(f"  Low quality + not conserved:                   {len(mixed_low_qual_not_conserved):,} ({100*len(mixed_low_qual_not_conserved)/len(merged):.1f}%)")

# Statistics by group
print(f"\nMean statistics:")
print(f"  Ancient Fossils:    quality={np.mean([d['quality'] for d in ancient_fossils]):.0f}, phyloP={np.mean([d['phyloP'] for d in ancient_fossils]):.2f}")
if recent_insertions:
    print(f"  Recent Insertions:  quality={np.mean([d['quality'] for d in recent_insertions]):.0f}, phyloP={np.mean([d['phyloP'] for d in recent_insertions]):.2f}")

# Correlation within groups
all_qual = [d['quality'] for d in merged]
all_cons = [d['phyloP'] for d in merged]
corr_all = np.corrcoef(all_qual, all_cons)[0, 1]

syn_qual = [d['quality'] for d in syntenic]
syn_cons = [d['phyloP'] for d in syntenic]
corr_syn = np.corrcoef(syn_qual, syn_cons)[0, 1]

nonsyn_qual = [d['quality'] for d in non_syntenic]
nonsyn_cons = [d['phyloP'] for d in non_syntenic]
corr_nonsyn = np.corrcoef(nonsyn_qual, nonsyn_cons)[0, 1]

print(f"\nCorrelation (quality vs conservation):")
print(f"  All hits:      r = {corr_all:.3f}")
print(f"  Syntenic:      r = {corr_syn:.3f}")
print(f"  Non-syntenic:  r = {corr_nonsyn:.3f}")

# Create second figure: 2D density plot for clearer visualization
fig2, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: Syntenic hits density
ax = axes[0]
x = [d['quality'] for d in syntenic]
y = [d['phyloP'] for d in syntenic]
h = ax.hist2d(x, y, bins=50, cmap='Greens', cmin=1)
ax.axhline(1, color='white', linestyle='--', alpha=0.8, linewidth=2)
ax.axvline(5000, color='white', linestyle='--', alpha=0.8, linewidth=2)
ax.set_xlabel('BLAST Quality (Identity % × Length)', fontsize=11)
ax.set_ylabel('Conservation Score (phyloP)', fontsize=11)
ax.set_title(f'Syntenic Hits (n={len(syntenic):,})\n"Ancient Fossils"', fontsize=12, fontweight='bold')
plt.colorbar(h[3], ax=ax, label='Count')

# Right: Non-syntenic hits density
ax = axes[1]
x = [d['quality'] for d in non_syntenic]
y = [d['phyloP'] for d in non_syntenic]
h = ax.hist2d(x, y, bins=50, cmap='Reds', cmin=1)
ax.axhline(1, color='white', linestyle='--', alpha=0.8, linewidth=2)
ax.axvline(5000, color='white', linestyle='--', alpha=0.8, linewidth=2)
ax.set_xlabel('BLAST Quality (Identity % × Length)', fontsize=11)
ax.set_ylabel('Conservation Score (phyloP)', fontsize=11)
ax.set_title(f'Non-Syntenic Hits (n={len(non_syntenic):,})\n"Recent Insertions"', fontsize=12, fontweight='bold')
plt.colorbar(h[3], ax=ax, label='Count')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/23_quality_conservation_density.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/23_quality_conservation_density.png")

plt.close('all')
