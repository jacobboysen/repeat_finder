#!/usr/bin/env python3
"""
Create quality vs conservation figure with better metrics.

Options:
1. Identical bases = (pident/100) × length  [interpretable: actual matching bp]
2. Length alone (simpler)
3. Identity alone (simpler)
"""

import matplotlib.pyplot as plt
import numpy as np

def load_data():
    """Load merged conservation and synteny data."""
    # Load conservation
    conservation = {}
    with open('results/repeatmasker_analysis/te_hits_all_conservation.tab') as f:
        for line in f:
            parts = line.strip().split('\t')
            conservation[parts[0]] = float(parts[5])

    # Load synteny with BLAST stats
    merged = []
    with open('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            name = parts[3]
            if name not in conservation:
                continue

            pident = float(parts[5])
            length = float(parts[6])
            sim_cov = float(parts[7])
            yak_cov = float(parts[8])
            ere_cov = float(parts[9])

            syntenic_count = sum(1 for c in [sim_cov, yak_cov, ere_cov] if c >= 0.5)

            merged.append({
                'pident': pident,
                'length': length,
                'identical_bp': (pident / 100) * length,  # Actual matching bases
                'phyloP': conservation[name],
                'syntenic_species': syntenic_count,
                'is_syntenic': syntenic_count >= 2
            })

    return merged

print("Loading data...")
data = load_data()
print(f"Loaded {len(data):,} hits")

syntenic = [d for d in data if d['is_syntenic']]
non_syntenic = [d for d in data if not d['is_syntenic']]

print(f"  Syntenic: {len(syntenic):,}")
print(f"  Non-syntenic: {len(non_syntenic):,}")

# Check value ranges
print(f"\nValue ranges:")
print(f"  pident: {min(d['pident'] for d in data):.1f} - {max(d['pident'] for d in data):.1f}")
print(f"  length: {min(d['length'] for d in data):.0f} - {max(d['length'] for d in data):.0f}")
print(f"  identical_bp: {min(d['identical_bp'] for d in data):.0f} - {max(d['identical_bp'] for d in data):.0f}")
print(f"  phyloP: {min(d['phyloP'] for d in data):.2f} - {max(d['phyloP'] for d in data):.2f}")

# Create 2x2 figure with different x-axis metrics
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

metrics = [
    ('identical_bp', 'Identical Bases (pident% × length / 100)'),
    ('length', 'Alignment Length (bp)'),
    ('pident', 'Percent Identity (%)'),
]

# Sample for plotting
max_pts = 8000
np.random.seed(42)

syn_idx = np.random.choice(len(syntenic), min(max_pts, len(syntenic)), replace=False)
syn_sample = [syntenic[i] for i in syn_idx]

nonsyn_idx = np.random.choice(len(non_syntenic), min(max_pts, len(non_syntenic)), replace=False)
nonsyn_sample = [non_syntenic[i] for i in nonsyn_idx]

for idx, (metric, label) in enumerate(metrics):
    ax = axes[idx // 2, idx % 2]

    x_syn = [d[metric] for d in syn_sample]
    y_syn = [d['phyloP'] for d in syn_sample]

    x_nonsyn = [d[metric] for d in nonsyn_sample]
    y_nonsyn = [d['phyloP'] for d in nonsyn_sample]

    ax.scatter(x_nonsyn, y_nonsyn, alpha=0.3, s=10, c='#e74c3c', label='Non-syntenic')
    ax.scatter(x_syn, y_syn, alpha=0.3, s=10, c='#27ae60', label='Syntenic')

    ax.axhline(1, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel(label)
    ax.set_ylabel('Conservation (phyloP)')
    ax.legend(loc='upper right')

    # Calculate correlation
    all_x = [d[metric] for d in data]
    all_y = [d['phyloP'] for d in data]
    corr = np.corrcoef(all_x, all_y)[0, 1]
    ax.set_title(f'{label}\nr = {corr:.3f}')

# 4th panel: 2D view of identical_bp
ax = axes[1, 1]

# Density plot with both populations
x_all = [d['identical_bp'] for d in data]
y_all = [d['phyloP'] for d in data]
synteny_status = [1 if d['is_syntenic'] else 0 for d in data]

# Create side-by-side comparison
h = ax.hist2d([d['identical_bp'] for d in syntenic],
              [d['phyloP'] for d in syntenic],
              bins=40, cmap='Greens', alpha=0.7, cmin=1)
ax.axhline(1, color='white', linestyle='--', linewidth=2)
ax.set_xlabel('Identical Bases')
ax.set_ylabel('Conservation (phyloP)')
ax.set_title('Syntenic Hits Density\n(Ancient Fossils)')
plt.colorbar(h[3], ax=ax, label='Count')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/24_quality_metrics_comparison.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/repeatmasker_comparison/24_quality_metrics_comparison.png")

# Now the clean figure for paper
fig2, ax = plt.subplots(figsize=(10, 8))

# Use identical_bp as x-axis
x_syn = [d['identical_bp'] for d in syn_sample]
y_syn = [d['phyloP'] for d in syn_sample]
x_nonsyn = [d['identical_bp'] for d in nonsyn_sample]
y_nonsyn = [d['phyloP'] for d in nonsyn_sample]

ax.scatter(x_nonsyn, y_nonsyn, alpha=0.4, s=15, c='#e74c3c',
           label=f'Non-syntenic (n={len(non_syntenic):,})', zorder=1)
ax.scatter(x_syn, y_syn, alpha=0.4, s=15, c='#27ae60',
           label=f'Syntenic ≥2 species (n={len(syntenic):,})', zorder=2)

ax.axhline(1, color='black', linestyle='--', alpha=0.6, linewidth=1.5, label='Conservation threshold')
ax.axvline(50, color='gray', linestyle=':', alpha=0.6, linewidth=1.5)

ax.set_xlabel('Identical Bases in Alignment', fontsize=12)
ax.set_ylabel('Conservation Score (phyloP)', fontsize=12)
ax.set_title('TE Hit Quality vs Conservation\nColored by Cross-Species Synteny', fontsize=14)
ax.legend(loc='upper right', fontsize=10)

# Add annotations
ax.text(0.02, 0.98, 'ANCIENT FOSSILS\nDiverged sequence\nHigh conservation\nSyntenic',
        transform=ax.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='#27ae60', alpha=0.2))
ax.text(0.98, 0.02, 'RECENT INSERTIONS\nHigh identity\nLow conservation\nMel-specific',
        transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='#e74c3c', alpha=0.2))

# Correlation annotation
corr = np.corrcoef([d['identical_bp'] for d in data], [d['phyloP'] for d in data])[0, 1]
ax.text(0.98, 0.98, f'r = {corr:.3f}', transform=ax.transAxes,
        fontsize=11, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/25_identical_bp_vs_conservation.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/25_identical_bp_vs_conservation.png")

# Statistics
print("\n" + "="*60)
print("STATISTICS BY METRIC")
print("="*60)

for metric, label in metrics:
    syn_vals = [d[metric] for d in syntenic]
    nonsyn_vals = [d[metric] for d in non_syntenic]
    print(f"\n{label}:")
    print(f"  Syntenic mean:     {np.mean(syn_vals):.1f}")
    print(f"  Non-syntenic mean: {np.mean(nonsyn_vals):.1f}")

print("\nConservation (phyloP):")
print(f"  Syntenic mean:     {np.mean([d['phyloP'] for d in syntenic]):.2f}")
print(f"  Non-syntenic mean: {np.mean([d['phyloP'] for d in non_syntenic]):.2f}")

plt.close('all')
