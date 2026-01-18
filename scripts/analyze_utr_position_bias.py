#!/usr/bin/env python3
"""
Analyze positional bias of TE hits on 3'UTRs.

Questions:
- Do TE hits preferentially occur at start/middle/end of 3'UTRs?
- Is there a difference between real and shuffled?
- Does this vary by UTR length?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

print("=" * 70)
print("3'UTR POSITIONAL BIAS ANALYSIS")
print("=" * 70)

# Load real hits
print("\nLoading real hits...")
real_positions = []  # (relative_start, relative_end, qlen)
with open('results/genome_wide_all_3utrs.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        qstart = int(parts[6])
        qend = int(parts[7])
        qlen = int(parts[12])

        # Relative positions (0-1 scale)
        rel_start = min(qstart, qend) / qlen
        rel_end = max(qstart, qend) / qlen
        rel_mid = (rel_start + rel_end) / 2

        real_positions.append({
            'rel_start': rel_start,
            'rel_end': rel_end,
            'rel_mid': rel_mid,
            'qlen': qlen,
            'hit_len': abs(qend - qstart) + 1
        })

print(f"  Loaded {len(real_positions):,} real hits")

# Load shuffled hits (all 10 reps)
print("\nLoading shuffled hits...")
shuf_positions = []
for rep in range(1, 11):
    count = 0
    with open(f'results/shuffled_full/replicate_{rep:02d}_blast.tsv') as f:
        for line in f:
            parts = line.strip().split('\t')
            qstart = int(parts[6])
            qend = int(parts[7])
            qlen = int(parts[12])

            rel_start = min(qstart, qend) / qlen
            rel_end = max(qstart, qend) / qlen
            rel_mid = (rel_start + rel_end) / 2

            shuf_positions.append({
                'rel_start': rel_start,
                'rel_end': rel_end,
                'rel_mid': rel_mid,
                'qlen': qlen,
                'hit_len': abs(qend - qstart) + 1
            })
            count += 1
    print(f"  Rep {rep}: {count:,} hits")

print(f"  Total shuffled: {len(shuf_positions):,}")

# Convert to arrays
real_mids = np.array([p['rel_mid'] for p in real_positions])
real_starts = np.array([p['rel_start'] for p in real_positions])
real_ends = np.array([p['rel_end'] for p in real_positions])
real_qlens = np.array([p['qlen'] for p in real_positions])

shuf_mids = np.array([p['rel_mid'] for p in shuf_positions])
shuf_starts = np.array([p['rel_start'] for p in shuf_positions])
shuf_ends = np.array([p['rel_end'] for p in shuf_positions])
shuf_qlens = np.array([p['qlen'] for p in shuf_positions])

# Summary statistics
print("\n" + "=" * 70)
print("POSITIONAL SUMMARY (relative position 0=5' end, 1=3' end of UTR)")
print("=" * 70)

print(f"\n{'Metric':<30} {'Real':>15} {'Shuffled':>15}")
print("-" * 60)
print(f"{'Mean midpoint':<30} {np.mean(real_mids):>15.3f} {np.mean(shuf_mids):>15.3f}")
print(f"{'Median midpoint':<30} {np.median(real_mids):>15.3f} {np.median(shuf_mids):>15.3f}")
print(f"{'Std midpoint':<30} {np.std(real_mids):>15.3f} {np.std(shuf_mids):>15.3f}")

# Quartile analysis
print("\n" + "=" * 70)
print("QUARTILE ANALYSIS (by hit midpoint)")
print("=" * 70)

quartile_edges = [0, 0.25, 0.5, 0.75, 1.0]
quartile_names = ['Q1 (0-25%)', 'Q2 (25-50%)', 'Q3 (50-75%)', 'Q4 (75-100%)']

print(f"\n{'Quartile':<20} {'Real %':>12} {'Shuffled %':>12} {'Ratio':>10}")
print("-" * 55)

for i, name in enumerate(quartile_names):
    lo, hi = quartile_edges[i], quartile_edges[i+1]
    real_pct = 100 * np.mean((real_mids >= lo) & (real_mids < hi))
    shuf_pct = 100 * np.mean((shuf_mids >= lo) & (shuf_mids < hi))
    ratio = real_pct / shuf_pct if shuf_pct > 0 else float('inf')
    print(f"{name:<20} {real_pct:>11.2f}% {shuf_pct:>11.2f}% {ratio:>10.2f}x")

# Decile analysis
print("\n" + "=" * 70)
print("DECILE ANALYSIS (by hit midpoint)")
print("=" * 70)

decile_edges = np.linspace(0, 1, 11)
print(f"\n{'Decile':<15} {'Real %':>10} {'Shuf %':>10} {'Ratio':>8} {'Bar'}")
print("-" * 60)

real_decile_pcts = []
shuf_decile_pcts = []

for i in range(10):
    lo, hi = decile_edges[i], decile_edges[i+1]
    real_pct = 100 * np.mean((real_mids >= lo) & (real_mids < hi))
    shuf_pct = 100 * np.mean((shuf_mids >= lo) & (shuf_mids < hi))
    ratio = real_pct / shuf_pct if shuf_pct > 0 else float('inf')

    real_decile_pcts.append(real_pct)
    shuf_decile_pcts.append(shuf_pct)

    # Visual bar
    bar_len = int(real_pct * 3)
    bar = '█' * bar_len

    print(f"D{i+1} ({lo:.1f}-{hi:.1f})" + f"{real_pct:>10.2f}% {shuf_pct:>10.2f}% {ratio:>7.2f}x {bar}")

# Analysis by UTR length
print("\n" + "=" * 70)
print("POSITIONAL BIAS BY UTR LENGTH")
print("=" * 70)

len_bins = [(0, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 10000)]
len_names = ['<200bp', '200-500bp', '500-1000bp', '1000-2000bp', '>2000bp']

print(f"\n{'UTR Length':<15} {'Real Mean':>12} {'Shuf Mean':>12} {'Real n':>12} {'Shuf n':>12}")
print("-" * 65)

for (lo, hi), name in zip(len_bins, len_names):
    real_mask = (real_qlens >= lo) & (real_qlens < hi)
    shuf_mask = (shuf_qlens >= lo) & (shuf_qlens < hi)

    real_mean = np.mean(real_mids[real_mask]) if np.sum(real_mask) > 0 else 0
    shuf_mean = np.mean(shuf_mids[shuf_mask]) if np.sum(shuf_mask) > 0 else 0

    print(f"{name:<15} {real_mean:>12.3f} {shuf_mean:>12.3f} {np.sum(real_mask):>12,} {np.sum(shuf_mask):>12,}")

# Create figures
print("\nGenerating figures...")
fig_dir = Path('figures/shuffle_convergence')
fig_dir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel A: Midpoint distribution
ax = axes[0, 0]
bins = np.linspace(0, 1, 51)
ax.hist(real_mids, bins=bins, alpha=0.7, density=True, label='Real', color='steelblue', edgecolor='white')
ax.hist(shuf_mids, bins=bins, alpha=0.5, density=True, label='Shuffled', color='coral', edgecolor='white')
ax.axvline(0.5, color='gray', linestyle='--', linewidth=1)
ax.set_xlabel('Relative Position on UTR (0=5\' end, 1=3\' end)')
ax.set_ylabel('Density')
ax.set_title('A. Hit Midpoint Distribution')
ax.legend()

# Panel B: Decile comparison
ax = axes[0, 1]
x = np.arange(10)
width = 0.35
ax.bar(x - width/2, real_decile_pcts, width, label='Real', color='steelblue', edgecolor='black')
ax.bar(x + width/2, shuf_decile_pcts, width, label='Shuffled', color='coral', edgecolor='black')
ax.axhline(10, color='gray', linestyle='--', linewidth=1, label='Expected (uniform)')
ax.set_xlabel('Decile (1=5\' end, 10=3\' end)')
ax.set_ylabel('% of Hits')
ax.set_title('B. Decile Distribution')
ax.set_xticks(x)
ax.set_xticklabels([f'D{i+1}' for i in range(10)], fontsize=8)
ax.legend(fontsize=9)

# Panel C: Real/Shuffled ratio by decile
ax = axes[0, 2]
ratios = [r/s if s > 0 else 1 for r, s in zip(real_decile_pcts, shuf_decile_pcts)]
colors = ['steelblue' if r > 1 else 'coral' for r in ratios]
ax.bar(x, ratios, color=colors, edgecolor='black')
ax.axhline(1, color='gray', linestyle='--', linewidth=2)
ax.set_xlabel('Decile (1=5\' end, 10=3\' end)')
ax.set_ylabel('Real/Shuffled Ratio')
ax.set_title('C. Enrichment by Position')
ax.set_xticks(x)
ax.set_xticklabels([f'D{i+1}' for i in range(10)], fontsize=8)

# Panel D: By UTR length - heatmap style
ax = axes[1, 0]
# Calculate mean position for each length bin and decile
heatmap_data = np.zeros((len(len_bins), 10))
for i, (lo, hi) in enumerate(len_bins):
    real_mask = (real_qlens >= lo) & (real_qlens < hi)
    mids_subset = real_mids[real_mask]
    for j in range(10):
        d_lo, d_hi = decile_edges[j], decile_edges[j+1]
        heatmap_data[i, j] = 100 * np.mean((mids_subset >= d_lo) & (mids_subset < d_hi))

im = ax.imshow(heatmap_data, cmap='YlOrRd', aspect='auto', vmin=5, vmax=15)
ax.set_xticks(range(10))
ax.set_xticklabels([f'D{i+1}' for i in range(10)], fontsize=8)
ax.set_yticks(range(len(len_names)))
ax.set_yticklabels(len_names)
ax.set_xlabel('Position Decile')
ax.set_ylabel('UTR Length')
ax.set_title('D. Real Hits: Position by UTR Length')
plt.colorbar(im, ax=ax, label='% of hits')

# Panel E: Start vs End positions
ax = axes[1, 1]
# Sample for visualization
n_sample = min(50000, len(real_positions))
idx = np.random.choice(len(real_starts), n_sample, replace=False)
ax.scatter(real_starts[idx], real_ends[idx], alpha=0.1, s=1, c='steelblue', label='Real')
idx = np.random.choice(len(shuf_starts), n_sample, replace=False)
ax.scatter(shuf_starts[idx], shuf_ends[idx], alpha=0.1, s=1, c='coral', label='Shuffled')
ax.plot([0, 1], [0, 1], 'k--', linewidth=1)
ax.set_xlabel('Hit Start (relative)')
ax.set_ylabel('Hit End (relative)')
ax.set_title('E. Hit Start vs End Position')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Panel F: Summary text
ax = axes[1, 2]
ax.axis('off')

# Calculate key stats
real_5prime = 100 * np.mean(real_mids < 0.25)
shuf_5prime = 100 * np.mean(shuf_mids < 0.25)
real_3prime = 100 * np.mean(real_mids >= 0.75)
shuf_3prime = 100 * np.mean(shuf_mids >= 0.75)

# Pre-compute conditional text (can't use backslash in f-string expressions)
real_mean = np.mean(real_mids)
shuf_mean = np.mean(shuf_mids)
if real_mean > shuf_mean + 0.01:
    bias_text = "-> Real biased toward 3' end"
elif real_mean < shuf_mean - 0.01:
    bias_text = "-> Real biased toward 5' end"
else:
    bias_text = "-> Similar distribution"

if real_5prime/shuf_5prime > 1.1:
    interp_text = "Real hits show 5' bias (enriched near stop codon)"
elif real_3prime/shuf_3prime > 1.1:
    interp_text = "Real hits show 3' bias (enriched near poly-A)"
else:
    interp_text = "No strong positional bias detected"

summary = f"""
POSITIONAL BIAS SUMMARY

Position scale: 0 = 5' end (near stop codon)
                1 = 3' end (near poly-A)

MEAN MIDPOINT:
  Real:     {real_mean:.3f}
  Shuffled: {shuf_mean:.3f}
  {bias_text}

5' QUARTER (0-25%):
  Real:     {real_5prime:.1f}%
  Shuffled: {shuf_5prime:.1f}%
  Ratio:    {real_5prime/shuf_5prime:.2f}x

3' QUARTER (75-100%):
  Real:     {real_3prime:.1f}%
  Shuffled: {shuf_3prime:.1f}%
  Ratio:    {real_3prime/shuf_3prime:.2f}x

INTERPRETATION:
{interp_text}
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_dir / 'utr_position_bias.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"Saved: {fig_dir}/utr_position_bias.png")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
