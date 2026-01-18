#!/usr/bin/env python3
"""
Plot real vs shuffled hit quality comparison.
Focus on identity and length distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

print("Loading data...")

# Load real hits
real_pidents = []
real_lengths = []
with open('results/genome_wide_all_3utrs.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        real_pidents.append(float(parts[2]))
        real_lengths.append(int(parts[3]))

real_pidents = np.array(real_pidents)
real_lengths = np.array(real_lengths)
print(f"  Real: {len(real_pidents):,} hits")

# Load shuffled (sample for plotting efficiency)
shuf_pidents = []
shuf_lengths = []
for rep in range(1, 11):
    with open(f'results/shuffled_full/replicate_{rep:02d}_blast.tsv') as f:
        for line in f:
            parts = line.strip().split('\t')
            shuf_pidents.append(float(parts[2]))
            shuf_lengths.append(int(parts[3]))

shuf_pidents = np.array(shuf_pidents)
shuf_lengths = np.array(shuf_lengths)
print(f"  Shuffled: {len(shuf_pidents):,} hits")

# Create figure
fig_dir = Path('figures/shuffle_convergence')
fig_dir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel A: Identity distribution
ax = axes[0, 0]
bins = np.linspace(60, 100, 41)
ax.hist(real_pidents, bins=bins, alpha=0.7, density=True, label=f'Real (n={len(real_pidents):,})',
        color='steelblue', edgecolor='white')
ax.hist(shuf_pidents, bins=bins, alpha=0.5, density=True, label=f'Shuffled (n={len(shuf_pidents):,})',
        color='coral', edgecolor='white')
ax.axvline(np.median(real_pidents), color='steelblue', linestyle='--', linewidth=2)
ax.axvline(np.median(shuf_pidents), color='coral', linestyle='--', linewidth=2)
ax.set_xlabel('Percent Identity (%)')
ax.set_ylabel('Density')
ax.set_title('A. Identity Distribution')
ax.legend(fontsize=9)

# Panel B: Length distribution
ax = axes[0, 1]
bins = np.linspace(0, 200, 41)
ax.hist(real_lengths, bins=bins, alpha=0.7, density=True, label='Real',
        color='steelblue', edgecolor='white')
ax.hist(shuf_lengths, bins=bins, alpha=0.5, density=True, label='Shuffled',
        color='coral', edgecolor='white')
ax.axvline(np.median(real_lengths), color='steelblue', linestyle='--', linewidth=2,
           label=f'Real median: {np.median(real_lengths):.0f}bp')
ax.axvline(np.median(shuf_lengths), color='coral', linestyle='--', linewidth=2,
           label=f'Shuf median: {np.median(shuf_lengths):.0f}bp')
ax.set_xlabel('Alignment Length (bp)')
ax.set_ylabel('Density')
ax.set_title('B. Length Distribution')
ax.legend(fontsize=9)

# Panel C: 2D histogram - Real
ax = axes[0, 2]
bins_pid = np.linspace(60, 100, 41)
bins_len = np.linspace(0, 200, 41)
h = ax.hist2d(real_pidents, np.minimum(real_lengths, 200), bins=[bins_pid, bins_len],
              cmap='Blues', density=True)
ax.axhline(50, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
ax.axvline(80, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
ax.set_xlabel('Percent Identity (%)')
ax.set_ylabel('Alignment Length (bp)')
ax.set_title('C. Real Hits: Identity vs Length')
plt.colorbar(h[3], ax=ax, label='Density')

# Panel D: 2D histogram - Shuffled
ax = axes[1, 0]
h = ax.hist2d(shuf_pidents, np.minimum(shuf_lengths, 200), bins=[bins_pid, bins_len],
              cmap='Oranges', density=True)
ax.axhline(50, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
ax.axvline(80, color='red', linestyle=':', linewidth=1.5, alpha=0.8)
ax.set_xlabel('Percent Identity (%)')
ax.set_ylabel('Alignment Length (bp)')
ax.set_title('D. Shuffled Hits: Identity vs Length')
plt.colorbar(h[3], ax=ax, label='Density')

# Panel E: High quality comparison
ax = axes[1, 1]

# Calculate HQ counts at different thresholds
thresholds = [(70, 30), (75, 40), (80, 50), (85, 75), (90, 100)]
real_counts = []
shuf_counts = []
labels = []

for min_id, min_len in thresholds:
    real_c = np.sum((real_pidents >= min_id) & (real_lengths >= min_len))
    shuf_c = np.sum((shuf_pidents >= min_id) & (shuf_lengths >= min_len)) / 10  # per replicate
    real_counts.append(real_c)
    shuf_counts.append(shuf_c)
    labels.append(f'≥{min_id}%\n≥{min_len}bp')

x = np.arange(len(labels))
width = 0.35
bars1 = ax.bar(x - width/2, real_counts, width, label='Real', color='steelblue', edgecolor='black')
bars2 = ax.bar(x + width/2, shuf_counts, width, label='Shuffled (per rep)', color='coral', edgecolor='black')

ax.set_ylabel('Number of Hits')
ax.set_title('E. High-Quality Hits by Threshold')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.legend()
ax.set_yscale('log')

# Add fold-change labels
for i, (r, s) in enumerate(zip(real_counts, shuf_counts)):
    if s > 0:
        fold = r / s
        ax.text(i, max(r, s) * 1.5, f'{fold:.0f}x', ha='center', fontsize=9, fontweight='bold')
    else:
        ax.text(i, r * 1.5, '∞', ha='center', fontsize=12, fontweight='bold')

# Panel F: Summary text
ax = axes[1, 2]
ax.axis('off')

# Calculate stats
real_hq = np.sum((real_pidents >= 80) & (real_lengths >= 50))
shuf_hq = np.sum((shuf_pidents >= 80) & (shuf_lengths >= 50))
real_vhq = np.sum((real_pidents >= 85) & (real_lengths >= 100))
shuf_vhq = np.sum((shuf_pidents >= 85) & (shuf_lengths >= 100))

summary = f"""
REAL vs SHUFFLED QUALITY COMPARISON

                          Real        Shuffled
Total hits:         {len(real_pidents):>10,}    {len(shuf_pidents):>10,}
Per replicate:      {len(real_pidents):>10,}    {len(shuf_pidents)//10:>10,}

IDENTITY:
  Mean:                 {np.mean(real_pidents):>6.1f}%        {np.mean(shuf_pidents):>6.1f}%
  Median:               {np.median(real_pidents):>6.1f}%        {np.median(shuf_pidents):>6.1f}%

LENGTH:
  Mean:                 {np.mean(real_lengths):>6.1f}bp       {np.mean(shuf_lengths):>6.1f}bp
  Median:               {np.median(real_lengths):>6.1f}bp       {np.median(shuf_lengths):>6.1f}bp

HIGH QUALITY:
  ≥80% & ≥50bp:      {real_hq:>10,}    {shuf_hq:>10,}
    (% of total)         {100*real_hq/len(real_pidents):>5.2f}%         {100*shuf_hq/len(shuf_pidents):>5.2f}%
    Fold enrichment:     {real_hq/(shuf_hq/10):.0f}x

  ≥85% & ≥100bp:     {real_vhq:>10,}    {shuf_vhq:>10,}
    Fold enrichment:     {"∞" if shuf_vhq == 0 else f"{real_vhq/(shuf_vhq/10):.0f}x"}

KEY FINDING:
Real hits are LONGER at similar identity.
The signal is in alignment LENGTH, not identity.
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_dir / 'real_vs_shuffled_quality.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"Saved: {fig_dir}/real_vs_shuffled_quality.png")

# Second figure: Length focus
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel 1: Length CDF
ax = axes[0]
real_sorted = np.sort(real_lengths)
shuf_sorted = np.sort(shuf_lengths)
ax.plot(real_sorted, np.linspace(0, 1, len(real_sorted)), label='Real', color='steelblue', linewidth=2)
ax.plot(shuf_sorted, np.linspace(0, 1, len(shuf_sorted)), label='Shuffled', color='coral', linewidth=2)
ax.axvline(50, color='gray', linestyle='--', alpha=0.7)
ax.axhline(0.5, color='gray', linestyle=':', alpha=0.7)
ax.set_xlabel('Alignment Length (bp)')
ax.set_ylabel('Cumulative Fraction')
ax.set_title('A. Length CDF')
ax.set_xlim(0, 300)
ax.legend()

# Calculate % above thresholds
for thresh in [50, 100, 150]:
    real_pct = 100 * np.mean(real_lengths >= thresh)
    shuf_pct = 100 * np.mean(shuf_lengths >= thresh)
    print(f"  ≥{thresh}bp: Real {real_pct:.1f}%, Shuffled {shuf_pct:.1f}%")

# Panel 2: Length difference at each identity level
ax = axes[1]
id_bins = np.arange(65, 100, 5)
real_len_by_id = []
shuf_len_by_id = []

for i, id_min in enumerate(id_bins[:-1]):
    id_max = id_bins[i+1]
    real_mask = (real_pidents >= id_min) & (real_pidents < id_max)
    shuf_mask = (shuf_pidents >= id_min) & (shuf_pidents < id_max)

    real_len_by_id.append(np.median(real_lengths[real_mask]) if np.sum(real_mask) > 0 else 0)
    shuf_len_by_id.append(np.median(shuf_lengths[shuf_mask]) if np.sum(shuf_mask) > 0 else 0)

x = (id_bins[:-1] + id_bins[1:]) / 2
ax.plot(x, real_len_by_id, 'o-', label='Real', color='steelblue', linewidth=2, markersize=8)
ax.plot(x, shuf_len_by_id, 's-', label='Shuffled', color='coral', linewidth=2, markersize=8)
ax.set_xlabel('Identity Bin (%)')
ax.set_ylabel('Median Length (bp)')
ax.set_title('B. Length by Identity Level')
ax.legend()

# Panel 3: Tail comparison (high quality)
ax = axes[2]
# Focus on lengths > 50bp
real_long = real_lengths[real_lengths >= 50]
shuf_long = shuf_lengths[shuf_lengths >= 50]

bins = np.linspace(50, 300, 26)
ax.hist(real_long, bins=bins, alpha=0.7, density=True, label=f'Real ≥50bp (n={len(real_long):,})',
        color='steelblue', edgecolor='white')
ax.hist(shuf_long, bins=bins, alpha=0.5, density=True, label=f'Shuffled ≥50bp (n={len(shuf_long):,})',
        color='coral', edgecolor='white')
ax.axvline(100, color='red', linestyle='--', linewidth=1.5)
ax.set_xlabel('Alignment Length (bp)')
ax.set_ylabel('Density')
ax.set_title('C. Long Hits Only (≥50bp)')
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(fig_dir / 'real_vs_shuffled_length_focus.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"Saved: {fig_dir}/real_vs_shuffled_length_focus.png")

print("\nDone!")
