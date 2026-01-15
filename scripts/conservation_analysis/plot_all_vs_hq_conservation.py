#!/usr/bin/env python3
"""Compare conservation between ALL hits and high-confidence hits."""

import matplotlib.pyplot as plt
import numpy as np

def parse_conservation_file(filepath, filter_func=None):
    """Parse bigWigAverageOverBed output with optional filtering."""
    scores = []
    metadata = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            mean_score = float(parts[5])

            # Parse name: prefix_num|transcript|te|pident|length
            name_parts = name.split('|')
            category = name_parts[0].split('_')[0]
            pident = float(name_parts[3])
            length = int(name_parts[4])

            if filter_func and not filter_func(pident, length):
                continue

            scores.append(mean_score)
            metadata.append({'pident': pident, 'length': length, 'category': category})

    return scores, metadata

# Load data with different filters
print("Loading data...")
all_scores, all_meta = parse_conservation_file('results/repeatmasker_analysis/te_hits_all_conservation.tab')
hq_scores, hq_meta = parse_conservation_file('results/repeatmasker_analysis/te_hits_all_conservation.tab',
                                              filter_func=lambda p, l: p >= 80 and l >= 50)
short_scores, short_meta = parse_conservation_file('results/repeatmasker_analysis/te_hits_all_conservation.tab',
                                                    filter_func=lambda p, l: l < 30)
long_scores, long_meta = parse_conservation_file('results/repeatmasker_analysis/te_hits_all_conservation.tab',
                                                  filter_func=lambda p, l: l >= 100)

print(f"All hits: {len(all_scores):,}")
print(f"HQ (≥80%, ≥50bp): {len(hq_scores):,}")
print(f"Short (<30bp): {len(short_scores):,}")
print(f"Long (≥100bp): {len(long_scores):,}")

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1. Distribution comparison: All vs HQ
ax1 = axes[0, 0]
bins = np.linspace(-2, 5, 50)
ax1.hist(all_scores, bins=bins, alpha=0.6, label=f'All hits (n={len(all_scores):,})', density=True, color='blue')
ax1.hist(hq_scores, bins=bins, alpha=0.6, label=f'HQ ≥80%/≥50bp (n={len(hq_scores):,})', density=True, color='red')
ax1.axvline(0, color='black', linestyle='--', linewidth=1, label='Neutral')
ax1.axvline(np.median(all_scores), color='blue', linestyle=':', linewidth=2, label=f'All median ({np.median(all_scores):.2f})')
ax1.axvline(np.median(hq_scores), color='red', linestyle=':', linewidth=2, label=f'HQ median ({np.median(hq_scores):.2f})')
ax1.set_xlabel('phyloP Conservation Score')
ax1.set_ylabel('Density')
ax1.set_title('Conservation: All Hits vs High-Quality Hits')
ax1.legend(fontsize=8)
ax1.set_xlim(-2, 5)

# 2. Conservation by length threshold
ax2 = axes[0, 1]
length_thresholds = [20, 30, 40, 50, 75, 100, 150, 200]
means = []
medians = []
pct_conserved = []
n_hits = []

for thresh in length_thresholds:
    subset = [s for s, m in zip(all_scores, all_meta) if m['length'] >= thresh]
    if len(subset) >= 100:
        means.append(np.mean(subset))
        medians.append(np.median(subset))
        pct_conserved.append(100 * sum(1 for s in subset if s >= 1) / len(subset))
        n_hits.append(len(subset))
    else:
        means.append(np.nan)
        medians.append(np.nan)
        pct_conserved.append(np.nan)
        n_hits.append(0)

ax2.plot(length_thresholds, means, 'o-', label='Mean', markersize=8, linewidth=2)
ax2.plot(length_thresholds, medians, 's-', label='Median', markersize=8, linewidth=2)
ax2.axhline(1, color='green', linestyle='--', alpha=0.7, label='Conserved threshold')
ax2.set_xlabel('Minimum Alignment Length (bp)')
ax2.set_ylabel('Conservation Score')
ax2.set_title('Conservation Decreases with Hit Length')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Add N hits as text
for i, (t, n) in enumerate(zip(length_thresholds, n_hits)):
    if n > 0:
        ax2.annotate(f'n={n:,}', (t, means[i]), textcoords="offset points",
                    xytext=(0, 10), ha='center', fontsize=7)

# 3. Conservation by identity threshold
ax3 = axes[1, 0]
id_thresholds = [60, 65, 70, 75, 80, 85, 90, 95]
id_means = []
id_medians = []
id_n_hits = []

for thresh in id_thresholds:
    subset = [s for s, m in zip(all_scores, all_meta) if m['pident'] >= thresh]
    if len(subset) >= 100:
        id_means.append(np.mean(subset))
        id_medians.append(np.median(subset))
        id_n_hits.append(len(subset))
    else:
        id_means.append(np.nan)
        id_medians.append(np.nan)
        id_n_hits.append(0)

ax3.plot(id_thresholds, id_means, 'o-', label='Mean', markersize=8, linewidth=2, color='purple')
ax3.plot(id_thresholds, id_medians, 's-', label='Median', markersize=8, linewidth=2, color='orange')
ax3.axhline(1, color='green', linestyle='--', alpha=0.7, label='Conserved threshold')
ax3.set_xlabel('Minimum Percent Identity')
ax3.set_ylabel('Conservation Score')
ax3.set_title('Conservation Decreases with Higher Identity')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Short vs Long hits distribution
ax4 = axes[1, 1]
ax4.hist(short_scores, bins=bins, alpha=0.6, label=f'Short <30bp (n={len(short_scores):,})', density=True, color='green')
ax4.hist(long_scores, bins=bins, alpha=0.6, label=f'Long ≥100bp (n={len(long_scores):,})', density=True, color='purple')
ax4.axvline(0, color='black', linestyle='--', linewidth=1)
ax4.axvline(np.median(short_scores), color='green', linestyle=':', linewidth=2, label=f'Short median ({np.median(short_scores):.2f})')
ax4.axvline(np.median(long_scores), color='purple', linestyle=':', linewidth=2, label=f'Long median ({np.median(long_scores):.2f})')
ax4.set_xlabel('phyloP Conservation Score')
ax4.set_ylabel('Density')
ax4.set_title('Short Hits Are MORE Conserved Than Long Hits')
ax4.legend(fontsize=8)
ax4.set_xlim(-2, 5)

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/15_all_vs_hq_conservation.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/15_all_vs_hq_conservation.png")

# Create summary statistics figure
fig2, ax = plt.subplots(figsize=(12, 8))

# Stacked bar showing conservation categories by quality tier
quality_tiers = [
    ('All hits', lambda p, l: True),
    ('≥70% id', lambda p, l: p >= 70),
    ('≥80% id', lambda p, l: p >= 80),
    ('≥90% id', lambda p, l: p >= 90),
    ('≥80%/≥30bp', lambda p, l: p >= 80 and l >= 30),
    ('≥80%/≥50bp', lambda p, l: p >= 80 and l >= 50),
    ('≥80%/≥100bp', lambda p, l: p >= 80 and l >= 100),
]

categories = ['Fast-evolving (<0)', 'Weakly conserved (0-1)', 'Conserved (1-2)', 'Highly conserved (≥2)']
colors = ['#e74c3c', '#f39c12', '#27ae60', '#2980b9']

bar_data = []
tier_labels = []
tier_counts = []

for label, filt in quality_tiers:
    subset = [s for s, m in zip(all_scores, all_meta) if filt(m['pident'], m['length'])]
    n = len(subset)
    tier_labels.append(label)
    tier_counts.append(n)

    fast = sum(1 for s in subset if s < 0) / n * 100
    weak = sum(1 for s in subset if 0 <= s < 1) / n * 100
    cons = sum(1 for s in subset if 1 <= s < 2) / n * 100
    high = sum(1 for s in subset if s >= 2) / n * 100
    bar_data.append([fast, weak, cons, high])

bar_data = np.array(bar_data)
x = np.arange(len(tier_labels))
width = 0.7

bottom = np.zeros(len(tier_labels))
for i, (cat, color) in enumerate(zip(categories, colors)):
    ax.bar(x, bar_data[:, i], width, bottom=bottom, label=cat, color=color)
    bottom += bar_data[:, i]

ax.set_xticks(x)
ax.set_xticklabels([f'{l}\n(n={n:,})' for l, n in zip(tier_labels, tier_counts)], fontsize=9)
ax.set_ylabel('Percentage')
ax.set_title('Conservation Categories by Quality Threshold\n(Stricter filters = LESS conserved hits)')
ax.legend(loc='upper right')
ax.set_ylim(0, 100)

# Add percentage labels for conserved (>1)
for i, (label, n) in enumerate(zip(tier_labels, tier_counts)):
    pct_cons = bar_data[i, 2] + bar_data[i, 3]
    ax.text(i, pct_cons + bottom[i] - bar_data[i, 3] - bar_data[i, 2] + 2,
            f'{pct_cons:.0f}%', ha='center', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/16_conservation_by_quality_tier.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/16_conservation_by_quality_tier.png")

plt.show()
