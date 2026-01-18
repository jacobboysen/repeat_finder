#!/usr/bin/env python3
"""
Analyze overlap between real and shuffled hits on TE sequences.

Identifies which TE regions are hit by both real and shuffled sequences,
suggesting these regions may be due to sequence composition rather than
genuine TE fossils.
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

COLORS = {
    'real_only': '#0077BB',      # Blue
    'shuffled_only': '#EE7733',  # Orange
    'both': '#CC3311',           # Red
    'neutral': '#BBBBBB',        # Gray
}

shuffled_dir = Path('results/shuffled_controls/deduplicated')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("REAL vs SHUFFLED: TE Region Overlap Analysis")
print("=" * 70)

def load_te_regions(filepath):
    """Load TE hit regions from BLAST file.

    Returns dict: te_id -> list of (start, end, pident, length) tuples
    """
    te_regions = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                try:
                    te_id = parts[1]
                    pident = float(parts[2])
                    length = int(parts[3])
                    sstart = int(parts[8])
                    send = int(parts[9])
                    # Normalize to min/max
                    te_start = min(sstart, send)
                    te_end = max(sstart, send)
                    te_regions[te_id].append((te_start, te_end, pident, length))
                except (ValueError, IndexError):
                    continue
    return te_regions

def regions_overlap(r1, r2):
    """Check if two regions overlap."""
    return r1[0] <= r2[1] and r2[0] <= r1[1]

def merge_overlapping_regions(regions):
    """Merge overlapping regions into non-overlapping intervals."""
    if not regions:
        return []
    # Sort by start position
    sorted_regions = sorted(regions, key=lambda x: x[0])
    merged = [sorted_regions[0][:2]]  # Just start, end

    for start, end, _, _ in sorted_regions[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:  # Overlapping
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged

# Load real data
print("\nLoading real data...")
real_file = shuffled_dir / 'real_deduplicated_hits.tsv'
real_regions = load_te_regions(real_file)
real_hit_count = sum(len(v) for v in real_regions.values())
print(f"  Real: {real_hit_count:,} hits across {len(real_regions):,} TEs")

# Load and merge all shuffled data
print("\nLoading shuffled data (all 10 replicates)...")
shuffled_regions = defaultdict(list)
total_shuffled_hits = 0

for i in range(1, 11):
    shuf_file = shuffled_dir / f'shuffled_rep{i}_deduplicated_hits.tsv'
    if shuf_file.exists():
        regions = load_te_regions(shuf_file)
        for te_id, region_list in regions.items():
            shuffled_regions[te_id].extend(region_list)
        hit_count = sum(len(v) for v in regions.values())
        total_shuffled_hits += hit_count
        print(f"  Shuffled {i}: {hit_count:,} hits")

print(f"  Total shuffled: {total_shuffled_hits:,} hits across {len(shuffled_regions):,} TEs")

# Analyze overlap at TE level
print("\n" + "=" * 70)
print("TE-LEVEL OVERLAP ANALYSIS")
print("=" * 70)

real_only_tes = set(real_regions.keys()) - set(shuffled_regions.keys())
shuffled_only_tes = set(shuffled_regions.keys()) - set(real_regions.keys())
shared_tes = set(real_regions.keys()) & set(shuffled_regions.keys())

print(f"\nTEs with hits in REAL only: {len(real_only_tes):,}")
print(f"TEs with hits in SHUFFLED only: {len(shuffled_only_tes):,}")
print(f"TEs with hits in BOTH: {len(shared_tes):,}")

# For shared TEs, analyze region-level overlap
print("\n" + "=" * 70)
print("REGION-LEVEL OVERLAP ANALYSIS (on shared TEs)")
print("=" * 70)

overlap_stats = {
    'real_only_hits': 0,
    'shuffled_only_hits': 0,
    'overlapping_real_hits': 0,
    'overlapping_shuffled_hits': 0,
}

te_overlap_details = []

for te_id in shared_tes:
    real_hits = real_regions[te_id]
    shuf_hits = shuffled_regions[te_id]

    # Check each real hit for overlap with any shuffled hit
    real_overlapping = 0
    real_non_overlapping = 0

    for r_start, r_end, r_pident, r_length in real_hits:
        has_overlap = False
        for s_start, s_end, _, _ in shuf_hits:
            if regions_overlap((r_start, r_end), (s_start, s_end)):
                has_overlap = True
                break
        if has_overlap:
            real_overlapping += 1
        else:
            real_non_overlapping += 1

    # Check each shuffled hit for overlap with any real hit
    shuf_overlapping = 0
    shuf_non_overlapping = 0

    for s_start, s_end, s_pident, s_length in shuf_hits:
        has_overlap = False
        for r_start, r_end, _, _ in real_hits:
            if regions_overlap((r_start, r_end), (s_start, s_end)):
                has_overlap = True
                break
        if has_overlap:
            shuf_overlapping += 1
        else:
            shuf_non_overlapping += 1

    overlap_stats['real_only_hits'] += real_non_overlapping
    overlap_stats['overlapping_real_hits'] += real_overlapping
    overlap_stats['shuffled_only_hits'] += shuf_non_overlapping
    overlap_stats['overlapping_shuffled_hits'] += shuf_overlapping

    if real_overlapping > 0:
        te_overlap_details.append({
            'te_id': te_id,
            'real_total': len(real_hits),
            'real_overlapping': real_overlapping,
            'shuf_total': len(shuf_hits),
            'shuf_overlapping': shuf_overlapping,
        })

# Add hits from TEs only in real or only in shuffled
for te_id in real_only_tes:
    overlap_stats['real_only_hits'] += len(real_regions[te_id])

for te_id in shuffled_only_tes:
    overlap_stats['shuffled_only_hits'] += len(shuffled_regions[te_id])

print(f"\nReal hits with NO shuffled overlap: {overlap_stats['real_only_hits']:,}")
print(f"Real hits WITH shuffled overlap: {overlap_stats['overlapping_real_hits']:,}")
print(f"Shuffled hits with NO real overlap: {overlap_stats['shuffled_only_hits']:,}")
print(f"Shuffled hits WITH real overlap: {overlap_stats['overlapping_shuffled_hits']:,}")

pct_real_overlapping = 100 * overlap_stats['overlapping_real_hits'] / real_hit_count
print(f"\nPercent of REAL hits overlapping shuffled regions: {pct_real_overlapping:.2f}%")

# ============================================================
# FIGURE 1: Overview of TE and Hit Overlap
# ============================================================
print("\nGenerating figures...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Real vs Shuffled: TE Region Overlap Analysis', fontsize=14, fontweight='bold', y=0.98)

# Panel A: TE-level Venn-like bar chart
ax = axes[0, 0]
categories = ['Real Only', 'Both', 'Shuffled Only']
counts = [len(real_only_tes), len(shared_tes), len(shuffled_only_tes)]
colors = [COLORS['real_only'], COLORS['both'], COLORS['shuffled_only']]

bars = ax.bar(categories, counts, color=colors, edgecolor='black', linewidth=1)
ax.set_ylabel('Number of TEs', fontsize=11)
ax.set_title('A. TEs Hit by Real vs Shuffled', fontsize=12, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.02,
            f'{count:,}', ha='center', fontsize=10, fontweight='bold')

# Panel B: Hit-level overlap for shared TEs
ax = axes[0, 1]
categories = ['Real\n(no overlap)', 'Real\n(overlaps shuffled)', 'Shuffled\n(overlaps real)', 'Shuffled\n(no overlap)']
counts = [overlap_stats['real_only_hits'], overlap_stats['overlapping_real_hits'],
          overlap_stats['overlapping_shuffled_hits'], overlap_stats['shuffled_only_hits']]
colors = [COLORS['real_only'], COLORS['both'], COLORS['both'], COLORS['shuffled_only']]

bars = ax.bar(categories, counts, color=colors, edgecolor='black', linewidth=1)
ax.set_ylabel('Number of Hits', fontsize=11)
ax.set_title('B. Hit-Level Overlap', fontsize=12, fontweight='bold', loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.setp(ax.xaxis.get_majorticklabels(), fontsize=9)

for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.02,
            f'{count:,}', ha='center', fontsize=9, fontweight='bold')

# Panel C: Pie chart of real hit categories
ax = axes[1, 0]

real_categories = {
    'Real-only TEs': sum(len(real_regions[te]) for te in real_only_tes),
    'Shared TEs,\nno overlap': overlap_stats['real_only_hits'] - sum(len(real_regions[te]) for te in real_only_tes),
    'Overlaps shuffled': overlap_stats['overlapping_real_hits'],
}

# Recalculate to ensure correct
real_in_real_only_tes = sum(len(real_regions[te]) for te in real_only_tes)
real_in_shared_no_overlap = overlap_stats['real_only_hits'] - real_in_real_only_tes
real_overlapping = overlap_stats['overlapping_real_hits']

pie_data = [real_in_real_only_tes, real_in_shared_no_overlap, real_overlapping]
pie_labels = [f'Real-only TEs\n({real_in_real_only_tes:,})',
              f'Shared TEs, no overlap\n({real_in_shared_no_overlap:,})',
              f'Overlaps shuffled\n({real_overlapping:,})']
pie_colors = [COLORS['real_only'], COLORS['neutral'], COLORS['both']]

wedges, texts, autotexts = ax.pie(pie_data, labels=pie_labels, colors=pie_colors,
                                   autopct='%1.1f%%', startangle=90,
                                   wedgeprops=dict(edgecolor='black', linewidth=1),
                                   textprops={'fontsize': 9})

ax.set_title('C. Classification of Real Hits', fontsize=12, fontweight='bold', loc='left')

# Panel D: Top TEs with most overlap
ax = axes[1, 1]

# Sort by number of overlapping real hits
te_overlap_details.sort(key=lambda x: -x['real_overlapping'])
top_overlap = te_overlap_details[:20]

if top_overlap:
    te_names = [d['te_id'] for d in top_overlap]
    real_overlap_counts = [d['real_overlapping'] for d in top_overlap]
    real_total_counts = [d['real_total'] for d in top_overlap]

    y_pos = np.arange(len(te_names))

    # Background: total real hits
    ax.barh(y_pos, real_total_counts, color=COLORS['neutral'], edgecolor='black',
            linewidth=0.5, label='Real (total)', alpha=0.5)
    # Foreground: overlapping hits
    ax.barh(y_pos, real_overlap_counts, color=COLORS['both'], edgecolor='black',
            linewidth=0.5, label='Real (overlaps shuffled)')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(te_names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Number of Hits', fontsize=11)
    ax.set_title('D. Top 20 TEs with Most Real-Shuffled Overlap', fontsize=12, fontweight='bold', loc='left')
    ax.legend(fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(figures_dir / '31_real_vs_shuffled_te_overlap.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/31_real_vs_shuffled_te_overlap.png")
plt.close()

# ============================================================
# FIGURE 2: Quality comparison of overlapping vs non-overlapping hits
# ============================================================
print("\nAnalyzing quality of overlapping vs non-overlapping hits...")

# Collect quality metrics for different categories
real_overlapping_quality = []
real_non_overlapping_quality = []

for te_id in shared_tes:
    real_hits = real_regions[te_id]
    shuf_hits = shuffled_regions[te_id]

    for r_start, r_end, r_pident, r_length in real_hits:
        has_overlap = any(regions_overlap((r_start, r_end), (s_start, s_end))
                         for s_start, s_end, _, _ in shuf_hits)
        if has_overlap:
            real_overlapping_quality.append((r_pident, r_length))
        else:
            real_non_overlapping_quality.append((r_pident, r_length))

# Add hits from real-only TEs
for te_id in real_only_tes:
    for r_start, r_end, r_pident, r_length in real_regions[te_id]:
        real_non_overlapping_quality.append((r_pident, r_length))

print(f"  Overlapping hits: {len(real_overlapping_quality):,}")
print(f"  Non-overlapping hits: {len(real_non_overlapping_quality):,}")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Quality Comparison: Overlapping vs Non-Overlapping Real Hits', fontsize=13, fontweight='bold')

# Panel A: Identity distribution
ax = axes[0]
if real_overlapping_quality:
    ax.hist([q[0] for q in real_non_overlapping_quality], bins=40, alpha=0.7,
            color=COLORS['real_only'], edgecolor='white', label=f'Non-overlapping (n={len(real_non_overlapping_quality):,})')
    ax.hist([q[0] for q in real_overlapping_quality], bins=40, alpha=0.7,
            color=COLORS['both'], edgecolor='white', label=f'Overlapping (n={len(real_overlapping_quality):,})')
ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Number of Hits', fontsize=11)
ax.set_title('A. Identity Distribution', fontsize=12, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Length distribution
ax = axes[1]
if real_overlapping_quality:
    ax.hist([min(q[1], 300) for q in real_non_overlapping_quality], bins=40, alpha=0.7,
            color=COLORS['real_only'], edgecolor='white', label='Non-overlapping')
    ax.hist([min(q[1], 300) for q in real_overlapping_quality], bins=40, alpha=0.7,
            color=COLORS['both'], edgecolor='white', label='Overlapping')
ax.set_xlabel('Alignment Length (bp)', fontsize=11)
ax.set_ylabel('Number of Hits', fontsize=11)
ax.set_title('B. Length Distribution', fontsize=12, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: 2D comparison
ax = axes[2]
if real_overlapping_quality and real_non_overlapping_quality:
    # Sample for plotting
    max_pts = 5000

    non_overlap_sample = real_non_overlapping_quality
    if len(non_overlap_sample) > max_pts:
        idx = np.random.choice(len(non_overlap_sample), max_pts, replace=False)
        non_overlap_sample = [real_non_overlapping_quality[i] for i in idx]

    overlap_sample = real_overlapping_quality
    if len(overlap_sample) > max_pts:
        idx = np.random.choice(len(overlap_sample), max_pts, replace=False)
        overlap_sample = [real_overlapping_quality[i] for i in idx]

    ax.scatter([q[0] for q in non_overlap_sample], [min(q[1], 300) for q in non_overlap_sample],
               alpha=0.3, s=10, c=COLORS['real_only'], label='Non-overlapping')
    ax.scatter([q[0] for q in overlap_sample], [min(q[1], 300) for q in overlap_sample],
               alpha=0.3, s=10, c=COLORS['both'], label='Overlapping')

    ax.axhline(50, color='black', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline(80, color='black', linestyle=':', linewidth=1, alpha=0.7)

ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Alignment Length (bp)', fontsize=11)
ax.set_title('C. Identity vs Length', fontsize=12, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(figures_dir / '32_overlap_quality_comparison.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/32_overlap_quality_comparison.png")
plt.close()

# Summary statistics
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

if real_overlapping_quality and real_non_overlapping_quality:
    overlap_pidents = [q[0] for q in real_overlapping_quality]
    overlap_lengths = [q[1] for q in real_overlapping_quality]
    non_overlap_pidents = [q[0] for q in real_non_overlapping_quality]
    non_overlap_lengths = [q[1] for q in real_non_overlapping_quality]

    print(f"\n{'Metric':<25} {'Overlapping':>15} {'Non-overlapping':>15}")
    print("-" * 55)
    print(f"{'Count':<25} {len(real_overlapping_quality):>15,} {len(real_non_overlapping_quality):>15,}")
    print(f"{'Median Identity':<25} {np.median(overlap_pidents):>15.1f} {np.median(non_overlap_pidents):>15.1f}")
    print(f"{'Mean Identity':<25} {np.mean(overlap_pidents):>15.1f} {np.mean(non_overlap_pidents):>15.1f}")
    print(f"{'Median Length':<25} {np.median(overlap_lengths):>15.0f} {np.median(non_overlap_lengths):>15.0f}")
    print(f"{'Mean Length':<25} {np.mean(overlap_lengths):>15.0f} {np.mean(non_overlap_lengths):>15.0f}")

    # High-quality hits
    overlap_hq = sum(1 for p, l in real_overlapping_quality if p >= 80 and l >= 50)
    non_overlap_hq = sum(1 for p, l in real_non_overlapping_quality if p >= 80 and l >= 50)
    print(f"{'High-quality (≥80%, ≥50bp)':<25} {overlap_hq:>15,} {non_overlap_hq:>15,}")
    print(f"{'HQ percentage':<25} {100*overlap_hq/len(real_overlapping_quality):>14.2f}% {100*non_overlap_hq/len(real_non_overlapping_quality):>14.2f}%")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
