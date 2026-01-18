#!/usr/bin/env python3
"""
Plot the ~6.7% of real hits that are strong TE fossil candidates:
- Non-overlapping hits (no shuffled hit on same TE region)
- Hits where real contains shuffled (real alignment is larger)

Compare their quality (identity vs length) to the rest.
"""

from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

COLORS = {
    'fossil': '#CC3311',      # Red - fossil candidates
    'overlap': '#BBBBBB',     # Gray - overlapping hits
    'primary': '#0077BB',
}

shuffled_dir = Path('results/shuffled_controls/deduplicated')
figures_dir = Path('figures/repeatmasker_comparison')

print("=" * 70)
print("FOSSIL CANDIDATE QUALITY ANALYSIS")
print("=" * 70)


def load_blast_hits(filepath):
    """Load BLAST hits with key columns."""
    hits = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                try:
                    hits.append({
                        'sseqid': parts[1],
                        'pident': float(parts[2]),
                        'length': int(parts[3]),
                        'sstart': int(parts[8]),
                        'send': int(parts[9]),
                    })
                except (ValueError, IndexError):
                    continue
    return hits


def get_te_range(hit):
    """Get normalized TE range (start, end) where start < end."""
    return (min(hit['sstart'], hit['send']), max(hit['sstart'], hit['send']))


# Load real hits
print("\nLoading real hits...")
real_file = shuffled_dir / 'real_deduplicated_hits.tsv'
real_hits = load_blast_hits(real_file)
print(f"  Real hits: {len(real_hits):,}")

# Build index of shuffled hits by TE
print("\nLoading shuffled hits and building TE index...")
shuffled_by_te = defaultdict(list)

for i in range(1, 11):
    shuf_file = shuffled_dir / f'shuffled_rep{i}_deduplicated_hits.tsv'
    if not shuf_file.exists():
        continue

    shuf_hits = load_blast_hits(shuf_file)
    print(f"  Shuffled rep {i}: {len(shuf_hits):,} hits")

    for hit in shuf_hits:
        te = hit['sseqid']
        s_start, s_end = get_te_range(hit)
        shuffled_by_te[te].append((s_start, s_end, hit['length']))

print(f"  Unique TEs in shuffled: {len(shuffled_by_te):,}")

# Classify each real hit
print("\nClassifying real hits...")
non_overlapping = []
real_contains_shuffled = []
other_hits = []

for hit in real_hits:
    te = hit['sseqid']
    r_start, r_end = get_te_range(hit)
    r_len = hit['length']

    # Check against all shuffled hits on this TE
    shuf_hits_on_te = shuffled_by_te.get(te, [])

    if not shuf_hits_on_te:
        # No shuffled hits on this TE at all
        non_overlapping.append(hit)
        continue

    # Check for overlap with any shuffled hit
    has_overlap = False
    real_is_bigger = True  # Assume true until proven false

    for s_start, s_end, s_len in shuf_hits_on_te:
        # Check overlap
        if max(r_start, s_start) <= min(r_end, s_end):
            has_overlap = True
            # Check containment - does real fully contain shuffled?
            if not (r_start <= s_start and r_end >= s_end):
                real_is_bigger = False
            # If shuffled is bigger, real doesn't contain it
            if s_len > r_len:
                real_is_bigger = False

    if not has_overlap:
        non_overlapping.append(hit)
    elif real_is_bigger and has_overlap:
        real_contains_shuffled.append(hit)
    else:
        other_hits.append(hit)

# Combine fossil candidates
fossil_candidates = non_overlapping + real_contains_shuffled

print(f"\n{'='*70}")
print("CLASSIFICATION RESULTS")
print(f"{'='*70}")
print(f"Non-overlapping (no shuffled on TE region): {len(non_overlapping):,} ({100*len(non_overlapping)/len(real_hits):.1f}%)")
print(f"Real contains shuffled: {len(real_contains_shuffled):,} ({100*len(real_contains_shuffled)/len(real_hits):.1f}%)")
print(f"TOTAL FOSSIL CANDIDATES: {len(fossil_candidates):,} ({100*len(fossil_candidates)/len(real_hits):.1f}%)")
print(f"Other (shuffled overlap/contains): {len(other_hits):,} ({100*len(other_hits)/len(real_hits):.1f}%)")

# Extract metrics
fossil_pidents = np.array([h['pident'] for h in fossil_candidates])
fossil_lengths = np.array([h['length'] for h in fossil_candidates])
other_pidents = np.array([h['pident'] for h in other_hits])
other_lengths = np.array([h['length'] for h in other_hits])

# Quality statistics
print(f"\n{'='*70}")
print("QUALITY COMPARISON")
print(f"{'='*70}")

def quality_stats(pidents, lengths, label):
    hq = np.sum((pidents >= 80) & (lengths >= 50))
    vhq = np.sum((pidents >= 85) & (lengths >= 100))
    print(f"\n{label}:")
    print(f"  Median identity: {np.median(pidents):.1f}%")
    print(f"  Median length: {np.median(lengths):.0f} bp")
    print(f"  Mean identity: {np.mean(pidents):.1f}%")
    print(f"  Mean length: {np.mean(lengths):.0f} bp")
    print(f"  High-quality (≥80% id, ≥50bp): {hq:,} ({100*hq/len(pidents):.1f}%)")
    print(f"  Very high-quality (≥85% id, ≥100bp): {vhq:,} ({100*vhq/len(pidents):.1f}%)")
    return hq, vhq

fossil_hq, fossil_vhq = quality_stats(fossil_pidents, fossil_lengths, "Fossil Candidates")
other_hq, other_vhq = quality_stats(other_pidents, other_lengths, "Other Hits (shuffled overlap)")

# Create figure
print("\nGenerating figure...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: 2D histogram of fossil candidates
ax = axes[0, 0]
bins_pident = np.linspace(60, 100, 41)
bins_length = np.linspace(0, 300, 31)

h = ax.hist2d(fossil_pidents, np.minimum(fossil_lengths, 300),
              bins=[bins_pident, bins_length], cmap='Reds', vmin=1)
ax.axhline(50, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax.axvline(80, color='black', linestyle=':', linewidth=1.5, alpha=0.7)
ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Alignment Length (bp)', fontsize=11)
ax.set_title(f'A. Fossil Candidates (n={len(fossil_candidates):,}, {100*len(fossil_candidates)/len(real_hits):.1f}%)\n'
             f'Median: {np.median(fossil_pidents):.1f}% identity, {np.median(fossil_lengths):.0f}bp',
             fontsize=11, fontweight='bold', loc='left')
plt.colorbar(h[3], ax=ax, label='Count')

# Panel B: 2D histogram of other hits (for comparison)
ax = axes[0, 1]
# Sample if too large
if len(other_hits) > 100000:
    idx = np.random.choice(len(other_hits), 100000, replace=False)
    plot_other_pidents = other_pidents[idx]
    plot_other_lengths = other_lengths[idx]
else:
    plot_other_pidents = other_pidents
    plot_other_lengths = other_lengths

h = ax.hist2d(plot_other_pidents, np.minimum(plot_other_lengths, 300),
              bins=[bins_pident, bins_length], cmap='Greys', vmin=1)
ax.axhline(50, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax.axvline(80, color='black', linestyle=':', linewidth=1.5, alpha=0.7)
ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Alignment Length (bp)', fontsize=11)
ax.set_title(f'B. Other Hits - Shuffled Overlap (n={len(other_hits):,}, {100*len(other_hits)/len(real_hits):.1f}%)\n'
             f'Median: {np.median(other_pidents):.1f}% identity, {np.median(other_lengths):.0f}bp',
             fontsize=11, fontweight='bold', loc='left')
plt.colorbar(h[3], ax=ax, label='Count')

# Panel C: Identity distribution comparison
ax = axes[1, 0]
bins = np.linspace(60, 100, 41)
ax.hist(fossil_pidents, bins=bins, alpha=0.7, color=COLORS['fossil'],
        label=f'Fossil candidates (n={len(fossil_candidates):,})', density=True, edgecolor='white')
ax.hist(other_pidents, bins=bins, alpha=0.5, color=COLORS['overlap'],
        label=f'Other hits (n={len(other_hits):,})', density=True, edgecolor='white')
ax.axvline(80, color='black', linestyle=':', linewidth=1.5, alpha=0.7)
ax.set_xlabel('Percent Identity (%)', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title('C. Identity Distribution Comparison', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Quality summary
ax = axes[1, 1]
ax.axis('off')

summary_text = f"""
FOSSIL CANDIDATE QUALITY SUMMARY

Total real hits: {len(real_hits):,}

FOSSIL CANDIDATES: {len(fossil_candidates):,} ({100*len(fossil_candidates)/len(real_hits):.1f}%)
  - Non-overlapping: {len(non_overlapping):,}
  - Real contains shuffled: {len(real_contains_shuffled):,}

                        Fossil        Other
                        Candidates    (overlap)
Median Identity:        {np.median(fossil_pidents):>6.1f}%      {np.median(other_pidents):>6.1f}%
Median Length:          {np.median(fossil_lengths):>6.0f}bp     {np.median(other_lengths):>6.0f}bp
High-quality (≥80%,≥50bp): {100*fossil_hq/len(fossil_candidates):>5.1f}%      {100*other_hq/len(other_hits):>6.2f}%
Very HQ (≥85%,≥100bp):  {100*fossil_vhq/len(fossil_candidates):>6.2f}%      {100*other_vhq/len(other_hits):>6.2f}%

INTERPRETATION:
Fossil candidates show {'HIGHER' if np.median(fossil_pidents) > np.median(other_pidents) else 'SIMILAR'} identity
({np.median(fossil_pidents):.1f}% vs {np.median(other_pidents):.1f}%) and
{'SIMILAR' if abs(np.median(fossil_lengths) - np.median(other_lengths)) < 5 else 'DIFFERENT'} length distributions.

High-quality rate is {fossil_hq/len(fossil_candidates) / (other_hq/len(other_hits)):.1f}x higher
in fossil candidates.
"""

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig(figures_dir / '34_fossil_candidates_quality.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/34_fossil_candidates_quality.png")
plt.close()

print(f"\n{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
