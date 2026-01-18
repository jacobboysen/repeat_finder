#!/usr/bin/env python3
"""
Analyze and visualize alignment pattern similarity between real and shuffled TE hits.

This script examines whether overlapping real/shuffled hits represent the same
alignment pattern or distinct patterns hitting the same TE region.
"""

from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Colorblind-friendly palette
COLORS = {
    'primary': '#0077BB',
    'secondary': '#EE7733',
    'tertiary': '#009988',
    'accent': '#CC3311',
    'neutral': '#BBBBBB',
}

# Directories
shuffled_dir = Path('results/shuffled_controls/deduplicated')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("ALIGNMENT PATTERN COMPARISON: Real vs Shuffled")
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
                        'qseqid': parts[0],
                        'sseqid': parts[1],
                        'pident': float(parts[2]),
                        'length': int(parts[3]),
                        'qstart': int(parts[6]),
                        'qend': int(parts[7]),
                        'sstart': int(parts[8]),
                        'send': int(parts[9]),
                    })
                except (ValueError, IndexError):
                    continue
    return hits


def extract_transcript_id(qseqid):
    """Extract FBtr ID from query sequence ID."""
    if qseqid.startswith('FBtr'):
        return qseqid.split('_')[0]
    return qseqid


def ranges_overlap(s1, e1, s2, e2):
    """Check if two ranges overlap."""
    return max(s1, s2) <= min(e1, e2)


# Load real and shuffled data
print("\nLoading data...")
real_file = shuffled_dir / 'real_deduplicated_hits.tsv'
real_hits = load_blast_hits(real_file)
print(f"  Real hits: {len(real_hits):,}")

# Build lookup for real hits by transcript and TE
real_by_transcript_te = defaultdict(list)
for hit in real_hits:
    transcript = extract_transcript_id(hit['qseqid'])
    te = hit['sseqid']
    real_by_transcript_te[(transcript, te)].append(hit)

# Load shuffled data and find overlapping pairs
print("\nFinding overlapping hit pairs...")
overlapping_pairs = []

for i in range(1, 11):
    shuf_file = shuffled_dir / f'shuffled_rep{i}_deduplicated_hits.tsv'
    if not shuf_file.exists():
        continue

    shuf_hits = load_blast_hits(shuf_file)
    print(f"  Processing shuffled rep {i}: {len(shuf_hits):,} hits")

    for shuf_hit in shuf_hits:
        transcript = extract_transcript_id(shuf_hit['qseqid'])
        te = shuf_hit['sseqid']

        # Find matching real hits for this transcript/TE pair
        real_matches = real_by_transcript_te.get((transcript, te), [])

        for real_hit in real_matches:
            # Check for TE coordinate overlap
            s_start = min(shuf_hit['sstart'], shuf_hit['send'])
            s_end = max(shuf_hit['sstart'], shuf_hit['send'])
            r_start = min(real_hit['sstart'], real_hit['send'])
            r_end = max(real_hit['sstart'], real_hit['send'])

            if ranges_overlap(s_start, s_end, r_start, r_end):
                # Check query overlap
                q_overlap = ranges_overlap(
                    real_hit['qstart'], real_hit['qend'],
                    shuf_hit['qstart'], shuf_hit['qend']
                )

                overlapping_pairs.append({
                    'real_pident': real_hit['pident'],
                    'shuf_pident': shuf_hit['pident'],
                    'real_length': real_hit['length'],
                    'shuf_length': shuf_hit['length'],
                    'query_overlap': q_overlap,
                    'real_qstart': real_hit['qstart'],
                    'real_qend': real_hit['qend'],
                    'shuf_qstart': shuf_hit['qstart'],
                    'shuf_qend': shuf_hit['qend'],
                    'exact_pos': (r_start == s_start and r_end == s_end),
                })

print(f"\nTotal overlapping pairs: {len(overlapping_pairs):,}")

if len(overlapping_pairs) == 0:
    print("No overlapping pairs found!")
    exit(1)

# Calculate statistics
real_pidents = np.array([p['real_pident'] for p in overlapping_pairs])
shuf_pidents = np.array([p['shuf_pident'] for p in overlapping_pairs])
identity_diffs = np.abs(real_pidents - shuf_pidents)
query_overlaps = np.array([p['query_overlap'] for p in overlapping_pairs])
exact_matches = sum(1 for p in overlapping_pairs if p['exact_pos'])

# Correlation
correlation = np.corrcoef(real_pidents, shuf_pidents)[0, 1]

print(f"\n{'='*70}")
print("STATISTICS")
print(f"{'='*70}")
print(f"Overlapping pairs: {len(overlapping_pairs):,}")
print(f"Exact position matches: {exact_matches} ({100*exact_matches/len(overlapping_pairs):.2f}%)")
print(f"Identity correlation (r): {correlation:.3f}")
print(f"Mean identity difference: {np.mean(identity_diffs):.2f}%")
print(f"Median identity difference: {np.median(identity_diffs):.2f}%")
print(f"Query overlap rate: {100*np.mean(query_overlaps):.1f}%")

# Create visualization
print("\nGenerating Figure 33: Alignment Pattern Comparison...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Identity scatter plot
ax = axes[0, 0]
# Sample for plotting if too many
if len(overlapping_pairs) > 5000:
    idx = np.random.choice(len(overlapping_pairs), 5000, replace=False)
    plot_real = real_pidents[idx]
    plot_shuf = shuf_pidents[idx]
else:
    plot_real = real_pidents
    plot_shuf = shuf_pidents

ax.scatter(plot_real, plot_shuf, alpha=0.2, s=10, c=COLORS['primary'])
ax.plot([60, 100], [60, 100], 'k--', linewidth=1.5, alpha=0.7, label='y=x')
ax.set_xlabel('Real Hit Identity (%)', fontsize=11)
ax.set_ylabel('Shuffled Hit Identity (%)', fontsize=11)
ax.set_title(f'A. Identity Correlation (r={correlation:.3f})\n(same transcript, overlapping TE region)',
             fontsize=11, fontweight='bold', loc='left')
ax.set_xlim(60, 100)
ax.set_ylim(60, 100)
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Identity difference histogram
ax = axes[0, 1]
bins = np.linspace(0, 30, 31)
ax.hist(identity_diffs, bins=bins, color=COLORS['tertiary'], edgecolor='white', alpha=0.8)
ax.axvline(np.mean(identity_diffs), color=COLORS['accent'], linestyle='--', linewidth=2,
           label=f'Mean: {np.mean(identity_diffs):.1f}%')
ax.axvline(np.median(identity_diffs), color=COLORS['primary'], linestyle=':', linewidth=2,
           label=f'Median: {np.median(identity_diffs):.1f}%')
ax.set_xlabel('|Real Identity - Shuffled Identity| (%)', fontsize=11)
ax.set_ylabel('Number of Pairs', fontsize=11)
ax.set_title('B. Identity Difference Distribution', fontsize=11, fontweight='bold', loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add text annotation
diff_lt5 = np.sum(identity_diffs < 5)
diff_5to10 = np.sum((identity_diffs >= 5) & (identity_diffs < 10))
diff_ge10 = np.sum(identity_diffs >= 10)
total = len(identity_diffs)
ax.text(0.97, 0.97, f'<5%: {100*diff_lt5/total:.1f}%\n5-10%: {100*diff_5to10/total:.1f}%\n≥10%: {100*diff_ge10/total:.1f}%',
        transform=ax.transAxes, va='top', ha='right', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Panel C: Length comparison (real vs shuffled)
ax = axes[1, 0]
real_lengths = np.array([p['real_length'] for p in overlapping_pairs])
shuf_lengths = np.array([p['shuf_length'] for p in overlapping_pairs])

# Sample for plotting
if len(overlapping_pairs) > 5000:
    plot_real_len = real_lengths[idx]
    plot_shuf_len = shuf_lengths[idx]
else:
    plot_real_len = real_lengths
    plot_shuf_len = shuf_lengths

ax.scatter(plot_real_len, plot_shuf_len, alpha=0.2, s=10, c=COLORS['tertiary'])
ax.plot([0, 300], [0, 300], 'k--', linewidth=1.5, alpha=0.7, label='y=x')
ax.set_xlabel('Real Hit Length (bp)', fontsize=11)
ax.set_ylabel('Shuffled Hit Length (bp)', fontsize=11)
len_corr = np.corrcoef(real_lengths, shuf_lengths)[0, 1]
ax.set_title(f'C. Length Comparison (r={len_corr:.3f})', fontsize=11, fontweight='bold', loc='left')
ax.set_xlim(0, 300)
ax.set_ylim(0, 300)
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Summary diagram
ax = axes[1, 1]
ax.axis('off')

summary_text = f"""
ALIGNMENT PATTERN ANALYSIS SUMMARY

Overlapping Pairs Analyzed: {len(overlapping_pairs):,}

Key Findings:

1. MODEST IDENTITY CORRELATION
   Correlation: r = {correlation:.3f} (r² = {correlation**2:.2f})
   → ~{100*correlation**2:.0f}% of variance explained
   → Certain TE regions produce similar scores from any AU-rich input

2. IDENTITY DIFFERENCES
   Mean difference: {np.mean(identity_diffs):.1f}%
   {100*diff_lt5/total:.0f}% have <5% difference
   {100*diff_ge10/total:.0f}% have ≥10% difference

3. LENGTH CORRELATION
   r = {len_corr:.3f}
   → Same TE regions produce similar-length alignments

4. ALMOST NO EXACT MATCHES
   Only {exact_matches} exact TE position matches ({100*exact_matches/total:.2f}%)

CRITICAL CAVEAT:
AU-rich 3'UTR sequences contain functional regulatory
elements (AREs, RBP sites). Shuffled controls preserve
AU-richness, so overlap may reflect:
  - Compositional artifacts (null)
  - TE-derived regulatory elements (exaptation)
  - Convergent evolution of AU-rich motifs
"""

ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig(figures_dir / '33_alignment_pattern_comparison.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/33_alignment_pattern_comparison.png")
plt.close()

print(f"\n{'='*70}")
print("ALIGNMENT PATTERN ANALYSIS COMPLETE")
print(f"{'='*70}")
