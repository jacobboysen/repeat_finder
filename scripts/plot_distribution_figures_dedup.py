#!/usr/bin/env python3
"""
Generate distribution figures (01-14) using DEDUPLICATED 3'UTR data.

Creates publication-quality figures analyzing:
- Identity and length distributions
- TE family comparisons
- Conservation analyses
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Colorblind-friendly palette
COLORS = {
    'primary': '#0077BB',   # Blue
    'secondary': '#EE7733', # Orange
    'tertiary': '#009988',  # Teal
    'accent': '#CC3311',    # Red
    'neutral': '#BBBBBB',   # Gray
}

# Directories
blast_file = Path('results/3utr_deduplicated/3utr_deduplicated_hits.tsv')
figures_dir = Path('figures/repeatmasker_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("DISTRIBUTION FIGURES - DEDUPLICATED 3'UTR DATA")
print("=" * 70)

# Load BLAST hits
print("\nLoading deduplicated 3'UTR hits...")

hits = []
with open(blast_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 14:
            hits.append({
                'qseqid': parts[0],
                'te_id': parts[1],
                'pident': float(parts[2]),
                'length': int(parts[3]),
                'mismatch': int(parts[4]),
                'gapopen': int(parts[5]),
                'qstart': int(parts[6]),
                'qend': int(parts[7]),
                'sstart': int(parts[8]),
                'send': int(parts[9]),
                'evalue': float(parts[10]),
                'bitscore': float(parts[11]),
                'qlen': int(parts[12]) if parts[12].isdigit() else 0,
                'slen': int(parts[13]) if parts[13].isdigit() else 0,
            })

print(f"  Loaded {len(hits):,} deduplicated hits")

# Extract metrics
pidents = np.array([h['pident'] for h in hits])
lengths = np.array([h['length'] for h in hits])

# Parse TE families
def parse_te_family(te_id):
    """Extract TE family from FlyBase TE ID."""
    if '{' in te_id and '}' in te_id:
        return te_id.split('{')[1].split('}')[0]
    return te_id

te_families = [parse_te_family(h['te_id']) for h in hits]

# Aggregate counts
family_counts = defaultdict(int)
for fam in te_families:
    family_counts[fam] += 1

# ============================================================
# FIGURE 01: Identity Distribution
# ============================================================
print("\nGenerating Figure 01: Identity Distribution...")

fig, ax = plt.subplots(figsize=(10, 6))

bins = np.linspace(60, 100, 41)
n, bins_out, patches = ax.hist(pidents, bins=bins, color=COLORS['primary'],
                                edgecolor='white', linewidth=0.5, alpha=0.8)

# Add statistics
ax.axvline(np.median(pidents), color=COLORS['accent'], linestyle='--',
           linewidth=2, label=f'Median: {np.median(pidents):.1f}%')
ax.axvline(np.mean(pidents), color=COLORS['secondary'], linestyle=':',
           linewidth=2, label=f'Mean: {np.mean(pidents):.1f}%')

ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title("3'UTR TE Hit Identity Distribution\n(Deduplicated, n={:,})".format(len(hits)),
             fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '01_identity_distribution.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/01_identity_distribution.png")
plt.close()

# ============================================================
# FIGURE 02: Length Distribution
# ============================================================
print("\nGenerating Figure 02: Length Distribution...")

fig, ax = plt.subplots(figsize=(10, 6))

# Cap at 500bp for visualization
lengths_capped = np.minimum(lengths, 500)
bins = np.arange(0, 510, 10)
ax.hist(lengths_capped, bins=bins, color=COLORS['tertiary'],
        edgecolor='white', linewidth=0.5, alpha=0.8)

ax.axvline(np.median(lengths), color=COLORS['accent'], linestyle='--',
           linewidth=2, label=f'Median: {np.median(lengths):.0f}bp')
ax.axvline(np.mean(lengths), color=COLORS['secondary'], linestyle=':',
           linewidth=2, label=f'Mean: {np.mean(lengths):.0f}bp')

ax.set_xlabel('Alignment Length (bp)', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title("3'UTR TE Hit Length Distribution\n(Deduplicated, n={:,})".format(len(hits)),
             fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '02_length_distribution.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/02_length_distribution.png")
plt.close()

# ============================================================
# FIGURE 03: Cumulative by Threshold
# ============================================================
print("\nGenerating Figure 03: Cumulative by Threshold...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: By identity
ax = axes[0]
thresholds = np.arange(60, 101, 5)
counts_above = [np.sum(pidents >= t) for t in thresholds]

ax.fill_between(thresholds, counts_above, alpha=0.3, color=COLORS['primary'])
ax.plot(thresholds, counts_above, 'o-', color=COLORS['primary'], linewidth=2, markersize=8)
ax.set_xlabel('Minimum Identity (%)', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title('Hits Above Identity Threshold', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add annotations for key thresholds
for t, c in [(70, counts_above[2]), (80, counts_above[4]), (90, counts_above[6])]:
    ax.annotate(f'{c:,}', (t, c), textcoords='offset points',
                xytext=(0, 10), ha='center', fontsize=9)

# Right: By length
ax = axes[1]
len_thresholds = np.arange(20, 201, 20)
len_counts = [np.sum(lengths >= t) for t in len_thresholds]

ax.fill_between(len_thresholds, len_counts, alpha=0.3, color=COLORS['tertiary'])
ax.plot(len_thresholds, len_counts, 'o-', color=COLORS['tertiary'], linewidth=2, markersize=8)
ax.set_xlabel('Minimum Length (bp)', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title('Hits Above Length Threshold', fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add annotations
for t, c in [(50, len_counts[1]), (100, len_counts[4])]:
    ax.annotate(f'{c:,}', (t, c), textcoords='offset points',
                xytext=(0, 10), ha='center', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / '03_cumulative_by_threshold.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/03_cumulative_by_threshold.png")
plt.close()

# ============================================================
# FIGURE 04: TE Families Comparison
# ============================================================
print("\nGenerating Figure 04: TE Families Comparison...")

fig, ax = plt.subplots(figsize=(12, 8))

# Top 20 families
top_families = sorted(family_counts.items(), key=lambda x: -x[1])[:20]
families = [f[0] for f in top_families]
counts = [f[1] for f in top_families]

colors = plt.cm.tab20(np.linspace(0, 1, len(families)))
bars = ax.barh(range(len(families)), counts, color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(families)))
ax.set_yticklabels(families, fontsize=10)
ax.invert_yaxis()
ax.set_xlabel('Number of Hits (Deduplicated)', fontsize=12)
ax.set_title("Top 20 TE Families in 3'UTR Hits", fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add count labels
for bar, count in zip(bars, counts):
    ax.text(count + max(counts)*0.01, bar.get_y() + bar.get_height()/2,
            f'{count:,}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / '04_te_families_comparison.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/04_te_families_comparison.png")
plt.close()

# ============================================================
# FIGURE 05: High-Quality Hits by Family
# ============================================================
print("\nGenerating Figure 05: High-Quality Hits by Family...")

fig, ax = plt.subplots(figsize=(12, 8))

# Count high-quality hits (>=80% identity, >=50bp) by family
hq_family_counts = defaultdict(int)
for h, fam in zip(hits, te_families):
    if h['pident'] >= 80 and h['length'] >= 50:
        hq_family_counts[fam] += 1

# Calculate percent HQ
top_families_data = []
for fam, total in sorted(family_counts.items(), key=lambda x: -x[1])[:20]:
    hq = hq_family_counts.get(fam, 0)
    pct = 100 * hq / total if total > 0 else 0
    top_families_data.append((fam, total, hq, pct))

families = [f[0] for f in top_families_data]
pcts = [f[3] for f in top_families_data]

bars = ax.barh(range(len(families)), pcts, color=COLORS['secondary'], edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(families)))
ax.set_yticklabels(families, fontsize=10)
ax.invert_yaxis()
ax.set_xlabel('Percent High-Quality Hits (>=80% id, >=50bp)', fontsize=12)
ax.set_title("High-Quality Hit Fraction by TE Family", fontsize=13, fontweight='bold')
ax.set_xlim(0, 100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add labels
for bar, (fam, total, hq, pct) in zip(bars, top_families_data):
    ax.text(pct + 1, bar.get_y() + bar.get_height()/2,
            f'{pct:.1f}% ({hq:,}/{total:,})', va='center', fontsize=8)

plt.tight_layout()
plt.savefig(figures_dir / '05_percent_novel_by_family.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/05_percent_novel_by_family.png")
plt.close()

# ============================================================
# FIGURE 06: Identity vs Length 2D
# ============================================================
print("\nGenerating Figure 06: Identity vs Length 2D...")

fig, ax = plt.subplots(figsize=(10, 8))

# Sample for plotting if too many
if len(hits) > 50000:
    idx = np.random.choice(len(hits), 50000, replace=False)
    plot_pidents = pidents[idx]
    plot_lengths = lengths[idx]
else:
    plot_pidents = pidents
    plot_lengths = lengths

# 2D histogram
h = ax.hist2d(plot_pidents, np.minimum(plot_lengths, 300), bins=[40, 30],
               cmap='YlOrRd', cmin=1)
plt.colorbar(h[3], ax=ax, label='Number of Hits')

ax.set_xlabel('Percent Identity (%)', fontsize=12)
ax.set_ylabel('Alignment Length (bp)', fontsize=12)
ax.set_title("Identity vs Length Distribution\n(Deduplicated 3'UTR Hits)", fontsize=13, fontweight='bold')

# Add quality threshold lines
ax.axhline(50, color='black', linestyle='--', linewidth=1.5, alpha=0.7, label='50bp threshold')
ax.axvline(80, color='black', linestyle=':', linewidth=1.5, alpha=0.7, label='80% threshold')
ax.legend(loc='upper left', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / '06_identity_vs_length_2d.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/06_identity_vs_length_2d.png")
plt.close()

# ============================================================
# FIGURE 07: Overall Pie Chart
# ============================================================
print("\nGenerating Figure 07: Overall Pie Chart...")

fig, ax = plt.subplots(figsize=(10, 10))

# Categorize hits by quality
categories = {
    'High-quality\n(>=80% id, >=50bp)': sum(1 for h in hits if h['pident'] >= 80 and h['length'] >= 50),
    'Medium-quality\n(70-80% id or 30-50bp)': sum(1 for h in hits if (70 <= h['pident'] < 80 or 30 <= h['length'] < 50)
                                                   and not (h['pident'] >= 80 and h['length'] >= 50)),
    'Low-quality\n(<70% id and <30bp)': sum(1 for h in hits if h['pident'] < 70 and h['length'] < 30),
}
# Recategorize properly
hq = sum(1 for h in hits if h['pident'] >= 80 and h['length'] >= 50)
lq = sum(1 for h in hits if h['pident'] < 70 or h['length'] < 30)
mq = len(hits) - hq - lq
if mq < 0:
    mq = 0
    lq = len(hits) - hq

categories = {
    'High-quality\n(>=80% id, >=50bp)': hq,
    'Medium-quality': mq,
    'Low-quality\n(<70% id or <30bp)': lq,
}

colors = [COLORS['primary'], COLORS['secondary'], COLORS['neutral']]
explode = (0.05, 0.02, 0.02)

wedges, texts, autotexts = ax.pie(categories.values(), labels=categories.keys(),
                                   colors=colors, explode=explode,
                                   autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*len(hits)):,})',
                                   startangle=90, textprops={'fontsize': 11},
                                   wedgeprops=dict(edgecolor='black', linewidth=1))

ax.set_title("Hit Quality Distribution\n(Deduplicated 3'UTR Hits)", fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / '07_overall_pie_chart.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/07_overall_pie_chart.png")
plt.close()

# ============================================================
# FIGURE 08: Stacked by Identity
# ============================================================
print("\nGenerating Figure 08: Stacked by Identity...")

fig, ax = plt.subplots(figsize=(12, 7))

# Bin by identity
identity_bins = [(60, 65), (65, 70), (70, 75), (75, 80), (80, 85), (85, 90), (90, 95), (95, 100)]
bin_labels = [f'{lo}-{hi}%' for lo, hi in identity_bins]

# Count by length categories within each identity bin
length_cats = [('>=100bp', lambda l: l >= 100),
               ('50-100bp', lambda l: 50 <= l < 100),
               ('30-50bp', lambda l: 30 <= l < 50),
               ('<30bp', lambda l: l < 30)]

data = {cat: [] for cat, _ in length_cats}

for lo, hi in identity_bins:
    bin_hits = [h for h in hits if lo <= h['pident'] < hi]
    for cat, cond in length_cats:
        data[cat].append(sum(1 for h in bin_hits if cond(h['length'])))

# Stacked bar chart
x = np.arange(len(bin_labels))
width = 0.7

bottom = np.zeros(len(bin_labels))
colors_len = [COLORS['primary'], COLORS['tertiary'], COLORS['secondary'], COLORS['neutral']]

for (cat, _), color in zip(length_cats, colors_len):
    ax.bar(x, data[cat], width, bottom=bottom, label=cat, color=color, edgecolor='black', linewidth=0.3)
    bottom += np.array(data[cat])

ax.set_xticks(x)
ax.set_xticklabels(bin_labels, fontsize=10)
ax.set_xlabel('Identity Bin', fontsize=12)
ax.set_ylabel('Number of Hits', fontsize=12)
ax.set_title("Hit Distribution by Identity and Length\n(Deduplicated 3'UTR)", fontsize=13, fontweight='bold')
ax.legend(title='Alignment Length', fontsize=9, title_fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '08_stacked_by_identity.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/08_stacked_by_identity.png")
plt.close()

# ============================================================
# FIGURE 09: Top Transcripts
# ============================================================
print("\nGenerating Figure 09: Top Transcripts...")

fig, ax = plt.subplots(figsize=(12, 8))

# Count hits per transcript
transcript_hits = defaultdict(int)
transcript_hq = defaultdict(int)
for h in hits:
    transcript_hits[h['qseqid']] += 1
    if h['pident'] >= 80 and h['length'] >= 50:
        transcript_hq[h['qseqid']] += 1

# Top 30 transcripts by total hits
top_transcripts = sorted(transcript_hits.items(), key=lambda x: -x[1])[:30]

y_pos = np.arange(len(top_transcripts))
totals = [t[1] for t in top_transcripts]
hqs = [transcript_hq.get(t[0], 0) for t in top_transcripts]
names = [t[0] for t in top_transcripts]

ax.barh(y_pos, totals, color=COLORS['neutral'], edgecolor='black', linewidth=0.5, label='All hits')
ax.barh(y_pos, hqs, color=COLORS['primary'], edgecolor='black', linewidth=0.5, label='High-quality')

ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Number of Hits', fontsize=12)
ax.set_title("Top 30 Transcripts by TE Hits\n(Deduplicated 3'UTR)", fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / '09_top_transcripts_novel.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/09_top_transcripts_novel.png")
plt.close()

# ============================================================
# FIGURE 10: TE Class Breakdown
# ============================================================
print("\nGenerating Figure 10: TE Class Breakdown...")

# Define TE class mapping based on known family classifications
te_class_map = {
    'LTR': ['gypsy', 'copia', 'bel', 'pao', 'mdg1', 'mdg3', 'roo', '412', '297', 'blood',
            'Tirant', 'stalker', 'springer', 'Tabor', 'gtwin', 'idefix', 'ZAM', 'Dm88',
            'opus', 'flea', 'aurora', 'accord', 'HMS-Beagle', 'accord2', 'diver', 'invader'],
    'LINE': ['jockey', 'R1', 'R2', 'I-element', 'F-element', 'Doc', 'G-element', 'HeT-A',
             'TART', 'TAHRE', 'BS', 'Fw', 'Rt1a', 'Juan'],
    'DNA': ['hobo', 'P-element', 'mariner', 'pogo', 'S-element', 'Bari', 'transib', 'Tc1',
            'hAT', 'piggyBac', 'Helitron'],
    'Other': []  # Default category
}

# Classify families
family_class = {}
for fam in family_counts.keys():
    fam_lower = fam.lower()
    assigned = False
    for cls, members in te_class_map.items():
        for m in members:
            if m.lower() in fam_lower or fam_lower in m.lower():
                family_class[fam] = cls
                assigned = True
                break
        if assigned:
            break
    if not assigned:
        family_class[fam] = 'Other'

# Aggregate by class
class_counts = defaultdict(int)
for fam, count in family_counts.items():
    cls = family_class.get(fam, 'Other')
    class_counts[cls] += count

fig, ax = plt.subplots(figsize=(10, 8))

classes = ['LTR', 'LINE', 'DNA', 'Other']
counts = [class_counts[c] for c in classes]
colors = [COLORS['primary'], COLORS['secondary'], COLORS['tertiary'], COLORS['neutral']]

wedges, texts, autotexts = ax.pie(counts, labels=classes, colors=colors,
                                   autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(counts)):,})',
                                   startangle=90, textprops={'fontsize': 11},
                                   explode=(0.03, 0.03, 0.03, 0.03),
                                   wedgeprops=dict(edgecolor='black', linewidth=1))

ax.set_title("TE Class Distribution\n(Deduplicated 3'UTR Hits)", fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / '10_te_class_breakdown.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/10_te_class_breakdown.png")
plt.close()

# ============================================================
# FIGURES 11-14: Conservation Analysis
# ============================================================
print("\nNote: Figures 11-14 require conservation scores (phyloP).")
print("Checking for conservation data...")

# Load conservation scores from dedicated file
conservation_tab = Path('results/repeatmasker_analysis/te_hits_all_conservation.tab')
synteny_file = Path('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')

conservation_scores = {}
if conservation_tab.exists():
    print(f"  Loading conservation from {conservation_tab.name}...")
    with open(conservation_tab) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                name = parts[0]
                phyloP = float(parts[5]) if parts[5] not in ['NA', ''] else 0
                conservation_scores[name] = phyloP
    print(f"  Loaded {len(conservation_scores):,} conservation scores")

if synteny_file.exists() and conservation_scores:
    print("  Merging with synteny data, generating figures 11-14...")

    # Load and merge conservation with synteny
    cons_hits = []
    with open(synteny_file) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                name = parts[3]
                pident = float(parts[5])
                length = int(float(parts[6]))
                phyloP = conservation_scores.get(name, 0)

                cons_hits.append({
                    'pident': pident,
                    'length': length,
                    'phyloP': phyloP,
                })

    print(f"  Merged {len(cons_hits):,} hits with conservation data")

    phyloP_scores = np.array([h['phyloP'] for h in cons_hits])
    cons_pidents = np.array([h['pident'] for h in cons_hits])
    cons_lengths = np.array([h['length'] for h in cons_hits])

    # FIGURE 11: Conservation Distribution
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(phyloP_scores, bins=50, color=COLORS['primary'], edgecolor='white', alpha=0.8)
    ax.axvline(np.median(phyloP_scores), color=COLORS['accent'], linestyle='--',
               linewidth=2, label=f'Median: {np.median(phyloP_scores):.2f}')
    ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax.set_xlabel('phyloP Conservation Score', fontsize=12)
    ax.set_ylabel('Number of Hits', fontsize=12)
    ax.set_title('Conservation Distribution of TE Hits', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(figures_dir / '11_conservation_distribution.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {figures_dir}/11_conservation_distribution.png")
    plt.close()

    # FIGURE 12: Conservation vs Quality
    fig, ax = plt.subplots(figsize=(10, 8))

    # Sample for plotting
    if len(cons_hits) > 20000:
        idx = np.random.choice(len(cons_hits), 20000, replace=False)
        plot_cons = phyloP_scores[idx]
        plot_pident = cons_pidents[idx]
    else:
        plot_cons = phyloP_scores
        plot_pident = cons_pidents

    ax.scatter(plot_pident, plot_cons, alpha=0.3, s=5, c=COLORS['primary'])
    ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.axhline(1, color=COLORS['accent'], linestyle='--', linewidth=1.5, label='High conservation')
    ax.set_xlabel('Percent Identity (%)', fontsize=12)
    ax.set_ylabel('phyloP Conservation Score', fontsize=12)
    ax.set_title('Conservation vs Hit Quality', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(figures_dir / '12_conservation_vs_quality.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {figures_dir}/12_conservation_vs_quality.png")
    plt.close()

    # FIGURE 13: Conservation by Identity Bins
    fig, ax = plt.subplots(figsize=(10, 7))

    bin_data = []
    bin_labels = ['60-70%', '70-80%', '80-90%', '90-100%']
    for lo, hi in [(60, 70), (70, 80), (80, 90), (90, 100)]:
        subset = [phyloP_scores[i] for i in range(len(phyloP_scores))
                  if lo <= cons_pidents[i] < hi]
        bin_data.append(subset if subset else [0])

    bp = ax.boxplot(bin_data, labels=bin_labels, patch_artist=True)
    for patch, color in zip(bp['boxes'], [COLORS['neutral'], COLORS['secondary'],
                                           COLORS['tertiary'], COLORS['primary']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xlabel('Identity Bin', fontsize=12)
    ax.set_ylabel('phyloP Conservation Score', fontsize=12)
    ax.set_title('Conservation by Hit Identity', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(figures_dir / '13_conservation_by_family.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {figures_dir}/13_conservation_by_family.png")
    plt.close()

    # FIGURE 14: Conservation Categories Pie
    fig, ax = plt.subplots(figsize=(10, 10))

    # Filter out NaN values
    valid_scores = phyloP_scores[~np.isnan(phyloP_scores)]

    categories = {
        'Conserved\n(phyloP >= 1)': sum(1 for p in valid_scores if p >= 1),
        'Moderately conserved\n(0 <= phyloP < 1)': sum(1 for p in valid_scores if 0 <= p < 1),
        'Neutral/Divergent\n(phyloP < 0)': sum(1 for p in valid_scores if p < 0),
    }

    # Only create pie if there's data
    if sum(categories.values()) > 0:
        colors = [COLORS['primary'], COLORS['secondary'], COLORS['neutral']]

        wedges, texts, autotexts = ax.pie(categories.values(), labels=categories.keys(),
                                           colors=colors,
                                           autopct=lambda pct: f'{pct:.1f}%',
                                           startangle=90, textprops={'fontsize': 11},
                                           wedgeprops=dict(edgecolor='black', linewidth=1))

        ax.set_title('Conservation Categories\n(TE Hits)', fontsize=14, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No valid conservation data', ha='center', va='center', fontsize=14)
        ax.axis('off')

    plt.tight_layout()
    plt.savefig(figures_dir / '14_conservation_pie.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {figures_dir}/14_conservation_pie.png")
    plt.close()

else:
    print("  Conservation data not found. Skipping figures 11-14.")
    print("  Run conservation analysis first to generate these figures.")

print("\n" + "=" * 70)
print("DISTRIBUTION FIGURES COMPLETE")
print("=" * 70)
print(f"\nGenerated figures in {figures_dir}/")
