#!/usr/bin/env python3
"""
Visualize gene-level TE comparison using DEDUPLICATED data.

Creates publication-quality figures:
1. Gene ranking by TE hits
2. TE family breakdown by top genes
3. Region comparison (3'UTR vs 5'UTR vs Exon)
"""

import sys
import re
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
    'utr3': '#0077BB',     # Blue
    'utr5': '#EE7733',     # Orange
    'exon': '#009988',     # Teal
    'accent': '#CC3311',   # Red
    'neutral': '#BBBBBB',  # Gray
}

def load_transcript_to_gene():
    """Load FBtr -> FBgn mapping from GFF."""
    mapping = {}
    gff = Path('data/references/dmel-all-r6.66.gff')

    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'mRNA':
                id_m = re.search(r'ID=([^;]+)', parts[8])
                parent_m = re.search(r'Parent=([^;]+)', parts[8])
                if id_m and parent_m:
                    mapping[id_m.group(1)] = parent_m.group(1)

    return mapping

def parse_te_family(te_id):
    """Extract TE family from FlyBase TE ID."""
    # Format: FBti0019149{mdg1}
    if '{' in te_id and '}' in te_id:
        return te_id.split('{')[1].split('}')[0]
    return te_id

# Directories
results_dir = Path('results/integrated_te_analysis')
blast_dir = Path('results/3utr_deduplicated')
figures_dir = Path('figures/gene_comparison')
figures_dir.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("GENE COMPARISON VISUALIZATION (DEDUPLICATED)")
print("=" * 60)

# Load transcript-to-gene mapping
print("\nLoading transcript-gene mapping from GFF...")
fbtr_to_fbgn = load_transcript_to_gene()
print(f"  Loaded {len(fbtr_to_fbgn)} mappings")

# Load integrated gene data
print("\nLoading gene data...")
gene_df = pd.read_csv(results_dir / 'gene_te_all_regions.tsv', sep='\t')
print(f"  Loaded {len(gene_df)} genes with TE hits")

# Load deduplicated 3'UTR BLAST hits for TE family analysis
print("\nLoading 3'UTR deduplicated hits for family analysis...")
blast_file = blast_dir / '3utr_deduplicated_hits.tsv'

# Parse TE families from BLAST hits
gene_te_families = defaultdict(lambda: defaultdict(int))

if blast_file.exists():
    with open(blast_file) as f:
        # No header in raw BLAST output
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            qseqid = parts[0]  # FBtr ID
            sseqid = parts[1]  # TE ID with family

            # Map transcript to gene
            fbgn = fbtr_to_fbgn.get(qseqid)
            if not fbgn:
                continue

            # Extract TE family
            te_family = parse_te_family(sseqid)
            gene_te_families[fbgn][te_family] += 1

    print(f"  Parsed TE families for {len(gene_te_families)} genes")
else:
    print(f"  Warning: {blast_file} not found, skipping family analysis")

# ============================================================
# FIGURE 1: Gene Ranking by Total TE Hits
# ============================================================
print("\nGenerating Figure 1: Gene Ranking...")

fig, axes = plt.subplots(1, 2, figsize=(14, 8))

# Top 30 genes by total hits
top_genes = gene_df.nlargest(30, 'total_hits').copy()

# Left panel: Stacked bar by region
ax = axes[0]
y_pos = np.arange(len(top_genes))

# Plot stacked bars
bar_width = 0.7
bars_utr3 = ax.barh(y_pos, top_genes['utr3_hits'], bar_width,
                    label="3'UTR", color=COLORS['utr3'], edgecolor='black', linewidth=0.5)
bars_utr5 = ax.barh(y_pos, top_genes['utr5_hits'], bar_width,
                    left=top_genes['utr3_hits'], label="5'UTR",
                    color=COLORS['utr5'], edgecolor='black', linewidth=0.5)
bars_exon = ax.barh(y_pos, top_genes['exon_hits'], bar_width,
                    left=top_genes['utr3_hits'] + top_genes['utr5_hits'],
                    label="Exon", color=COLORS['exon'], edgecolor='black', linewidth=0.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(top_genes['symbol'], fontsize=9)
ax.set_xlabel('Number of TE Hits (Deduplicated)', fontsize=11)
ax.set_title('Top 30 Genes by Total TE Hits\n(Stacked by Region)', fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Right panel: Proportion by region
ax = axes[1]

# Calculate proportions
top_genes['utr3_pct'] = top_genes['utr3_hits'] / top_genes['total_hits'] * 100
top_genes['utr5_pct'] = top_genes['utr5_hits'] / top_genes['total_hits'] * 100
top_genes['exon_pct'] = top_genes['exon_hits'] / top_genes['total_hits'] * 100

bars_utr3 = ax.barh(y_pos, top_genes['utr3_pct'], bar_width,
                    label="3'UTR", color=COLORS['utr3'], edgecolor='black', linewidth=0.5)
bars_utr5 = ax.barh(y_pos, top_genes['utr5_pct'], bar_width,
                    left=top_genes['utr3_pct'], label="5'UTR",
                    color=COLORS['utr5'], edgecolor='black', linewidth=0.5)
bars_exon = ax.barh(y_pos, top_genes['exon_pct'], bar_width,
                    left=top_genes['utr3_pct'] + top_genes['utr5_pct'],
                    label="Exon", color=COLORS['exon'], edgecolor='black', linewidth=0.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(top_genes['symbol'], fontsize=9)
ax.set_xlabel('Percentage of TE Hits', fontsize=11)
ax.set_xlim(0, 100)
ax.set_title('TE Hit Distribution by Region\n(Proportional)', fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / 'gene_ranking.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/gene_ranking.png")
plt.close()

# ============================================================
# FIGURE 2: TE Family Breakdown
# ============================================================
print("\nGenerating Figure 2: TE Family Breakdown...")

if gene_te_families:
    fig, ax = plt.subplots(figsize=(12, 8))

    # Get top 20 genes
    top20_genes = gene_df.nlargest(20, 'total_hits')

    # Aggregate TE families across all genes
    all_families = defaultdict(int)
    for fbgn in top20_genes['fbgn']:
        for family, count in gene_te_families.get(fbgn, {}).items():
            all_families[family] += count

    # Get top 8 families
    top_families = sorted(all_families.keys(), key=lambda x: -all_families[x])[:8]

    # Build stacked bar data
    gene_names = top20_genes['symbol'].tolist()
    y_pos = np.arange(len(gene_names))

    # Color palette for TE families
    family_colors = plt.cm.Set2(np.linspace(0, 1, len(top_families) + 1))

    # Initialize cumulative array
    cumulative = np.zeros(len(gene_names))

    # Plot each family
    for i, family in enumerate(top_families):
        counts = []
        for fbgn in top20_genes['fbgn']:
            counts.append(gene_te_families.get(fbgn, {}).get(family, 0))
        counts = np.array(counts)

        ax.barh(y_pos, counts, left=cumulative, label=family[:15],
                color=family_colors[i], edgecolor='black', linewidth=0.3)
        cumulative += counts

    # Add "Other" category
    other_counts = []
    for fbgn in top20_genes['fbgn']:
        total = sum(gene_te_families.get(fbgn, {}).values())
        top_sum = sum(gene_te_families.get(fbgn, {}).get(f, 0) for f in top_families)
        other_counts.append(total - top_sum)
    other_counts = np.array(other_counts)

    ax.barh(y_pos, other_counts, left=cumulative, label='Other',
            color=COLORS['neutral'], edgecolor='black', linewidth=0.3)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(gene_names, fontsize=9)
    ax.set_xlabel('Number of 3\'UTR TE Hits (Deduplicated)', fontsize=11)
    ax.set_title('TE Family Breakdown - Top 20 Genes\n(Deduplicated 3\'UTR Hits)', fontsize=12, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9, title='TE Family')
    ax.invert_yaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(figures_dir / 'te_family_breakdown.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"  Saved: {figures_dir}/te_family_breakdown.png")
    plt.close()
else:
    print("  Skipped: No TE family data available")

# ============================================================
# FIGURE 3: Region Comparison (3'UTR vs 5'UTR vs Exon)
# ============================================================
print("\nGenerating Figure 3: Region Comparison...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('TE Hit Distribution Across Transcript Regions\n(Deduplicated Data)',
             fontsize=14, fontweight='bold', y=0.98)

# Panel A: Total hits by region
ax = axes[0, 0]
regions = ["3'UTR", "5'UTR", "Exon"]
totals = [gene_df['utr3_hits'].sum(), gene_df['utr5_hits'].sum(), gene_df['exon_hits'].sum()]
colors = [COLORS['utr3'], COLORS['utr5'], COLORS['exon']]

bars = ax.bar(regions, totals, color=colors, edgecolor='black', linewidth=1.2)
ax.set_ylabel('Total TE Hits', fontsize=11)
ax.set_title('A. Total Hits by Region', fontsize=12, loc='left')
for bar, val in zip(bars, totals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20000,
            f'{val:,.0f}', ha='center', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(totals) * 1.15)

# Panel B: Genes with hits in each region
ax = axes[0, 1]
genes_with_hits = [
    len(gene_df[gene_df['utr3_hits'] > 0]),
    len(gene_df[gene_df['utr5_hits'] > 0]),
    len(gene_df[gene_df['exon_hits'] > 0])
]

bars = ax.bar(regions, genes_with_hits, color=colors, edgecolor='black', linewidth=1.2)
ax.set_ylabel('Number of Genes', fontsize=11)
ax.set_title('B. Genes with TE Hits', fontsize=12, loc='left')
for bar, val in zip(bars, genes_with_hits):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
            f'{val:,}', ha='center', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(genes_with_hits) * 1.15)

# Panel C: Distribution of hits per gene
ax = axes[1, 0]
for region, col, color in [("3'UTR", 'utr3_hits', COLORS['utr3']),
                            ("5'UTR", 'utr5_hits', COLORS['utr5']),
                            ("Exon", 'exon_hits', COLORS['exon'])]:
    data = gene_df[gene_df[col] > 0][col]
    ax.hist(data, bins=50, alpha=0.6, label=region, color=color, edgecolor='black', linewidth=0.3)

ax.set_xlabel('Hits per Gene', fontsize=11)
ax.set_ylabel('Number of Genes', fontsize=11)
ax.set_title('C. Distribution of Hits per Gene', fontsize=12, loc='left')
ax.legend(fontsize=9)
ax.set_xlim(0, 500)  # Cap for readability
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Scatter of 3'UTR vs 5'UTR hits
ax = axes[1, 1]
# Filter genes with hits in both regions
both = gene_df[(gene_df['utr3_hits'] > 0) & (gene_df['utr5_hits'] > 0)]
scatter = ax.scatter(both['utr3_hits'], both['utr5_hits'],
                     c=both['exon_hits'], cmap='YlOrRd',
                     s=30, alpha=0.6, edgecolors='black', linewidths=0.3)
ax.plot([0, 5000], [0, 5000], 'k--', alpha=0.5, label='1:1 line')
ax.set_xlabel("3'UTR Hits", fontsize=11)
ax.set_ylabel("5'UTR Hits", fontsize=11)
ax.set_title("D. 3'UTR vs 5'UTR Hits\n(color = exon hits)", fontsize=12, loc='left')
ax.set_xlim(0, 5000)
ax.set_ylim(0, 5000)
plt.colorbar(scatter, ax=ax, label='Exon Hits', shrink=0.8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(figures_dir / 'length_vs_signal.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"  Saved: {figures_dir}/length_vs_signal.png")
plt.close()

print("\n" + "=" * 60)
print("GENE COMPARISON VISUALIZATION COMPLETE")
print("=" * 60)
print(f"\nGenerated 3 figures in {figures_dir}/")
print(f"  - gene_ranking.png")
print(f"  - te_family_breakdown.png")
print(f"  - length_vs_signal.png")
