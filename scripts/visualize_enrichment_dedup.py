#!/usr/bin/env python3
"""
Visualize TE enrichment results using DEDUPLICATED data.

Creates publication-quality figures from integrated analysis results.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set up matplotlib for publication quality
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Color palette
COLORS = {
    'enriched': '#0077BB',
    'depleted': '#EE7733',
    'neutral': '#BBBBBB',
    'accent': '#009988',
    'highlight': '#CC3311'
}

# Directories
results_dir = Path('results/integrated_te_analysis')
figures_dir = Path('figures/enrichment')
figures_dir.mkdir(parents=True, exist_ok=True)

print("Loading deduplicated enrichment data...")

# Load gene data
gene_df = pd.read_csv(results_dir / 'gene_te_all_regions.tsv', sep='\t')
enrich_df = pd.read_csv(results_dir / 'functional_enrichment_results.tsv', sep='\t')

print(f"Loaded {len(gene_df)} genes, {len(enrich_df)} gene sets")

# ============================================================
# FIGURE 1: TE Density Heatmap by Region
# ============================================================
fig, ax = plt.subplots(figsize=(10, 8))

# Top 30 genes by total hits
top_genes = gene_df.nlargest(30, 'total_hits').copy()

# Create heatmap data
heatmap_data = top_genes[['symbol', 'utr3_hits', 'utr5_hits', 'exon_hits']].set_index('symbol')
heatmap_data.columns = ["3'UTR", "5'UTR", "Exon"]

# Log transform for better visualization
heatmap_log = np.log10(heatmap_data + 1)

# Plot
import matplotlib.colors as mcolors
cmap = plt.cm.YlOrRd

im = ax.imshow(heatmap_log.values, aspect='auto', cmap=cmap)
ax.set_xticks(range(len(heatmap_data.columns)))
ax.set_xticklabels(heatmap_data.columns, fontsize=11)
ax.set_yticks(range(len(heatmap_data)))
ax.set_yticklabels(heatmap_data.index, fontsize=9)

ax.set_title('TE Hit Distribution Across Transcript Regions\n(Top 30 Genes by Total Hits, Deduplicated)',
             fontsize=12, fontweight='bold')

# Colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('Log10(Hits + 1)', fontsize=10)

# Add actual values as text
for i in range(len(heatmap_data)):
    for j in range(len(heatmap_data.columns)):
        val = heatmap_data.iloc[i, j]
        text_color = 'white' if heatmap_log.iloc[i, j] > 2 else 'black'
        ax.text(j, i, f'{val:,.0f}', ha='center', va='center', fontsize=7, color=text_color)

plt.tight_layout()
plt.savefig(figures_dir / 'te_density_heatmap.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/te_density_heatmap.png")
plt.close()

# ============================================================
# FIGURE 2: Enrichment Forest Plot
# ============================================================
fig, ax = plt.subplots(figsize=(10, 10))

# Sort by fold enrichment
enrich_sorted = enrich_df.sort_values('fold_enrichment', ascending=True)

# Take top 15 enriched and top 15 depleted
top_enriched = enrich_sorted.tail(15)
top_depleted = enrich_sorted.head(15)
plot_df = pd.concat([top_depleted, top_enriched])

y_pos = np.arange(len(plot_df))
colors = [COLORS['enriched'] if x > 1 else COLORS['depleted'] for x in plot_df['fold_enrichment']]

# Plot bars
bars = ax.barh(y_pos, plot_df['fold_enrichment'], color=colors, edgecolor='black', linewidth=0.8, height=0.7)

# Reference line at 1
ax.axvline(1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7, zorder=0)

# Labels
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['set_name'], fontsize=9)
ax.set_xlabel('Fold Enrichment (vs Genome-Wide)', fontsize=11)
ax.set_title('Functional Gene Set TE Enrichment\n(Deduplicated Data, Integrated Analysis)',
             fontsize=12, fontweight='bold')

# Add fold values
for bar, val in zip(bars, plot_df['fold_enrichment']):
    x_pos = bar.get_width() + 0.02 if val > 1 else bar.get_width() - 0.02
    ha = 'left' if val > 1 else 'right'
    ax.text(x_pos, bar.get_y() + bar.get_height()/2, f'{val:.2f}x',
            va='center', ha=ha, fontsize=8)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLORS['enriched'], edgecolor='black', label='Enriched (>1x)'),
    Patch(facecolor=COLORS['depleted'], edgecolor='black', label='Depleted (<1x)')
]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig(figures_dir / 'te_enrichment_forest.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/te_enrichment_forest.png")
plt.close()

# ============================================================
# FIGURE 3: Volcano Plot (Effect Size vs Significance)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 8))

# Calculate -log10(p-value)
enrich_df['neg_log_p'] = -np.log10(enrich_df['p_value'].replace(0, 1e-10))
enrich_df['log_fold'] = np.log2(enrich_df['fold_enrichment'])

# Color by significance
sig_threshold = -np.log10(0.05)
colors = []
for _, row in enrich_df.iterrows():
    if row['neg_log_p'] > sig_threshold and row['log_fold'] > 0.5:
        colors.append(COLORS['enriched'])
    elif row['neg_log_p'] > sig_threshold and row['log_fold'] < -0.5:
        colors.append(COLORS['depleted'])
    else:
        colors.append(COLORS['neutral'])

ax.scatter(enrich_df['log_fold'], enrich_df['neg_log_p'], c=colors, s=50, alpha=0.7, edgecolors='black', linewidths=0.5)

# Reference lines
ax.axhline(sig_threshold, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.axvline(0, color='gray', linestyle='--', linewidth=1, alpha=0.7)
ax.axvline(0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax.axvline(-0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)

# Labels for significant gene sets
for _, row in enrich_df.iterrows():
    if row['neg_log_p'] > 2 and abs(row['log_fold']) > 0.8:
        ax.annotate(row['set_name'], (row['log_fold'], row['neg_log_p']),
                   fontsize=7, alpha=0.8,
                   xytext=(5, 5), textcoords='offset points')

ax.set_xlabel('Log2(Fold Enrichment)', fontsize=11)
ax.set_ylabel('-Log10(p-value)', fontsize=11)
ax.set_title('TE Enrichment Volcano Plot\n(Deduplicated Data)', fontsize=12, fontweight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add annotations
ax.text(0.95, 0.95, f'Enriched: {sum(1 for c in colors if c == COLORS["enriched"])}',
        transform=ax.transAxes, ha='right', va='top', fontsize=10, color=COLORS['enriched'])
ax.text(0.05, 0.95, f'Depleted: {sum(1 for c in colors if c == COLORS["depleted"])}',
        transform=ax.transAxes, ha='left', va='top', fontsize=10, color=COLORS['depleted'])

plt.tight_layout()
plt.savefig(figures_dir / 'te_enrichment_volcano.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/te_enrichment_volcano.png")
plt.close()

# ============================================================
# FIGURE 4: Top Sets Bar Chart by Category
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Parse category from set name
def get_category(name):
    if name.startswith('expr_'):
        return 'Expression'
    elif name.startswith('go_'):
        return 'GO Term'
    elif name.startswith('flyfish_'):
        return 'FlyFISH'
    else:
        return 'Other'

enrich_df['category'] = enrich_df['set_name'].apply(get_category)

# Left: Top enriched
ax = axes[0]
top_enrich = enrich_df.nlargest(15, 'fold_enrichment')
y_pos = np.arange(len(top_enrich))
colors = [COLORS['accent'] if c == 'FlyFISH' else COLORS['enriched'] for c in top_enrich['category']]

bars = ax.barh(y_pos, top_enrich['fold_enrichment'], color=colors, edgecolor='black', linewidth=0.8)
ax.axvline(1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(top_enrich['set_name'], fontsize=9)
ax.set_xlabel('Fold Enrichment', fontsize=11)
ax.set_title('Top 15 TE-Enriched Gene Sets', fontsize=12, fontweight='bold')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for bar, val in zip(bars, top_enrich['fold_enrichment']):
    ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2, f'{val:.2f}x',
            va='center', fontsize=8)

# Right: Top depleted
ax = axes[1]
top_deplete = enrich_df.nsmallest(15, 'fold_enrichment')
y_pos = np.arange(len(top_deplete))

bars = ax.barh(y_pos, top_deplete['fold_enrichment'], color=COLORS['depleted'], edgecolor='black', linewidth=0.8)
ax.axvline(1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(top_deplete['set_name'], fontsize=9)
ax.set_xlabel('Fold Enrichment', fontsize=11)
ax.set_title('Top 15 TE-Depleted Gene Sets', fontsize=12, fontweight='bold')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for bar, val in zip(bars, top_deplete['fold_enrichment']):
    ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.2f}x',
            va='center', fontsize=8)

plt.tight_layout()
plt.savefig(figures_dir / 'te_enrichment_by_category.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/te_enrichment_by_category.png")
plt.close()

# ============================================================
# FIGURE 5: Summary Overview
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('TE Enrichment Summary (Deduplicated Data)', fontsize=14, fontweight='bold', y=0.98)

# Panel A: Region distribution
ax = axes[0, 0]
regions = ["3'UTR", "5'UTR", "Exon"]
hits = [gene_df['utr3_hits'].sum(), gene_df['utr5_hits'].sum(), gene_df['exon_hits'].sum()]
colors_pie = [COLORS['enriched'], COLORS['depleted'], COLORS['accent']]

wedges, texts, autotexts = ax.pie(hits, labels=regions, colors=colors_pie,
                                   autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*sum(hits)):,})',
                                   startangle=90, explode=(0.02, 0.02, 0.02),
                                   wedgeprops=dict(edgecolor='black', linewidth=1))
ax.set_title('A. TE Hits by Region', fontsize=12, loc='left')

# Panel B: Genes by region combination
ax = axes[0, 1]
combos = {
    'All 3': len(gene_df[gene_df['regions_with_hits'] == 3]),
    '2 regions': len(gene_df[gene_df['regions_with_hits'] == 2]),
    '1 region': len(gene_df[gene_df['regions_with_hits'] == 1]),
}
ax.bar(combos.keys(), combos.values(), color=COLORS['enriched'], edgecolor='black')
ax.set_ylabel('Number of Genes', fontsize=11)
ax.set_title('B. Genes by Region Coverage', fontsize=12, loc='left')
for i, (k, v) in enumerate(combos.items()):
    ax.text(i, v + 100, f'{v:,}', ha='center', fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Enrichment distribution
ax = axes[1, 0]
ax.hist(enrich_df['fold_enrichment'], bins=20, color=COLORS['neutral'], edgecolor='black', alpha=0.7)
ax.axvline(1.0, color='black', linestyle='--', linewidth=2, label='No enrichment')
ax.axvline(enrich_df['fold_enrichment'].median(), color=COLORS['enriched'], linestyle='-',
           linewidth=2, label=f"Median: {enrich_df['fold_enrichment'].median():.2f}x")
ax.set_xlabel('Fold Enrichment', fontsize=11)
ax.set_ylabel('Number of Gene Sets', fontsize=11)
ax.set_title('C. Distribution of Enrichment Scores', fontsize=12, loc='left')
ax.legend(fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Category breakdown
ax = axes[1, 1]
cat_counts = enrich_df['category'].value_counts()
ax.bar(cat_counts.index, cat_counts.values, color=COLORS['accent'], edgecolor='black')
ax.set_ylabel('Number of Gene Sets', fontsize=11)
ax.set_title('D. Gene Sets by Category', fontsize=12, loc='left')
ax.set_xticklabels(cat_counts.index, rotation=30, ha='right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(figures_dir / 'te_enrichment_top_sets.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/te_enrichment_top_sets.png")
plt.close()

print("\n" + "="*60)
print("ENRICHMENT VISUALIZATION COMPLETE")
print("="*60)
print(f"Generated 5 figures in {figures_dir}/")
print(f"  - te_density_heatmap.png")
print(f"  - te_enrichment_forest.png")
print(f"  - te_enrichment_volcano.png")
print(f"  - te_enrichment_by_category.png")
print(f"  - te_enrichment_top_sets.png")
