#!/usr/bin/env python3
"""
FlyFISH localization vs TE content analysis using DEDUPLICATED data.

Analyzes relationship between mRNA subcellular localization and TE content.
"""

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir, get_references_dir

# Style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Color palette
COLORS = ['#0077BB', '#33BBEE', '#009988', '#EE7733', '#CC3311', '#EE3377', '#BBBBBB']

# Directories
results_dir = get_results_dir()
figures_dir = Path('figures/flyfish')
figures_dir.mkdir(parents=True, exist_ok=True)

def load_flyfish_gene_sets():
    """Load FlyFISH localization gene sets."""
    gene_sets = {}
    gene_sets_dir = Path('data/gene_lists/functional')

    for f in gene_sets_dir.glob('flyfish_*_fbgn_ids.txt'):
        name = f.stem.replace('_fbgn_ids', '')
        genes = set()
        with open(f) as fh:
            for line in fh:
                gene = line.strip()
                if gene.startswith('FBgn'):
                    genes.add(gene)
        if genes:
            gene_sets[name] = genes

    return gene_sets

def load_transcript_to_gene():
    """Load FBtr -> FBgn mapping."""
    import re
    mapping = {}
    gff = get_references_dir() / 'dmel-all-r6.66.gff'

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

def analyze_dedup_hits_by_localization(hits_file, gene_sets, fbtr_to_fbgn):
    """Analyze deduplicated hits by localization category."""
    results = defaultdict(lambda: {
        'n_genes': 0,
        'genes_with_hits': 0,
        'total_hits': 0,
        'sense_hits': 0,
        'anti_hits': 0,
        'hit_lengths': [],
        'hit_pidents': [],
        'te_families': defaultdict(int)
    })

    # Count hits per gene
    gene_hits = defaultdict(lambda: {'hits': 0, 'sense': 0, 'anti': 0, 'lengths': [], 'pidents': [], 'tes': defaultdict(int)})

    print(f"  Reading {hits_file}...")
    with open(hits_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            fbtr = parts[0].split('_')[0]  # Handle shuffled format
            fbgn = fbtr_to_fbgn.get(fbtr)
            if not fbgn:
                continue

            te_id = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            sstart = int(parts[8])
            send = int(parts[9])
            strand = '+' if sstart < send else '-'

            gene_hits[fbgn]['hits'] += 1
            gene_hits[fbgn]['lengths'].append(length)
            gene_hits[fbgn]['pidents'].append(pident)
            gene_hits[fbgn]['tes'][te_id] += 1
            if strand == '+':
                gene_hits[fbgn]['sense'] += 1
            else:
                gene_hits[fbgn]['anti'] += 1

    # Aggregate by localization
    for loc_name, genes in gene_sets.items():
        results[loc_name]['n_genes'] = len(genes)

        for gene in genes:
            if gene in gene_hits:
                gh = gene_hits[gene]
                results[loc_name]['genes_with_hits'] += 1
                results[loc_name]['total_hits'] += gh['hits']
                results[loc_name]['sense_hits'] += gh['sense']
                results[loc_name]['anti_hits'] += gh['anti']
                results[loc_name]['hit_lengths'].extend(gh['lengths'])
                results[loc_name]['hit_pidents'].extend(gh['pidents'])
                for te, count in gh['tes'].items():
                    results[loc_name]['te_families'][te] += count

    return dict(results)

print("="*70)
print("FlyFISH TE Analysis (DEDUPLICATED DATA)")
print("="*70)

# Load data
print("\nLoading reference data...")
gene_sets = load_flyfish_gene_sets()
print(f"  Loaded {len(gene_sets)} FlyFISH localization categories")

fbtr_to_fbgn = load_transcript_to_gene()
print(f"  Loaded {len(fbtr_to_fbgn)} transcript-gene mappings")

# Analyze 3'UTR deduplicated hits
print("\nAnalyzing 3'UTR deduplicated hits...")
utr3_hits_file = results_dir / '3utr_deduplicated' / '3utr_deduplicated_hits.tsv'
utr3_results = analyze_dedup_hits_by_localization(utr3_hits_file, gene_sets, fbtr_to_fbgn)

# Create summary dataframe
summary_data = []
for loc, data in utr3_results.items():
    if data['n_genes'] == 0:
        continue

    pct_with_hits = 100 * data['genes_with_hits'] / data['n_genes']
    sense_pct = 100 * data['sense_hits'] / data['total_hits'] if data['total_hits'] > 0 else 50
    mean_length = np.mean(data['hit_lengths']) if data['hit_lengths'] else 0
    mean_pident = np.mean(data['hit_pidents']) if data['hit_pidents'] else 0
    hits_per_gene = data['total_hits'] / data['n_genes'] if data['n_genes'] > 0 else 0

    summary_data.append({
        'localization': loc.replace('flyfish_', ''),
        'n_genes': data['n_genes'],
        'genes_with_hits': data['genes_with_hits'],
        'pct_with_hits': pct_with_hits,
        'total_hits': data['total_hits'],
        'hits_per_gene': hits_per_gene,
        'sense_pct': sense_pct,
        'mean_length': mean_length,
        'mean_pident': mean_pident
    })

summary_df = pd.DataFrame(summary_data).sort_values('total_hits', ascending=False)

print("\nFlyFISH Summary (Deduplicated 3'UTR):")
print(summary_df.to_string(index=False))

# ============================================================
# FIGURE 1: TE Presence by Localization
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

plot_df = summary_df.sort_values('pct_with_hits', ascending=True)
y_pos = np.arange(len(plot_df))

bars = ax.barh(y_pos, plot_df['pct_with_hits'], color=COLORS[0], edgecolor='black', linewidth=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['localization'], fontsize=10)
ax.set_xlabel('% Genes with TE Hits', fontsize=11)
ax.set_title('TE Content by mRNA Localization\n(Deduplicated 3\'UTR Data)', fontsize=13, fontweight='bold')
ax.set_xlim(0, 105)

for bar, val, n in zip(bars, plot_df['pct_with_hits'], plot_df['n_genes']):
    ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
            f'{val:.1f}% (n={n})', va='center', fontsize=9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / 'te_presence_by_localization.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"\nSaved: {figures_dir}/te_presence_by_localization.png")
plt.close()

# ============================================================
# FIGURE 2: Strand Bias by Localization
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

plot_df = summary_df.sort_values('sense_pct', ascending=True)
y_pos = np.arange(len(plot_df))

# Stacked bar chart
sense_vals = plot_df['sense_pct']
anti_vals = 100 - plot_df['sense_pct']

bars1 = ax.barh(y_pos, sense_vals, color=COLORS[0], edgecolor='black', linewidth=0.5, label='Sense')
bars2 = ax.barh(y_pos, anti_vals, left=sense_vals, color=COLORS[3], edgecolor='black', linewidth=0.5, label='Antisense')

ax.axvline(50, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['localization'], fontsize=10)
ax.set_xlabel('Strand Bias (%)', fontsize=11)
ax.set_title('TE Strand Orientation by mRNA Localization\n(Deduplicated 3\'UTR Data)', fontsize=13, fontweight='bold')
ax.legend(loc='lower right')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / 'strand_bias_by_localization.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/strand_bias_by_localization.png")
plt.close()

# ============================================================
# FIGURE 3: TE Family Heatmap by Localization
# ============================================================
# Build TE family matrix
te_family_data = []
for loc, data in utr3_results.items():
    loc_clean = loc.replace('flyfish_', '')
    for te, count in data['te_families'].items():
        te_family_data.append({'localization': loc_clean, 'te_family': te, 'count': count})

te_df = pd.DataFrame(te_family_data)
if len(te_df) > 0:
    # Get top 15 TE families overall
    top_tes = te_df.groupby('te_family')['count'].sum().nlargest(15).index.tolist()
    te_subset = te_df[te_df['te_family'].isin(top_tes)]

    pivot = te_subset.pivot_table(values='count', index='localization', columns='te_family', fill_value=0)

    fig, ax = plt.subplots(figsize=(14, 8))

    # Log transform for better visualization
    pivot_log = np.log10(pivot + 1)

    im = ax.imshow(pivot_log.values, aspect='auto', cmap='YlOrRd')
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=10)

    ax.set_title('TE Family Distribution by mRNA Localization\n(Top 15 Families, Deduplicated 3\'UTR)',
                 fontsize=13, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Log10(Hits + 1)', fontsize=10)

    plt.tight_layout()
    plt.savefig(figures_dir / 'te_family_heatmap.png', dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"Saved: {figures_dir}/te_family_heatmap.png")
    plt.close()

# ============================================================
# FIGURE 4: Maternal vs Non-Maternal Comparison
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Maternal vs Degraded mRNA: TE Content Comparison\n(Deduplicated Data)',
             fontsize=13, fontweight='bold', y=1.02)

# Get maternal and degraded data
maternal_data = summary_df[summary_df['localization'].str.contains('maternal', case=False)]
degraded_data = summary_df[summary_df['localization'].str.contains('degraded', case=False)]

# Panel A: Hits per gene
ax = axes[0]
if len(maternal_data) > 0 and len(degraded_data) > 0:
    categories = ['Maternal\nmRNAs', 'Degraded\nmRNAs']
    maternal_hpg = maternal_data['hits_per_gene'].values[0] if len(maternal_data) > 0 else 0
    degraded_hpg = degraded_data['hits_per_gene'].values[0] if len(degraded_data) > 0 else 0
    values = [maternal_hpg, degraded_hpg]

    bars = ax.bar(categories, values, color=[COLORS[0], COLORS[3]], edgecolor='black', width=0.5)
    ax.set_ylabel('Hits per Gene', fontsize=11)
    ax.set_title('A. TE Hit Density', fontsize=12, loc='left')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f'{val:.1f}', ha='center', fontsize=11, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Panel B: Strand bias
ax = axes[1]
if len(maternal_data) > 0 and len(degraded_data) > 0:
    maternal_sense = maternal_data['sense_pct'].values[0] if len(maternal_data) > 0 else 50
    degraded_sense = degraded_data['sense_pct'].values[0] if len(degraded_data) > 0 else 50
    values = [maternal_sense, degraded_sense]

    bars = ax.bar(categories, values, color=[COLORS[0], COLORS[3]], edgecolor='black', width=0.5)
    ax.axhline(50, color='black', linestyle='--', linewidth=1, alpha=0.7)
    ax.set_ylabel('% Sense Strand', fontsize=11)
    ax.set_title('B. Strand Bias', fontsize=12, loc='left')
    ax.set_ylim(0, 100)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                f'{val:.1f}%', ha='center', fontsize=11, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(figures_dir / 'degraded_vs_stable_maternal.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/degraded_vs_stable_maternal.png")
plt.close()

# ============================================================
# FIGURE 5: Summary Overview
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('FlyFISH Localization vs TE Content Summary\n(Deduplicated 3\'UTR Data)',
             fontsize=14, fontweight='bold', y=0.98)

# Panel A: Hits per gene by localization
ax = axes[0, 0]
plot_df = summary_df.nlargest(10, 'hits_per_gene')
y_pos = np.arange(len(plot_df))
bars = ax.barh(y_pos, plot_df['hits_per_gene'], color=COLORS[0], edgecolor='black')
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['localization'], fontsize=9)
ax.set_xlabel('Hits per Gene', fontsize=10)
ax.set_title('A. Top 10 by TE Hit Density', fontsize=11, loc='left')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel B: Total genes by category
ax = axes[0, 1]
plot_df = summary_df.nlargest(10, 'n_genes')
y_pos = np.arange(len(plot_df))
bars = ax.barh(y_pos, plot_df['n_genes'], color=COLORS[2], edgecolor='black')
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['localization'], fontsize=9)
ax.set_xlabel('Number of Genes', fontsize=10)
ax.set_title('B. Largest Gene Sets', fontsize=11, loc='left')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel C: Sense bias scatter
ax = axes[1, 0]
ax.scatter(summary_df['sense_pct'], summary_df['pct_with_hits'],
           s=summary_df['n_genes']/10, c=COLORS[0], alpha=0.7, edgecolors='black')
ax.axvline(50, color='gray', linestyle='--', linewidth=1)
ax.set_xlabel('% Sense Strand', fontsize=10)
ax.set_ylabel('% Genes with TE Hits', fontsize=10)
ax.set_title('C. Strand Bias vs TE Prevalence', fontsize=11, loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel D: Mean hit length by category
ax = axes[1, 1]
plot_df = summary_df.nlargest(10, 'mean_length')
y_pos = np.arange(len(plot_df))
bars = ax.barh(y_pos, plot_df['mean_length'], color=COLORS[4], edgecolor='black')
ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df['localization'], fontsize=9)
ax.set_xlabel('Mean Alignment Length (bp)', fontsize=10)
ax.set_title('D. Top 10 by Hit Length', fontsize=11, loc='left')
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(figures_dir / 'flyfish_summary_overview.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print(f"Saved: {figures_dir}/flyfish_summary_overview.png")
plt.close()

# Save summary table
summary_df.to_csv(results_dir / 'flyfish_hit_characteristics_dedup.tsv', sep='\t', index=False)
print(f"\nSaved: {results_dir}/flyfish_hit_characteristics_dedup.tsv")

print("\n" + "="*70)
print("FlyFISH Analysis Complete (Deduplicated)")
print("="*70)
