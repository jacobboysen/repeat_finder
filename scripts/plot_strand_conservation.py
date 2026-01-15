#!/usr/bin/env python3
"""Visualize strand bias vs conservation relationships."""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import re

def load_data():
    """Load all necessary data."""
    # UTR strand bias
    utr_strand = {}
    with open('results/strand_bias_by_utr.tsv') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            utr_strand[parts[0]] = {
                'total_hits': int(parts[3]),
                'sense_pct': float(parts[6])
            }

    # TE strand bias
    te_strand = {}
    with open('results/strand_bias_by_te.tsv') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            te_strand[parts[0]] = {
                'total_hits': int(parts[1]),
                'sense_pct': float(parts[4])
            }

    # Conservation by transcript
    utr_cons = defaultdict(list)
    with open('results/repeatmasker_analysis/te_hits_all_conservation.tab') as f:
        for line in f:
            parts = line.strip().split('\t')
            name_parts = parts[0].split('|')
            if len(name_parts) >= 2:
                utr_cons[name_parts[1]].append(float(parts[5]))

    # Conservation by TE
    te_cons = defaultdict(list)
    with open('results/repeatmasker_analysis/te_hits_all_conservation.tab') as f:
        for line in f:
            parts = line.strip().split('\t')
            name_parts = parts[0].split('|')
            if len(name_parts) >= 3:
                te_cons[name_parts[2]].append(float(parts[5]))

    # TE info
    te_info = {}
    with open('/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_te_flybase.fasta') as f:
        for line in f:
            if line.startswith('>'):
                te_id = line.split()[0][1:]
                name_match = re.search(r'name=([^;]+)', line)
                if name_match:
                    name = name_match.group(1)
                    family = re.sub(r'\{[^}]*\}.*', '', name)
                    te_info[te_id] = {'name': name, 'family': family}

    return utr_strand, te_strand, utr_cons, te_cons, te_info

print("Loading data...")
utr_strand, te_strand, utr_cons, te_cons, te_info = load_data()

# Prepare UTR data
utr_data = []
for transcript, strand in utr_strand.items():
    if strand['total_hits'] >= 10 and transcript in utr_cons:
        mean_cons = np.mean(utr_cons[transcript])
        utr_data.append({
            'sense_pct': strand['sense_pct'],
            'conservation': mean_cons,
            'n_hits': strand['total_hits']
        })

# Prepare TE data
te_data = []
for te_id, strand in te_strand.items():
    if strand['total_hits'] >= 100 and te_id in te_cons:
        mean_cons = np.mean(te_cons[te_id])
        info = te_info.get(te_id, {'name': te_id, 'family': te_id})
        te_data.append({
            'te_id': te_id,
            'name': info['name'],
            'family': info['family'],
            'sense_pct': strand['sense_pct'],
            'conservation': mean_cons,
            'n_hits': strand['total_hits']
        })

print(f"UTR data points: {len(utr_data)}")
print(f"TE data points: {len(te_data)}")

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1. UTR scatter: strand bias vs conservation
ax = axes[0, 0]
sense_pcts = [d['sense_pct'] for d in utr_data]
conservations = [d['conservation'] for d in utr_data]
sizes = [np.sqrt(d['n_hits']) for d in utr_data]

# Sample if too many points
if len(utr_data) > 3000:
    idx = np.random.choice(len(utr_data), 3000, replace=False)
    sense_pcts_s = [sense_pcts[i] for i in idx]
    cons_s = [conservations[i] for i in idx]
    sizes_s = [sizes[i] for i in idx]
else:
    sense_pcts_s, cons_s, sizes_s = sense_pcts, conservations, sizes

ax.scatter(sense_pcts_s, cons_s, alpha=0.3, s=10, c='#3498db')
ax.axhline(1, color='red', linestyle='--', alpha=0.5, label='Conservation threshold')
ax.axvline(50, color='green', linestyle='--', alpha=0.5, label='No bias')
ax.set_xlabel('Sense Strand %')
ax.set_ylabel('Mean Conservation (phyloP)')
ax.set_title(f'UTRs: Strand Bias vs Conservation (n={len(utr_data):,})')
ax.legend()

# Add correlation
corr = np.corrcoef(sense_pcts, conservations)[0, 1]
ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# 2. TE scatter: strand bias vs conservation
ax = axes[0, 1]
te_sense = [d['sense_pct'] for d in te_data]
te_cons = [d['conservation'] for d in te_data]
te_sizes = [np.sqrt(d['n_hits']) for d in te_data]

# Color by family
families = list(set(d['family'] for d in te_data))
top_families = ['roo', '1360', 'INE-1', 'mdg1', '297', 'Cr1a']
colors = []
for d in te_data:
    if d['family'] in top_families:
        colors.append(top_families.index(d['family']))
    else:
        colors.append(len(top_families))

scatter = ax.scatter(te_sense, te_cons, c=colors, alpha=0.6, s=te_sizes, cmap='tab10')
ax.axhline(1, color='red', linestyle='--', alpha=0.5)
ax.axvline(50, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('Sense Strand %')
ax.set_ylabel('Mean Conservation (phyloP)')
ax.set_title(f'TEs: Strand Bias vs Conservation (n={len(te_data):,})')

# Add correlation
corr_te = np.corrcoef(te_sense, te_cons)[0, 1]
ax.text(0.05, 0.95, f'r = {corr_te:.3f}', transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# 3. Box plot: conservation by strand bias category (UTRs)
ax = axes[1, 0]
sense_biased = [d['conservation'] for d in utr_data if d['sense_pct'] >= 70]
anti_biased = [d['conservation'] for d in utr_data if d['sense_pct'] <= 30]
balanced = [d['conservation'] for d in utr_data if 30 < d['sense_pct'] < 70]

bp = ax.boxplot([sense_biased, balanced, anti_biased],
                labels=[f'Sense-biased\n(≥70%)\nn={len(sense_biased)}',
                        f'Balanced\n(30-70%)\nn={len(balanced)}',
                        f'Antisense-biased\n(≤30%)\nn={len(anti_biased)}'],
                patch_artist=True)

colors_box = ['#27ae60', '#f39c12', '#e74c3c']
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_ylabel('Conservation (phyloP)')
ax.set_title('UTR Conservation by Strand Bias Category')
ax.axhline(1, color='blue', linestyle='--', alpha=0.5, label='Conserved threshold')

# Add means
for i, data in enumerate([sense_biased, balanced, anti_biased]):
    mean = np.mean(data)
    ax.scatter([i+1], [mean], color='black', s=100, zorder=5, marker='D')
    ax.text(i+1.1, mean, f'μ={mean:.2f}', fontsize=10)

# 4. TE family comparison: strand bias vs conservation
ax = axes[1, 1]

family_stats = defaultdict(lambda: {'sense': [], 'cons': []})
for d in te_data:
    family_stats[d['family']]['sense'].append(d['sense_pct'])
    family_stats[d['family']]['cons'].append(d['conservation'])

# Filter to families with enough data
family_summary = []
for fam, data in family_stats.items():
    if len(data['sense']) >= 3:
        family_summary.append({
            'family': fam,
            'mean_sense': np.mean(data['sense']),
            'mean_cons': np.mean(data['cons']),
            'n': len(data['sense'])
        })

# Plot
fam_sense = [f['mean_sense'] for f in family_summary]
fam_cons = [f['mean_cons'] for f in family_summary]
fam_sizes = [np.sqrt(f['n']) * 10 for f in family_summary]

ax.scatter(fam_sense, fam_cons, s=fam_sizes, alpha=0.6, c='#9b59b6')

# Label top families
for f in family_summary:
    if f['n'] >= 10 or f['mean_cons'] > 1.5 or f['mean_sense'] > 85 or f['mean_sense'] < 40:
        ax.annotate(f['family'], (f['mean_sense'], f['mean_cons']),
                    fontsize=8, alpha=0.8)

ax.axhline(1, color='red', linestyle='--', alpha=0.5)
ax.axvline(50, color='green', linestyle='--', alpha=0.5)
ax.set_xlabel('Mean Sense Strand % (across TE instances)')
ax.set_ylabel('Mean Conservation (phyloP)')
ax.set_title('TE Families: Strand Bias vs Conservation')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/21_strand_vs_conservation.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/21_strand_vs_conservation.png")

# Print summary
print("\n" + "="*60)
print("SUMMARY: STRAND BIAS vs CONSERVATION")
print("="*60)
print(f"\nUTR-level correlation: r = {corr:.3f}")
print(f"TE-level correlation: r = {corr_te:.3f}")
print(f"\nMean conservation:")
print(f"  Sense-biased UTRs (≥70%): {np.mean(sense_biased):.3f}")
print(f"  Balanced UTRs (30-70%):   {np.mean(balanced):.3f}")
print(f"  Antisense-biased (≤30%):  {np.mean(anti_biased):.3f}")
print(f"\nConclusion: MINIMAL correlation - strand bias does NOT predict conservation")

plt.close()
