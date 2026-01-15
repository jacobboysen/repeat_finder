#!/usr/bin/env python3
"""Visualize the ancient functional TE candidate set."""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import re

def load_candidates(filepath):
    """Load ancient TE candidates."""
    candidates = []
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            for key in ['pident', 'length', 'phyloP', 'syntenic_species', 'sim_cov', 'yak_cov', 'ere_cov']:
                if key in row:
                    row[key] = float(row[key])
            row['syntenic_species'] = int(row['syntenic_species'])
            candidates.append(row)
    return candidates

def load_te_info(te_fasta):
    """Load TE family names."""
    te_info = {}
    with open(te_fasta) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.split()
                te_id = parts[0][1:]
                name_match = re.search(r'name=([^;]+)', line)
                if name_match:
                    te_info[te_id] = name_match.group(1)
                else:
                    te_info[te_id] = te_id
    return te_info

print("Loading data...")
candidates = load_candidates('results/repeatmasker_analysis/ancient_te_candidates.tsv')
te_info = load_te_info('/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_te_flybase.fasta')

print(f"Loaded {len(candidates):,} ancient TE candidates")

# Create figure
fig = plt.figure(figsize=(16, 14))

# 1. Conservation vs Identity scatter
ax1 = fig.add_subplot(2, 2, 1)
pidents = [c['pident'] for c in candidates]
phyloPs = [c['phyloP'] for c in candidates]

# Sample for plotting if too many
if len(candidates) > 5000:
    idx = np.random.choice(len(candidates), 5000, replace=False)
    pidents_sample = [pidents[i] for i in idx]
    phyloPs_sample = [phyloPs[i] for i in idx]
else:
    pidents_sample = pidents
    phyloPs_sample = phyloPs

ax1.scatter(pidents_sample, phyloPs_sample, alpha=0.3, s=10, c='#3498db')
ax1.axhline(1, color='red', linestyle='--', label='Conservation threshold')
ax1.set_xlabel('Percent Identity')
ax1.set_ylabel('phyloP Conservation Score')
ax1.set_title(f'Ancient TE Candidates: Conservation vs Identity\n(n={len(candidates):,})')
ax1.legend()

# Add marginal histograms
ax1_histx = ax1.inset_axes([0, 1.02, 1, 0.15])
ax1_histy = ax1.inset_axes([1.02, 0, 0.15, 1])
ax1_histx.hist(pidents, bins=30, color='#3498db', alpha=0.7)
ax1_histx.set_xlim(ax1.get_xlim())
ax1_histx.axis('off')
ax1_histy.hist(phyloPs, bins=30, orientation='horizontal', color='#3498db', alpha=0.7)
ax1_histy.set_ylim(ax1.get_ylim())
ax1_histy.axis('off')

# 2. Identity distribution comparison (all vs ancient)
ax2 = fig.add_subplot(2, 2, 2)

# Load all hits for comparison
all_pidents = []
with open('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        all_pidents.append(float(parts[5]))

bins = np.linspace(60, 100, 30)
ax2.hist(all_pidents, bins=bins, alpha=0.5, label=f'All hits (n={len(all_pidents):,})', density=True, color='gray')
ax2.hist(pidents, bins=bins, alpha=0.7, label=f'Ancient candidates (n={len(candidates):,})', density=True, color='#27ae60')
ax2.set_xlabel('Percent Identity')
ax2.set_ylabel('Density')
ax2.set_title('Identity Distribution: All Hits vs Ancient Candidates')
ax2.legend()

# 3. TE Family distribution (top 15)
ax3 = fig.add_subplot(2, 2, 3)

te_family_counts = defaultdict(int)
for c in candidates:
    te_name = te_info.get(c['te_id'], c['te_id'])
    family = re.sub(r'\{[^}]*\}.*', '', te_name)
    te_family_counts[family] += 1

top_families = sorted(te_family_counts.items(), key=lambda x: -x[1])[:15]
families = [f[0] for f in top_families]
counts = [f[1] for f in top_families]

bars = ax3.barh(range(len(families)), counts, color='#9b59b6')
ax3.set_yticks(range(len(families)))
ax3.set_yticklabels(families)
ax3.invert_yaxis()
ax3.set_xlabel('Number of Ancient Candidates')
ax3.set_title('Top 15 TE Families in Ancient Candidates')

# Add count labels
for i, (bar, count) in enumerate(zip(bars, counts)):
    ax3.text(count + 100, i, f'{count:,}', va='center', fontsize=9)

# 4. Synteny species distribution
ax4 = fig.add_subplot(2, 2, 4)

syn_counts = defaultdict(int)
for c in candidates:
    syn_counts[c['syntenic_species']] += 1

species_nums = sorted(syn_counts.keys())
species_counts = [syn_counts[s] for s in species_nums]

colors = ['#e74c3c', '#f39c12', '#27ae60', '#2980b9'][:len(species_nums)]
bars = ax4.bar(species_nums, species_counts, color=colors)
ax4.set_xlabel('Number of Species with Synteny')
ax4.set_ylabel('Number of Candidates')
ax4.set_title('Synteny Conservation Across Species')
ax4.set_xticks(species_nums)
ax4.set_xticklabels([f'{s} species' for s in species_nums])

# Add percentage labels
total = sum(species_counts)
for bar, count in zip(bars, species_counts):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height,
             f'{count:,}\n({100*count/total:.0f}%)',
             ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/19_ancient_candidates_overview.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/19_ancient_candidates_overview.png")

# Figure 2: Conservation and length details
fig2, axes = plt.subplots(2, 2, figsize=(14, 12))

# 2A: Conservation histogram
ax = axes[0, 0]
ax.hist(phyloPs, bins=50, color='#27ae60', edgecolor='white', alpha=0.8)
ax.axvline(np.median(phyloPs), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(phyloPs):.2f}')
ax.axvline(np.mean(phyloPs), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(phyloPs):.2f}')
ax.set_xlabel('phyloP Conservation Score')
ax.set_ylabel('Count')
ax.set_title('Conservation Distribution of Ancient Candidates')
ax.legend()

# 2B: Length histogram
ax = axes[0, 1]
lengths = [c['length'] for c in candidates]
ax.hist(lengths, bins=50, color='#3498db', edgecolor='white', alpha=0.8)
ax.axvline(np.median(lengths), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(lengths):.0f}bp')
ax.set_xlabel('Alignment Length (bp)')
ax.set_ylabel('Count')
ax.set_title('Length Distribution of Ancient Candidates')
ax.legend()

# 2C: Conservation by length bins
ax = axes[1, 0]
len_bins = [(0, 30), (30, 50), (50, 100), (100, 200), (200, 500)]
bin_labels = ['<30bp', '30-50bp', '50-100bp', '100-200bp', '200-500bp']
bin_phyloPs = []

for lo, hi in len_bins:
    subset = [c['phyloP'] for c in candidates if lo <= c['length'] < hi]
    bin_phyloPs.append(subset if subset else [0])

bp = ax.boxplot(bin_phyloPs, labels=bin_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('#3498db')
    patch.set_alpha(0.7)
ax.set_xlabel('Alignment Length')
ax.set_ylabel('phyloP Conservation Score')
ax.set_title('Conservation by Alignment Length')

# 2D: Category breakdown (novel vs known)
ax = axes[1, 1]
novel = [c for c in candidates if c['category'] == 'novel']
known = [c for c in candidates if c['category'] == 'known']

categories = ['Novel\n(not in RepeatMasker)', 'Known\n(in RepeatMasker)']
cat_counts = [len(novel), len(known)]
colors = ['#27ae60', '#3498db']

bars = ax.bar(categories, cat_counts, color=colors)
ax.set_ylabel('Number of Ancient Candidates')
ax.set_title('Ancient Candidates: Novel vs Known')

for bar, count in zip(bars, cat_counts):
    height = bar.get_height()
    pct = 100 * count / len(candidates)
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{count:,}\n({pct:.1f}%)',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/20_ancient_candidates_details.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/20_ancient_candidates_details.png")

# Summary stats
print("\n" + "="*70)
print("ANCIENT TE CANDIDATE STATISTICS")
print("="*70)
print(f"\nTotal candidates: {len(candidates):,}")
print(f"\nConservation (phyloP):")
print(f"  Mean: {np.mean(phyloPs):.2f}")
print(f"  Median: {np.median(phyloPs):.2f}")
print(f"  Range: {min(phyloPs):.2f} - {max(phyloPs):.2f}")
print(f"\nAlignment length:")
print(f"  Mean: {np.mean(lengths):.0f}bp")
print(f"  Median: {np.median(lengths):.0f}bp")
print(f"  Range: {min(lengths)} - {max(lengths)}bp")
print(f"\nIdentity:")
print(f"  Mean: {np.mean(pidents):.1f}%")
print(f"  Median: {np.median(pidents):.1f}%")
print(f"\nNovel vs Known:")
print(f"  Novel (not in RM): {len(novel):,} ({100*len(novel)/len(candidates):.1f}%)")
print(f"  Known (in RM): {len(known):,} ({100*len(known)/len(candidates):.1f}%)")

plt.show()
