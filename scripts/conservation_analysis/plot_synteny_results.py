#!/usr/bin/env python3
"""Visualize synteny analysis results."""

import matplotlib.pyplot as plt
import numpy as np

def load_synteny_results(filepath):
    """Load synteny TSV."""
    results = []
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            for key in ['pident', 'length', 'sim_coverage', 'yak_coverage',
                       'ere_coverage', 'sec_coverage', 'any_coverage']:
                if key in row:
                    row[key] = float(row[key])
            results.append(row)
    return results

print("Loading data...")
hq = load_synteny_results('results/repeatmasker_analysis/te_hits_hq_synteny.tsv')
all_hits = load_synteny_results('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')

# Figure 1: Comparison of All vs HQ
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1A: Pie charts comparing All vs HQ
ax1 = axes[0, 0]
hq_syn = sum(1 for r in hq if r['any_coverage'] >= 0.5)
hq_mel = len(hq) - hq_syn
all_syn = sum(1 for r in all_hits if r['any_coverage'] >= 0.5)
all_mel = len(all_hits) - all_syn

x = np.arange(2)
width = 0.35

syn_pcts = [100*all_syn/len(all_hits), 100*hq_syn/len(hq)]
mel_pcts = [100*all_mel/len(all_hits), 100*hq_mel/len(hq)]

bars1 = ax1.bar(x - width/2, syn_pcts, width, label='Syntenic (conserved)', color='#27ae60')
bars2 = ax1.bar(x + width/2, mel_pcts, width, label='Mel-specific', color='#e74c3c')

ax1.set_ylabel('Percentage')
ax1.set_title('Synteny: All Hits vs High-Quality Hits')
ax1.set_xticks(x)
ax1.set_xticklabels([f'All Hits\n(n={len(all_hits):,})', f'HQ (≥80%/≥50bp)\n(n={len(hq):,})'])
ax1.legend()
ax1.set_ylim(0, 110)

# Add percentage labels
for bar in bars1:
    height = bar.get_height()
    ax1.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', fontweight='bold')
for bar in bars2:
    height = bar.get_height()
    ax1.annotate(f'{height:.1f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', fontweight='bold')

# 1B: Synteny by identity threshold (HQ hits)
ax2 = axes[0, 1]
id_ranges = [(80, 85), (85, 90), (90, 95), (95, 100)]
id_labels = ['80-85%', '85-90%', '90-95%', '95-100%']
id_syntenic = []

for min_id, max_id in id_ranges:
    tier = [r for r in hq if min_id <= r['pident'] < max_id]
    n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
    id_syntenic.append(100 * n_syn / len(tier) if tier else 0)

bars = ax2.bar(id_labels, id_syntenic, color='#3498db')
ax2.set_ylabel('% Syntenic (in any species)')
ax2.set_xlabel('Percent Identity')
ax2.set_title('HQ Hits: Higher Identity = LESS Syntenic')
ax2.set_ylim(0, 100)
ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='50% threshold')

for bar in bars:
    height = bar.get_height()
    ax2.annotate(f'{height:.0f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', fontweight='bold')

# 1C: Synteny by length (HQ hits)
ax3 = axes[1, 0]
len_ranges = [(50, 100), (100, 200), (200, 500), (500, 1000)]
len_labels = ['50-100bp', '100-200bp', '200-500bp', '500-1000bp']
len_syntenic = []

for min_len, max_len in len_ranges:
    tier = [r for r in hq if min_len <= r['length'] < max_len]
    n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
    len_syntenic.append(100 * n_syn / len(tier) if tier else 0)

bars = ax3.bar(len_labels, len_syntenic, color='#9b59b6')
ax3.set_ylabel('% Syntenic (in any species)')
ax3.set_xlabel('Alignment Length')
ax3.set_title('HQ Hits: Longer Hits = LESS Syntenic')
ax3.set_ylim(0, 100)
ax3.axhline(50, color='red', linestyle='--', alpha=0.5)

for bar in bars:
    height = bar.get_height()
    ax3.annotate(f'{height:.0f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', fontweight='bold')

# 1D: Novel vs Known (HQ hits)
ax4 = axes[1, 1]
novel_hq = [r for r in hq if r['category'] == 'novel']
known_hq = [r for r in hq if r['category'] == 'known']

novel_syn = 100 * sum(1 for r in novel_hq if r['any_coverage'] >= 0.5) / len(novel_hq)
known_syn = 100 * sum(1 for r in known_hq if r['any_coverage'] >= 0.5) / len(known_hq)

bars = ax4.bar(['Novel\n(not in RepeatMasker)', 'Known\n(in RepeatMasker)'],
               [novel_syn, known_syn], color=['#27ae60', '#3498db'])
ax4.set_ylabel('% Syntenic (in any species)')
ax4.set_title('HQ Hits: Novel Hits are MORE Syntenic')
ax4.set_ylim(0, 100)

for bar in bars:
    height = bar.get_height()
    ax4.annotate(f'{height:.0f}%', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords="offset points", ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/17_synteny_analysis.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/17_synteny_analysis.png")

# Figure 2: Heatmap combining conservation and synteny
fig2, ax = plt.subplots(figsize=(10, 8))

# Load conservation data for HQ hits
cons_data = {}
with open('results/repeatmasker_analysis/te_hits_hq_synteny.tsv') as f:
    header = next(f).strip().split('\t')
    for line in f:
        parts = line.strip().split('\t')
        row = dict(zip(header, parts))
        cons_data[row['name']] = float(row['any_coverage'])

# Match with conservation scores
cons_scores = {}
with open('results/repeatmasker_analysis/te_hits_conservation_hq.tab') as f:
    for line in f:
        parts = line.strip().split('\t')
        name = parts[0]
        score = float(parts[5])
        cons_scores[name] = score

# Create 2D histogram
synteny_vals = []
conservation_vals = []

for name in cons_data:
    # Try to find matching conservation score
    for full_name in cons_scores:
        if name.split('|')[1:3] == full_name.split('|')[1:3]:  # Match transcript and TE
            synteny_vals.append(cons_data[name])
            conservation_vals.append(cons_scores[full_name])
            break

print(f"Matched {len(synteny_vals)} hits with both synteny and conservation")

h = ax.hist2d(synteny_vals, conservation_vals, bins=30, cmap='YlOrRd')
ax.set_xlabel('Syntenic Coverage (fraction aligned in other species)')
ax.set_ylabel('Conservation Score (phyloP)')
ax.set_title('Conservation vs Synteny for HQ Hits')
plt.colorbar(h[3], ax=ax, label='Count')
ax.axhline(1, color='blue', linestyle='--', alpha=0.7, label='Conserved threshold')
ax.axvline(0.5, color='green', linestyle='--', alpha=0.7, label='Syntenic threshold')
ax.legend()

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/18_conservation_vs_synteny.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/18_conservation_vs_synteny.png")

plt.show()
