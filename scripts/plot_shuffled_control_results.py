#!/usr/bin/env python3
"""Visualize shuffled control analysis results."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Color scheme: contrasting colors for real vs shuffled
COLOR_REAL = '#2980b9'      # Bold blue
COLOR_SHUFFLED = '#e67e22'  # Orange

# Load results
results_dir = Path('results/shuffled_controls')

# Load replicate data
replicates = []
with open(results_dir / 'shuffled_replicates.tsv') as f:
    next(f)  # header
    for line in f:
        parts = line.strip().split('\t')
        replicates.append({
            'total_hits': int(parts[1]),
            'unique_queries': int(parts[2]),
            'mean_pident': float(parts[3]),
            'mean_length': float(parts[4]),
            'hits_hq': int(parts[5])
        })

# Load comparison
comparison = {}
with open(results_dir / 'control_comparison.tsv') as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        comparison[parts[0]] = {
            'real': float(parts[1]),
            'shuffled_mean': float(parts[2]),
            'shuffled_std': float(parts[3])
        }

# Real values
real_total = comparison['total_hits']['real']
real_hq = comparison['hits_hq']['real']
real_pident = comparison['mean_pident']['real']
real_length = comparison['mean_length']['real']

# Shuffled values
shuf_total = [r['total_hits'] for r in replicates]
shuf_hq = [r['hits_hq'] for r in replicates]
shuf_pident = [r['mean_pident'] for r in replicates]
shuf_length = [r['mean_length'] for r in replicates]

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Total hits comparison
ax = axes[0, 0]
x = ['Real\n3\'UTRs', 'Shuffled\n(10 reps)']
heights = [real_total, np.mean(shuf_total)]
errors = [0, np.std(shuf_total)]
colors = [COLOR_REAL, COLOR_SHUFFLED]

bars = ax.bar(x, heights, yerr=errors, capsize=5, color=colors, edgecolor='black')
ax.set_ylabel('Total BLAST Hits')
ax.set_title('Total TE Hits: Real vs Shuffled')

# Add fold change annotation
fold = real_total / np.mean(shuf_total)
ax.text(0.5, 0.95, f'{fold:.1f}x enrichment\nZ = {(real_total - np.mean(shuf_total))/np.std(shuf_total):.0f}',
        transform=ax.transAxes, ha='center', va='top', fontsize=12,
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

# Add values on bars
for bar, val in zip(bars, heights):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5000,
            f'{val:,.0f}', ha='center', va='bottom', fontsize=11)

# 2. High-quality hits comparison
ax = axes[0, 1]
heights = [real_hq, np.mean(shuf_hq)]
errors = [0, np.std(shuf_hq)]

bars = ax.bar(x, heights, yerr=errors, capsize=5, color=colors, edgecolor='black')
ax.set_ylabel('High-Quality Hits (≥80%, ≥50bp)')
ax.set_title('High-Quality Hits: Real vs Shuffled')

fold_hq = real_hq / np.mean(shuf_hq) if np.mean(shuf_hq) > 0 else float('inf')
ax.text(0.5, 0.95, f'{fold_hq:.0f}x enrichment',
        transform=ax.transAxes, ha='center', va='top', fontsize=12,
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

for bar, val in zip(bars, heights):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(heights)*0.02,
            f'{val:,.0f}', ha='center', va='bottom', fontsize=11)

# 3. Identity distribution comparison
ax = axes[1, 0]
ax.bar(['Real'], [real_pident], color=COLOR_REAL, edgecolor='black', label='Real')
ax.bar(['Shuffled'], [np.mean(shuf_pident)], yerr=[np.std(shuf_pident)],
       capsize=5, color=COLOR_SHUFFLED, edgecolor='black', label='Shuffled')
ax.set_ylabel('Mean Percent Identity (%)')
ax.set_title('Mean Hit Identity')
ax.set_ylim(70, 80)

# 4. Length distribution comparison
ax = axes[1, 1]
ax.bar(['Real'], [real_length], color=COLOR_REAL, edgecolor='black')
ax.bar(['Shuffled'], [np.mean(shuf_length)], yerr=[np.std(shuf_length)],
       capsize=5, color=COLOR_SHUFFLED, edgecolor='black')
ax.set_ylabel('Mean Alignment Length (bp)')
ax.set_title('Mean Hit Length')

# Add interpretation text
diff_length = real_length - np.mean(shuf_length)
ax.text(0.5, 0.95, f'Real hits {diff_length:.0f}bp longer',
        transform=ax.transAxes, ha='center', va='top', fontsize=11,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/26_shuffled_control_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/26_shuffled_control_comparison.png")

# Create second figure: replicate distribution
fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: Distribution of total hits
ax = axes[0]
ax.hist(shuf_total, bins=10, color=COLOR_SHUFFLED, edgecolor='black', alpha=0.7, label='Shuffled replicates')
ax.axvline(real_total, color=COLOR_REAL, linewidth=3, linestyle='--', label=f'Real ({real_total:,})')
ax.axvline(np.mean(shuf_total), color='#e74c3c', linewidth=2, linestyle=':', label=f'Shuffled mean ({np.mean(shuf_total):,.0f})')
ax.set_xlabel('Total BLAST Hits')
ax.set_ylabel('Frequency')
ax.set_title('Distribution of Total Hits')
ax.legend()

# Right: Distribution of HQ hits
ax = axes[1]
ax.hist(shuf_hq, bins=10, color=COLOR_SHUFFLED, edgecolor='black', alpha=0.7, label='Shuffled replicates')
ax.axvline(real_hq, color=COLOR_REAL, linewidth=3, linestyle='--', label=f'Real ({real_hq:,})')
ax.axvline(np.mean(shuf_hq), color='#e74c3c', linewidth=2, linestyle=':', label=f'Shuffled mean ({np.mean(shuf_hq):.0f})')
ax.set_xlabel('High-Quality Hits (≥80%, ≥50bp)')
ax.set_ylabel('Frequency')
ax.set_title('Distribution of High-Quality Hits')
ax.legend()

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/27_shuffled_replicate_distributions.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/27_shuffled_replicate_distributions.png")

# Print summary
print("\n" + "="*60)
print("SHUFFLED CONTROL SUMMARY")
print("="*60)
print(f"\nSample: 10% of 3'UTRs ({len(shuf_total)} sequences)")
print(f"Replicates: {len(replicates)}")
print(f"\nTotal hits:")
print(f"  Real:     {real_total:,}")
print(f"  Shuffled: {np.mean(shuf_total):,.0f} ± {np.std(shuf_total):,.0f}")
print(f"  Fold:     {fold:.1f}x")
print(f"  Z-score:  {(real_total - np.mean(shuf_total))/np.std(shuf_total):.1f}")
print(f"\nHigh-quality hits (≥80%, ≥50bp):")
print(f"  Real:     {real_hq:,}")
print(f"  Shuffled: {np.mean(shuf_hq):.0f} ± {np.std(shuf_hq):.0f}")
print(f"  Fold:     {fold_hq:.0f}x")
print(f"\nConclusion:")
print(f"  ~46% of total hits ({np.mean(shuf_total)/real_total*100:.0f}%) are due to sequence composition")
print(f"  ~54% ({(real_total - np.mean(shuf_total))/real_total*100:.0f}%) represent genuine TE-derived content")
print(f"  High-quality hits are {fold_hq:.0f}x enriched - virtually all are genuine")

plt.close('all')
