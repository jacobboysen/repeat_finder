#!/usr/bin/env python3
"""
Compare TE family distributions between real and shuffled sequences.
Analyze e-value distributions and quality thresholds.
"""

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from pathlib import Path
import re

def load_te_info(te_fasta):
    """Load TE family names from FASTA."""
    te_info = {}
    with open(te_fasta) as f:
        for line in f:
            if line.startswith('>'):
                te_id = line.split()[0][1:]
                name_match = re.search(r'name=([^;]+)', line)
                if name_match:
                    name = name_match.group(1)
                    # Extract family (remove instance suffix)
                    family = re.sub(r'\{[^}]*\}.*', '', name)
                    te_info[te_id] = {'name': name, 'family': family}
                else:
                    te_info[te_id] = {'name': te_id, 'family': te_id}
    return te_info

def parse_blast(filepath):
    """Parse BLAST results."""
    hits = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                hits.append({
                    'qseqid': parts[0],
                    'sseqid': parts[1],
                    'pident': float(parts[2]),
                    'length': int(parts[3]),
                    'evalue': float(parts[10]),
                    'bitscore': float(parts[11])
                })
    return hits

def get_family_counts(hits, te_info):
    """Count hits per TE family."""
    counts = defaultdict(int)
    for hit in hits:
        te_id = hit['sseqid']
        family = te_info.get(te_id, {}).get('family', te_id)
        counts[family] += 1
    return dict(counts)

def filter_hits(hits, min_pident=0, min_length=0):
    """Filter hits by identity and length."""
    return [h for h in hits if h['pident'] >= min_pident and h['length'] >= min_length]

# Load data
print("Loading data...")
te_info = load_te_info('data/references/dmel_te_flybase.fasta')
print(f"  Loaded {len(te_info)} TE annotations")

results_dir = Path('results/shuffled_controls')
real_hits = parse_blast(results_dir / 'real_blast.tsv')
print(f"  Real hits: {len(real_hits):,}")

# Load all shuffled replicates
shuf_hits_all = []
for i in range(1, 11):
    shuf_file = results_dir / f'shuffled_rep{i}_blast.tsv'
    if shuf_file.exists():
        shuf_hits_all.extend(parse_blast(shuf_file))
print(f"  Shuffled hits (10 reps): {len(shuf_hits_all):,}")

# Average per replicate for fair comparison
shuf_hits_avg_count = len(shuf_hits_all) / 10

# =====================================================================
# 1. TE FAMILY COMPARISON
# =====================================================================
print("\n" + "="*70)
print("TE FAMILY COMPARISON: REAL vs SHUFFLED")
print("="*70)

real_families = get_family_counts(real_hits, te_info)
shuf_families = get_family_counts(shuf_hits_all, te_info)

# Normalize shuffled to per-replicate
shuf_families_norm = {k: v/10 for k, v in shuf_families.items()}

# Get top families by real count
top_families = sorted(real_families.items(), key=lambda x: -x[1])[:20]

print(f"\n{'TE Family':<20} {'Real':>12} {'Shuffled':>12} {'Fold':>8} {'Enriched':>10}")
print("-"*65)

enrichment_data = []
for family, real_count in top_families:
    shuf_count = shuf_families_norm.get(family, 0)
    fold = real_count / shuf_count if shuf_count > 0 else float('inf')
    enriched = "YES" if fold > 1.5 else ("depleted" if fold < 0.67 else "similar")
    print(f"{family:<20} {real_count:>12,} {shuf_count:>12,.0f} {fold:>8.1f}x {enriched:>10}")
    enrichment_data.append({
        'family': family,
        'real': real_count,
        'shuffled': shuf_count,
        'fold': fold
    })

# =====================================================================
# 2. E-VALUE DISTRIBUTION
# =====================================================================
print("\n" + "="*70)
print("E-VALUE DISTRIBUTION")
print("="*70)

real_evalues = [h['evalue'] for h in real_hits]
shuf_evalues = [h['evalue'] for h in shuf_hits_all]

# E-value bins
evalue_bins = [0, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 1, 10]
bin_labels = ['<1e-10', '1e-10 to 1e-5', '1e-5 to 1e-3', '1e-3 to 0.01', '0.01 to 0.1', '0.1 to 1', '1 to 10']

def bin_evalues(evalues, bins):
    counts = []
    for i in range(len(bins)-1):
        count = sum(1 for e in evalues if bins[i] <= e < bins[i+1])
        counts.append(count)
    return counts

real_evalue_counts = bin_evalues(real_evalues, evalue_bins)
shuf_evalue_counts = bin_evalues(shuf_evalues, evalue_bins)
shuf_evalue_norm = [c/10 for c in shuf_evalue_counts]

print(f"\n{'E-value Range':<20} {'Real':>12} {'Shuffled':>12} {'Fold':>8}")
print("-"*55)
for label, real_c, shuf_c in zip(bin_labels, real_evalue_counts, shuf_evalue_norm):
    fold = real_c / shuf_c if shuf_c > 0 else float('inf')
    print(f"{label:<20} {real_c:>12,} {shuf_c:>12,.0f} {fold:>8.1f}x")

# =====================================================================
# 3. QUALITY THRESHOLD ANALYSIS
# =====================================================================
print("\n" + "="*70)
print("QUALITY THRESHOLD ANALYSIS")
print("="*70)

thresholds = [
    (0, 0, "All hits"),
    (70, 0, "≥70% identity"),
    (80, 0, "≥80% identity"),
    (90, 0, "≥90% identity"),
    (0, 30, "≥30bp"),
    (0, 50, "≥50bp"),
    (0, 100, "≥100bp"),
    (70, 30, "≥70%, ≥30bp"),
    (80, 50, "≥80%, ≥50bp"),
    (90, 100, "≥90%, ≥100bp"),
]

print(f"\n{'Threshold':<20} {'Real':>12} {'Shuffled':>12} {'Fold':>8}")
print("-"*55)

threshold_results = []
for min_pid, min_len, label in thresholds:
    real_filt = filter_hits(real_hits, min_pid, min_len)
    shuf_filt = filter_hits(shuf_hits_all, min_pid, min_len)
    shuf_norm = len(shuf_filt) / 10
    fold = len(real_filt) / shuf_norm if shuf_norm > 0 else float('inf')
    print(f"{label:<20} {len(real_filt):>12,} {shuf_norm:>12,.0f} {fold:>8.1f}x")
    threshold_results.append({
        'threshold': label,
        'real': len(real_filt),
        'shuffled': shuf_norm,
        'fold': fold
    })

# =====================================================================
# 4. TOP ENRICHED/DEPLETED FAMILIES
# =====================================================================
print("\n" + "="*70)
print("MOST ENRICHED TE FAMILIES (Real vs Shuffled)")
print("="*70)

# Calculate enrichment for all families with enough hits
all_enrichments = []
for family in set(list(real_families.keys()) + list(shuf_families.keys())):
    real_c = real_families.get(family, 0)
    shuf_c = shuf_families_norm.get(family, 0)
    if real_c >= 100 or shuf_c >= 100:  # Require minimum hits
        fold = real_c / shuf_c if shuf_c > 0 else float('inf')
        all_enrichments.append({
            'family': family,
            'real': real_c,
            'shuffled': shuf_c,
            'fold': fold
        })

# Sort by fold enrichment
all_enrichments.sort(key=lambda x: -x['fold'])

print("\nMost ENRICHED in real sequences:")
print(f"{'Family':<20} {'Real':>10} {'Shuffled':>10} {'Fold':>8}")
print("-"*50)
for e in all_enrichments[:10]:
    if e['fold'] != float('inf'):
        print(f"{e['family']:<20} {e['real']:>10,} {e['shuffled']:>10,.0f} {e['fold']:>8.1f}x")

print("\nMost SIMILAR between real and shuffled:")
similar = [e for e in all_enrichments if 0.8 <= e['fold'] <= 1.2 and e['real'] >= 500]
similar.sort(key=lambda x: -x['real'])
for e in similar[:10]:
    print(f"{e['family']:<20} {e['real']:>10,} {e['shuffled']:>10,.0f} {e['fold']:>8.1f}x")

print("\nMost DEPLETED in real (enriched in shuffled):")
all_enrichments.sort(key=lambda x: x['fold'])
for e in all_enrichments[:5]:
    if e['fold'] > 0:
        print(f"{e['family']:<20} {e['real']:>10,} {e['shuffled']:>10,.0f} {e['fold']:>8.1f}x")

# =====================================================================
# 5. CREATE VISUALIZATIONS
# =====================================================================
print("\n" + "="*70)
print("Creating visualizations...")
print("="*70)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Color scheme: contrasting colors for real vs shuffled
COLOR_REAL = '#2980b9'      # Bold blue
COLOR_SHUFFLED = '#e67e22'  # Orange
COLOR_ENRICHED = '#27ae60'  # Green for enrichment bars

# 1. Top TE families comparison
ax = axes[0, 0]
top_n = 15
families = [e['family'] for e in enrichment_data[:top_n]]
real_counts = [e['real'] for e in enrichment_data[:top_n]]
shuf_counts = [e['shuffled'] for e in enrichment_data[:top_n]]

x = np.arange(len(families))
width = 0.35

bars1 = ax.bar(x - width/2, real_counts, width, label='Real', color=COLOR_REAL)
bars2 = ax.bar(x + width/2, shuf_counts, width, label='Shuffled', color=COLOR_SHUFFLED)

ax.set_ylabel('Hit Count')
ax.set_title('Top TE Families: Real vs Shuffled')
ax.set_xticks(x)
ax.set_xticklabels(families, rotation=45, ha='right')
ax.legend()

# 2. Fold enrichment by family - all green since all are enriched >1
ax = axes[0, 1]
folds = [e['fold'] for e in enrichment_data[:top_n]]
ax.barh(families, folds, color=COLOR_ENRICHED)
ax.axvline(1, color='black', linestyle='--', linewidth=2)
ax.set_xlabel('Fold Enrichment (Real/Shuffled)')
ax.set_title('TE Family Enrichment')
ax.invert_yaxis()

# 3. E-value distribution
ax = axes[1, 0]
x = np.arange(len(bin_labels))
ax.bar(x - width/2, real_evalue_counts, width, label='Real', color=COLOR_REAL)
ax.bar(x + width/2, shuf_evalue_norm, width, label='Shuffled', color=COLOR_SHUFFLED)
ax.set_ylabel('Hit Count')
ax.set_xlabel('E-value Range')
ax.set_title('E-value Distribution')
ax.set_xticks(x)
ax.set_xticklabels(bin_labels, rotation=45, ha='right')
ax.legend()

# 4. Fold enrichment by quality threshold - all green since all enriched
ax = axes[1, 1]
labels = [t['threshold'] for t in threshold_results]
folds = [t['fold'] for t in threshold_results]
bars = ax.barh(labels, folds, color=COLOR_ENRICHED)
ax.axvline(1, color='black', linestyle='--', linewidth=2)
ax.set_xlabel('Fold Enrichment (Real/Shuffled)')
ax.set_title('Enrichment by Quality Threshold')
ax.invert_yaxis()

# Add value labels
for bar, fold in zip(bars, folds):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
            f'{fold:.1f}x', va='center', fontsize=9)

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/28_te_family_enrichment.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/28_te_family_enrichment.png")

# Second figure: detailed e-value and length distributions
fig2, axes = plt.subplots(2, 2, figsize=(14, 10))

# Sample for plotting
np.random.seed(42)
real_sample = np.random.choice(len(real_hits), min(50000, len(real_hits)), replace=False)
shuf_sample = np.random.choice(len(shuf_hits_all), min(50000, len(shuf_hits_all)), replace=False)

# 1. Identity distribution
ax = axes[0, 0]
real_pidents = [real_hits[i]['pident'] for i in real_sample]
shuf_pidents = [shuf_hits_all[i]['pident'] for i in shuf_sample]
bins = np.linspace(60, 100, 40)
ax.hist(shuf_pidents, bins=bins, alpha=0.7, label='Shuffled', color=COLOR_SHUFFLED, density=True)
ax.hist(real_pidents, bins=bins, alpha=0.7, label='Real', color=COLOR_REAL, density=True)
ax.set_xlabel('Percent Identity')
ax.set_ylabel('Density')
ax.set_title('Identity Distribution')
ax.legend()

# 2. Length distribution
ax = axes[0, 1]
real_lengths = [real_hits[i]['length'] for i in real_sample]
shuf_lengths = [shuf_hits_all[i]['length'] for i in shuf_sample]
bins = np.linspace(0, 200, 50)
ax.hist(shuf_lengths, bins=bins, alpha=0.7, label='Shuffled', color=COLOR_SHUFFLED, density=True)
ax.hist(real_lengths, bins=bins, alpha=0.7, label='Real', color=COLOR_REAL, density=True)
ax.set_xlabel('Alignment Length (bp)')
ax.set_ylabel('Density')
ax.set_title('Length Distribution')
ax.legend()

# 3. E-value distribution (log scale)
ax = axes[1, 0]
real_evalues_sample = [max(real_hits[i]['evalue'], 1e-180) for i in real_sample]
shuf_evalues_sample = [max(shuf_hits_all[i]['evalue'], 1e-180) for i in shuf_sample]
bins = np.logspace(-180, 1, 50)
ax.hist(shuf_evalues_sample, bins=bins, alpha=0.7, label='Shuffled', color=COLOR_SHUFFLED, density=True)
ax.hist(real_evalues_sample, bins=bins, alpha=0.7, label='Real', color=COLOR_REAL, density=True)
ax.set_xscale('log')
ax.set_xlabel('E-value')
ax.set_ylabel('Density')
ax.set_title('E-value Distribution')
ax.legend()

# 4. 2D: Identity vs Length colored by source
ax = axes[1, 1]
real_sample_small = np.random.choice(len(real_hits), min(5000, len(real_hits)), replace=False)
shuf_sample_small = np.random.choice(len(shuf_hits_all), min(5000, len(shuf_hits_all)), replace=False)

ax.scatter([shuf_hits_all[i]['pident'] for i in shuf_sample_small],
           [shuf_hits_all[i]['length'] for i in shuf_sample_small],
           alpha=0.4, s=8, c=COLOR_SHUFFLED, label='Shuffled', edgecolors='none')
ax.scatter([real_hits[i]['pident'] for i in real_sample_small],
           [real_hits[i]['length'] for i in real_sample_small],
           alpha=0.4, s=8, c=COLOR_REAL, label='Real', edgecolors='none')
ax.set_xlabel('Percent Identity')
ax.set_ylabel('Alignment Length (bp)')
ax.set_title('Identity vs Length')
ax.legend()

plt.tight_layout()
plt.savefig('figures/repeatmasker_comparison/29_real_vs_shuffled_distributions.png', dpi=150, bbox_inches='tight')
print("Saved: figures/repeatmasker_comparison/29_real_vs_shuffled_distributions.png")

# Save detailed results
results_file = Path('results/shuffled_controls/te_family_enrichment.tsv')
with open(results_file, 'w') as f:
    f.write("family\treal_hits\tshuffled_hits\tfold_enrichment\n")
    for e in sorted(all_enrichments, key=lambda x: -x['real']):
        f.write(f"{e['family']}\t{e['real']}\t{e['shuffled']:.0f}\t{e['fold']:.2f}\n")
print(f"Saved: {results_file}")

plt.close('all')
