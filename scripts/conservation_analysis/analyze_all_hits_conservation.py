#!/usr/bin/env python3
"""Analyze conservation scores for ALL TE hits vs high-confidence subset."""

import sys
from collections import defaultdict

def parse_conservation_file(filepath):
    """Parse bigWigAverageOverBed output with hit metadata."""
    hits = []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            mean_score = float(parts[5])  # mean (not mean0)

            # Parse name: prefix_num|transcript|te|pident|length
            name_parts = name.split('|')
            category = name_parts[0].split('_')[0]  # 'novel' or 'known'
            transcript = name_parts[1]
            te_id = name_parts[2]
            pident = float(name_parts[3])
            length = int(name_parts[4])

            hits.append({
                'category': category,
                'transcript': transcript,
                'te_id': te_id,
                'pident': pident,
                'length': length,
                'conservation': mean_score
            })
    return hits

print("Loading ALL hits conservation data...", file=sys.stderr)
all_hits = parse_conservation_file('results/repeatmasker_analysis/te_hits_all_conservation.tab')
print(f"Loaded {len(all_hits):,} hits", file=sys.stderr)

# Overall statistics
cons_scores = [h['conservation'] for h in all_hits]
print(f"\n{'='*70}")
print("CONSERVATION ANALYSIS: ALL TE HITS (n={:,})".format(len(all_hits)))
print(f"{'='*70}\n")

# Basic stats
cons_scores_sorted = sorted(cons_scores)
n = len(cons_scores)
mean_cons = sum(cons_scores) / n
median_cons = cons_scores_sorted[n // 2]

print(f"Mean conservation:   {mean_cons:.3f}")
print(f"Median conservation: {median_cons:.3f}")
print(f"Min: {min(cons_scores):.3f}, Max: {max(cons_scores):.3f}")

# Distribution by conservation category
print(f"\n{'='*70}")
print("CONSERVATION CATEGORIES")
print(f"{'='*70}\n")

cats = {
    'Highly conserved (≥2)': lambda x: x >= 2,
    'Conserved (1-2)': lambda x: 1 <= x < 2,
    'Weakly conserved (0-1)': lambda x: 0 <= x < 1,
    'Fast-evolving (<0)': lambda x: x < 0
}

for cat_name, test in cats.items():
    count = sum(1 for s in cons_scores if test(s))
    pct = 100 * count / n
    print(f"{cat_name:<30} {count:>10,} ({pct:>5.1f}%)")

# Break down by quality thresholds
print(f"\n{'='*70}")
print("CONSERVATION BY QUALITY THRESHOLD")
print(f"{'='*70}\n")

thresholds = [
    ('All hits', lambda h: True),
    ('≥60% identity', lambda h: h['pident'] >= 60),
    ('≥70% identity', lambda h: h['pident'] >= 70),
    ('≥80% identity', lambda h: h['pident'] >= 80),
    ('≥90% identity', lambda h: h['pident'] >= 90),
    ('≥80% id, ≥30bp', lambda h: h['pident'] >= 80 and h['length'] >= 30),
    ('≥80% id, ≥50bp', lambda h: h['pident'] >= 80 and h['length'] >= 50),
    ('≥80% id, ≥100bp', lambda h: h['pident'] >= 80 and h['length'] >= 100),
]

print(f"{'Threshold':<25} {'N hits':>12} {'Mean cons':>12} {'Median cons':>12} {'% cons>1':>10}")
print("-" * 75)

for label, test in thresholds:
    subset = [h for h in all_hits if test(h)]
    if not subset:
        continue
    scores = [h['conservation'] for h in subset]
    scores_sorted = sorted(scores)
    m = sum(scores) / len(scores)
    med = scores_sorted[len(scores) // 2]
    pct_cons = 100 * sum(1 for s in scores if s >= 1) / len(scores)
    print(f"{label:<25} {len(subset):>12,} {m:>12.3f} {med:>12.3f} {pct_cons:>9.1f}%")

# Novel vs Known
print(f"\n{'='*70}")
print("NOVEL vs KNOWN (RepeatMasker)")
print(f"{'='*70}\n")

novel = [h for h in all_hits if h['category'] == 'novel']
known = [h for h in all_hits if h['category'] == 'known']

for label, subset in [('Novel (not in RM)', novel), ('Known (in RM)', known)]:
    scores = [h['conservation'] for h in subset]
    scores_sorted = sorted(scores)
    m = sum(scores) / len(scores)
    med = scores_sorted[len(scores) // 2]
    pct_cons = 100 * sum(1 for s in scores if s >= 1) / len(scores)
    print(f"{label:<25} n={len(subset):>10,}  mean={m:.3f}  median={med:.3f}  %cons>1={pct_cons:.1f}%")

# Top TE families by conservation
print(f"\n{'='*70}")
print("TOP 20 TE FAMILIES BY MEAN CONSERVATION (n≥100)")
print(f"{'='*70}\n")

from collections import defaultdict
family_scores = defaultdict(list)
for h in all_hits:
    # Extract family from TE ID (e.g., FBti0060993 -> look up)
    family_scores[h['te_id']].append(h['conservation'])

# Group by TE family (simple: use first 4 chars of te_id as proxy, or just use te_id)
# Actually let's aggregate all scores
family_means = []
for te_id, scores in family_scores.items():
    if len(scores) >= 100:
        mean_s = sum(scores) / len(scores)
        family_means.append((te_id, len(scores), mean_s))

family_means.sort(key=lambda x: -x[2])
print(f"{'TE ID':<20} {'N hits':>10} {'Mean cons':>12}")
print("-" * 45)
for te_id, n_hits, mean_s in family_means[:20]:
    print(f"{te_id:<20} {n_hits:>10,} {mean_s:>12.3f}")

# Conservation score distribution for plotting
print(f"\n{'='*70}")
print("DISTRIBUTION FOR PLOTTING")
print(f"{'='*70}\n")

bins = [(-5, -2), (-2, -1), (-1, 0), (0, 0.5), (0.5, 1), (1, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 10)]
print(f"{'Bin':<15} {'Count':>12} {'Percent':>10}")
print("-" * 40)
for low, high in bins:
    count = sum(1 for s in cons_scores if low <= s < high)
    pct = 100 * count / n
    print(f"[{low:>4}, {high:>4}){'':<5} {count:>12,} {pct:>9.1f}%")
