#!/usr/bin/env python3
"""
Analyze correlation between protein disorder and exon TE content.

Also compares exons that overlap UTRs vs those that don't.
"""

import sys
from collections import defaultdict
from pathlib import Path

def load_disorder_scores(path):
    """Load disorder predictions."""
    scores = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                fbgn = parts[0]
                disorder_fraction = float(parts[3])
                scores[fbgn] = disorder_fraction
    return scores

def load_gene_te_stats(path):
    """Load gene-level TE statistics."""
    stats = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                fbgn = parts[1]
                stats[fbgn] = {
                    'gene_symbol': parts[2],
                    'exons': int(parts[3]),
                    'hits': int(parts[5]),
                    'density': float(parts[8]),
                }
    return stats

def load_exon_metadata(path):
    """Load exon-level metadata with UTR overlap info."""
    exons = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            exons[row['exon_id']] = {
                'fbgn': row['fbgn'],
                'position': row['position'],
                'utr_overlap': row['utr_overlap'],
                'length': int(row['length'])
            }
    return exons

def load_exon_te_stats(path):
    """Load per-exon TE hit counts."""
    stats = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                exon_id = parts[0]
                stats[exon_id] = {
                    'hits': int(parts[7]),
                    'hit_bp': int(parts[8]),
                    'density': float(parts[9]) if parts[9] else 0,
                }
    return stats

def correlation(x, y):
    """Calculate Pearson correlation coefficient."""
    n = len(x)
    if n < 3:
        return 0, 1.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    den_x = sum((xi - mean_x) ** 2 for xi in x) ** 0.5
    den_y = sum((yi - mean_y) ** 2 for yi in y) ** 0.5

    if den_x == 0 or den_y == 0:
        return 0, 1.0

    r = num / (den_x * den_y)

    # Approximate p-value using t-distribution
    import math
    if abs(r) >= 1:
        return r, 0.0
    t = r * math.sqrt((n - 2) / (1 - r**2))
    # Very rough p-value approximation
    p = 2 * (1 - min(0.9999, abs(t) / (abs(t) + n)))

    return r, p

def main():
    # Paths
    disorder_path = Path('data/annotations/gene_disorder_scores.tsv')
    gene_te_path = Path('results/exon_analysis/gene_exon_te_summary.tsv')
    exon_metadata_path = Path('data/queries/genome_wide/exon_metadata.tsv')
    exon_te_path = Path('results/exon_analysis/exon_te_summary.tsv')

    # Load data
    print("Loading data...")
    disorder = load_disorder_scores(disorder_path)
    gene_te = load_gene_te_stats(gene_te_path)
    exon_meta = load_exon_metadata(exon_metadata_path)
    exon_te = load_exon_te_stats(exon_te_path)

    print(f"  Disorder scores: {len(disorder):,} genes")
    print(f"  Gene TE stats: {len(gene_te):,} genes")
    print(f"  Exon metadata: {len(exon_meta):,} exons")
    print(f"  Exon TE stats: {len(exon_te):,} exons with hits")

    # ========================================
    # PART 1: Disorder vs TE Density Correlation
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: DISORDER vs TE DENSITY CORRELATION")
    print("=" * 70)

    # Merge data
    merged = []
    for fbgn in gene_te:
        if fbgn in disorder:
            merged.append({
                'fbgn': fbgn,
                'symbol': gene_te[fbgn]['gene_symbol'],
                'disorder': disorder[fbgn],
                'te_density': gene_te[fbgn]['density'],
                'te_hits': gene_te[fbgn]['hits'],
            })

    print(f"\nGenes with both disorder and TE data: {len(merged):,}")

    # Calculate correlation
    disorder_vals = [m['disorder'] for m in merged]
    density_vals = [m['te_density'] for m in merged]
    hits_vals = [m['te_hits'] for m in merged]

    r_density, p_density = correlation(disorder_vals, density_vals)
    r_hits, p_hits = correlation(disorder_vals, hits_vals)

    print(f"\nCorrelations:")
    print(f"  Disorder vs TE Density: r = {r_density:.4f}")
    print(f"  Disorder vs TE Hits:    r = {r_hits:.4f}")

    # Bin by disorder level
    low_disorder = [m for m in merged if m['disorder'] < 0.3]
    med_disorder = [m for m in merged if 0.3 <= m['disorder'] < 0.6]
    high_disorder = [m for m in merged if m['disorder'] >= 0.6]

    def avg(lst, key):
        if not lst:
            return 0
        return sum(m[key] for m in lst) / len(lst)

    print(f"\nTE Density by Disorder Level:")
    print(f"  {'Disorder Level':<20} {'N Genes':>10} {'Mean Density':>15} {'Mean Hits':>12}")
    print(f"  {'-'*20} {'-'*10} {'-'*15} {'-'*12}")
    print(f"  {'Low (<30%)':<20} {len(low_disorder):>10,} {avg(low_disorder, 'te_density'):>15.1f} {avg(low_disorder, 'te_hits'):>12.1f}")
    print(f"  {'Medium (30-60%)':<20} {len(med_disorder):>10,} {avg(med_disorder, 'te_density'):>15.1f} {avg(med_disorder, 'te_hits'):>12.1f}")
    print(f"  {'High (>60%)':<20} {len(high_disorder):>10,} {avg(high_disorder, 'te_density'):>15.1f} {avg(high_disorder, 'te_hits'):>12.1f}")

    # Top disordered genes with high TE content
    high_both = sorted([m for m in merged if m['disorder'] > 0.5 and m['te_density'] > 100],
                       key=lambda x: -x['te_density'])[:15]

    print(f"\nTop 15 High-Disorder Genes with High TE Density:")
    print(f"  {'Gene':<15} {'Disorder%':>10} {'TE Density':>12} {'Hits':>8}")
    print(f"  {'-'*15} {'-'*10} {'-'*12} {'-'*8}")
    for m in high_both:
        print(f"  {m['symbol']:<15} {m['disorder']*100:>9.1f}% {m['te_density']:>12.1f} {m['te_hits']:>8,}")

    # ========================================
    # PART 2: UTR Overlap Analysis
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: UTR OVERLAP vs NO OVERLAP COMPARISON")
    print("=" * 70)

    # Categorize exons by UTR overlap
    overlap_stats = defaultdict(lambda: {'count': 0, 'with_hits': 0, 'total_hits': 0, 'total_bp': 0, 'total_len': 0})

    for exon_id, meta in exon_meta.items():
        utr_overlap = meta['utr_overlap']
        overlap_stats[utr_overlap]['count'] += 1
        overlap_stats[utr_overlap]['total_len'] += meta['length']

        if exon_id in exon_te:
            te = exon_te[exon_id]
            overlap_stats[utr_overlap]['with_hits'] += 1
            overlap_stats[utr_overlap]['total_hits'] += te['hits']
            overlap_stats[utr_overlap]['total_bp'] += te['hit_bp']

    print(f"\nExon Statistics by UTR Overlap:")
    print(f"  {'UTR Overlap':<12} {'Exons':>10} {'With Hits':>12} {'Total Hits':>12} {'Avg Hits':>10} {'Density':>10}")
    print(f"  {'-'*12} {'-'*10} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")

    for overlap in ['none', '5utr', '3utr', 'both']:
        s = overlap_stats[overlap]
        avg_hits = s['total_hits'] / s['count'] if s['count'] > 0 else 0
        density = (s['total_bp'] / s['total_len'] * 100) if s['total_len'] > 0 else 0
        pct_with = 100 * s['with_hits'] / s['count'] if s['count'] > 0 else 0
        print(f"  {overlap:<12} {s['count']:>10,} {s['with_hits']:>10,} ({pct_with:>4.1f}%) {s['total_hits']:>10,} {avg_hits:>10.1f} {density:>10.1f}")

    # Compare internal (no overlap) vs terminal (any overlap)
    print(f"\n  Summary:")
    internal = overlap_stats['none']
    terminal_hits = sum(overlap_stats[k]['total_hits'] for k in ['5utr', '3utr', 'both'])
    terminal_count = sum(overlap_stats[k]['count'] for k in ['5utr', '3utr', 'both'])
    terminal_len = sum(overlap_stats[k]['total_len'] for k in ['5utr', '3utr', 'both'])
    terminal_bp = sum(overlap_stats[k]['total_bp'] for k in ['5utr', '3utr', 'both'])

    internal_density = (internal['total_bp'] / internal['total_len'] * 100) if internal['total_len'] > 0 else 0
    terminal_density = (terminal_bp / terminal_len * 100) if terminal_len > 0 else 0

    print(f"  - Internal exons (no UTR overlap): {internal['count']:,} exons, {internal['total_hits']:,} hits, density {internal_density:.1f}")
    print(f"  - Terminal exons (any UTR overlap): {terminal_count:,} exons, {terminal_hits:,} hits, density {terminal_density:.1f}")
    print(f"  - Ratio (terminal/internal density): {terminal_density/internal_density:.2f}x")

    # By position type
    print(f"\n  By Exon Position:")
    position_stats = defaultdict(lambda: {'count': 0, 'hits': 0, 'bp': 0, 'len': 0})

    for exon_id, meta in exon_meta.items():
        pos = meta['position']
        position_stats[pos]['count'] += 1
        position_stats[pos]['len'] += meta['length']
        if exon_id in exon_te:
            position_stats[pos]['hits'] += exon_te[exon_id]['hits']
            position_stats[pos]['bp'] += exon_te[exon_id]['hit_bp']

    print(f"  {'Position':<15} {'Exons':>10} {'Hits':>12} {'Avg Hits':>10} {'Density':>10}")
    print(f"  {'-'*15} {'-'*10} {'-'*12} {'-'*10} {'-'*10}")
    for pos in ['first_exon', 'internal_exon', 'last_exon', 'single_exon']:
        s = position_stats[pos]
        avg_hits = s['hits'] / s['count'] if s['count'] > 0 else 0
        density = (s['bp'] / s['len'] * 100) if s['len'] > 0 else 0
        print(f"  {pos:<15} {s['count']:>10,} {s['hits']:>12,} {avg_hits:>10.1f} {density:>10.1f}")

    # ========================================
    # PART 3: Combined Analysis
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: INTERNAL EXONS ONLY - DISORDER CORRELATION")
    print("=" * 70)

    # Get genes that have internal exons
    internal_exon_genes = defaultdict(lambda: {'hits': 0, 'bp': 0, 'len': 0})
    for exon_id, meta in exon_meta.items():
        if meta['utr_overlap'] == 'none':  # No UTR overlap
            fbgn = meta['fbgn']
            internal_exon_genes[fbgn]['len'] += meta['length']
            if exon_id in exon_te:
                internal_exon_genes[fbgn]['hits'] += exon_te[exon_id]['hits']
                internal_exon_genes[fbgn]['bp'] += exon_te[exon_id]['hit_bp']

    # Merge with disorder
    internal_merged = []
    for fbgn, stats in internal_exon_genes.items():
        if fbgn in disorder and stats['len'] > 0:
            density = (stats['bp'] / stats['len']) * 100
            internal_merged.append({
                'fbgn': fbgn,
                'disorder': disorder[fbgn],
                'te_density': density,
                'te_hits': stats['hits'],
            })

    print(f"\nGenes with internal exons and disorder data: {len(internal_merged):,}")

    # Correlation for internal exons only
    if internal_merged:
        disorder_vals = [m['disorder'] for m in internal_merged]
        density_vals = [m['te_density'] for m in internal_merged]
        r, _ = correlation(disorder_vals, density_vals)
        print(f"Disorder vs TE Density (internal exons only): r = {r:.4f}")

        # Bin by disorder
        low = [m for m in internal_merged if m['disorder'] < 0.3]
        med = [m for m in internal_merged if 0.3 <= m['disorder'] < 0.6]
        high = [m for m in internal_merged if m['disorder'] >= 0.6]

        print(f"\nInternal Exon TE Density by Disorder Level:")
        print(f"  {'Disorder Level':<20} {'N Genes':>10} {'Mean Density':>15}")
        print(f"  {'-'*20} {'-'*10} {'-'*15}")
        print(f"  {'Low (<30%)':<20} {len(low):>10,} {avg(low, 'te_density'):>15.1f}")
        print(f"  {'Medium (30-60%)':<20} {len(med):>10,} {avg(med, 'te_density'):>15.1f}")
        print(f"  {'High (>60%)':<20} {len(high):>10,} {avg(high, 'te_density'):>15.1f}")


if __name__ == '__main__':
    main()
