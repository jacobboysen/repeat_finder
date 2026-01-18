#!/usr/bin/env python3
"""
Analyze overlap between exon BLAST TE hits and RepeatMasker annotations.

This script:
1. Loads exon metadata with genomic coordinates
2. Parses RepeatMasker .out file
3. For each exon BLAST hit, determines if it overlaps a known RM region
4. Calculates distance to nearest RM hit for all hits
5. Generates distance distribution plots
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from bisect import bisect_left

# Add scripts directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_references_dir, get_results_dir
from utils.blast_io import load_blast_results


def parse_repeatmasker_out(rm_file):
    """
    Parse RepeatMasker .out file.

    Returns list of dicts with repeat annotations.
    """
    repeats = []

    with open(rm_file) as f:
        # Skip header lines (first 3 lines)
        for _ in range(3):
            next(f)

        for line in f:
            line = line.strip()
            if not line:
                continue

            # RepeatMasker output is space-delimited with variable spacing
            parts = line.split()
            if len(parts) < 14:
                continue

            try:
                score = int(parts[0])
                pct_div = float(parts[1])
                pct_del = float(parts[2])
                pct_ins = float(parts[3])
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                strand = '+' if parts[8] == '+' else '-'
                repeat_name = parts[9]
                repeat_class = parts[10]

                # Normalize chromosome names (chr2L -> 2L)
                if chrom.startswith('chr'):
                    chrom = chrom[3:]

                repeats.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'repeat_name': repeat_name,
                    'repeat_class': repeat_class,
                    'score': score,
                    'pct_divergence': pct_div,
                })
            except (ValueError, IndexError):
                continue

    return repeats


def build_rm_index(repeats):
    """
    Build chromosome-indexed interval structure for fast overlap queries.

    Returns dict of {chrom: sorted list of (start, end, repeat_info)}
    """
    by_chrom = defaultdict(list)

    for r in repeats:
        by_chrom[r['chrom']].append({
            'start': r['start'],
            'end': r['end'],
            'repeat_name': r['repeat_name'],
            'repeat_class': r['repeat_class'],
            'divergence': r['pct_divergence'],
        })

    # Sort by start position
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x['start'])

    return dict(by_chrom)


def build_rm_starts_index(rm_index):
    """
    Build sorted list of start positions for binary search.
    """
    starts_by_chrom = {}
    for chrom, intervals in rm_index.items():
        starts_by_chrom[chrom] = [iv['start'] for iv in intervals]
    return starts_by_chrom


def load_exon_metadata(path):
    """
    Load exon metadata with genomic coordinates.
    """
    exons = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            exons[row['exon_id']] = {
                'chrom': row['chrom'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'length': int(row['length']),
                'fbgn': row['fbgn'],
                'gene_symbol': row['gene_symbol'],
                'position': row['position'],
                'utr_overlap': row['utr_overlap'],
            }
    return exons


def find_overlap(query_start, query_end, rm_intervals, rm_starts):
    """
    Find if query interval overlaps any RM intervals.

    Uses binary search for efficiency.
    Returns (is_overlap, matching_rm_info or None)
    """
    if not rm_intervals:
        return False, None

    # Binary search to find starting point
    idx = bisect_left(rm_starts, query_start)

    # Check intervals around this position
    # Look back in case we landed past an overlapping interval
    start_idx = max(0, idx - 10)
    end_idx = min(len(rm_intervals), idx + 10)

    for i in range(start_idx, end_idx):
        rm = rm_intervals[i]
        # Check overlap
        if query_start <= rm['end'] and query_end >= rm['start']:
            return True, rm
        # Can stop if we've gone past all possible overlaps
        if rm['start'] > query_end:
            break

    return False, None


def find_nearest_rm(query_start, query_end, rm_intervals, rm_starts):
    """
    Find distance to nearest RM interval.

    Returns (distance, nearest_rm_info)
    - distance = 0 means overlap
    - negative distance is not used (always positive or 0)
    """
    if not rm_intervals:
        return float('inf'), None

    # Binary search to find insertion point
    idx = bisect_left(rm_starts, query_start)

    min_dist = float('inf')
    nearest_rm = None

    # Check intervals around this position
    for i in range(max(0, idx - 1), min(len(rm_intervals), idx + 2)):
        rm = rm_intervals[i]

        # Calculate distance (0 if overlapping)
        if query_start <= rm['end'] and query_end >= rm['start']:
            return 0, rm
        elif query_end < rm['start']:
            # Query is entirely before RM
            dist = rm['start'] - query_end
        else:
            # Query is entirely after RM
            dist = query_start - rm['end']

        if dist < min_dist:
            min_dist = dist
            nearest_rm = rm

    return min_dist, nearest_rm


def analyze_hits(blast_df, exon_meta, rm_index, rm_starts_index):
    """
    Analyze each BLAST hit for RM overlap and distance.
    """
    results = []

    for _, hit in blast_df.iterrows():
        qseqid = hit['qseqid']
        qstart = hit['qstart']
        qend = hit['qend']

        if qseqid not in exon_meta:
            continue

        exon = exon_meta[qseqid]
        chrom = exon['chrom']

        # Convert query-relative coordinates to genomic coordinates
        if exon['strand'] == '+':
            genomic_start = exon['start'] + qstart - 1
            genomic_end = exon['start'] + qend - 1
        else:
            # For minus strand, coordinates are flipped
            genomic_end = exon['end'] - qstart + 1
            genomic_start = exon['end'] - qend + 1

        # Find overlap and nearest RM
        rm_intervals = rm_index.get(chrom, [])
        rm_starts = rm_starts_index.get(chrom, [])

        is_overlap, overlap_rm = find_overlap(genomic_start, genomic_end, rm_intervals, rm_starts)
        distance, nearest_rm = find_nearest_rm(genomic_start, genomic_end, rm_intervals, rm_starts)

        results.append({
            'qseqid': qseqid,
            'qstart': qstart,
            'qend': qend,
            'sseqid': hit['sseqid'],
            'pident': hit['pident'],
            'length': hit['length'],
            'evalue': hit['evalue'],
            'bitscore': hit['bitscore'],
            'chrom': chrom,
            'genomic_start': genomic_start,
            'genomic_end': genomic_end,
            'exon_strand': exon['strand'],
            'gene_symbol': exon['gene_symbol'],
            'position': exon['position'],
            'utr_overlap': exon['utr_overlap'],
            'rm_overlap': is_overlap,
            'rm_distance': distance,
            'rm_name': nearest_rm['repeat_name'] if nearest_rm else None,
            'rm_class': nearest_rm['repeat_class'] if nearest_rm else None,
            'rm_divergence': nearest_rm['divergence'] if nearest_rm else None,
        })

    return results


def generate_distance_histogram(results, output_path, max_distance=10000):
    """
    Generate distance distribution plot using matplotlib.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Get distances (excluding overlaps for distance plot)
    distances = [r['rm_distance'] for r in results if r['rm_distance'] > 0 and r['rm_distance'] < max_distance]
    overlap_count = sum(1 for r in results if r['rm_distance'] == 0)
    beyond_count = sum(1 for r in results if r['rm_distance'] >= max_distance)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Overall histogram with log scale
    ax1 = axes[0, 0]
    bins = np.logspace(0, np.log10(max_distance), 50)
    ax1.hist(distances, bins=bins, edgecolor='black', alpha=0.7)
    ax1.set_xscale('log')
    ax1.set_xlabel('Distance to Nearest RepeatMasker Hit (bp)')
    ax1.set_ylabel('Number of BLAST Hits')
    ax1.set_title(f'Distance Distribution (log scale)\nOverlapping: {overlap_count:,} | Beyond {max_distance/1000:.0f}kb: {beyond_count:,}')
    ax1.axvline(x=100, color='red', linestyle='--', alpha=0.5, label='100bp')
    ax1.axvline(x=1000, color='orange', linestyle='--', alpha=0.5, label='1kb')
    ax1.legend()

    # Plot 2: Linear scale for close hits (< 1kb)
    ax2 = axes[0, 1]
    close_distances = [d for d in distances if d < 1000]
    ax2.hist(close_distances, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.set_xlabel('Distance to Nearest RepeatMasker Hit (bp)')
    ax2.set_ylabel('Number of BLAST Hits')
    ax2.set_title(f'Distance Distribution (< 1kb)\nn = {len(close_distances):,} hits')

    # Plot 3: Cumulative distribution
    ax3 = axes[1, 0]
    sorted_dist = sorted(distances)
    cumulative = np.arange(1, len(sorted_dist) + 1) / len(sorted_dist)
    ax3.plot(sorted_dist, cumulative, linewidth=2)
    ax3.set_xscale('log')
    ax3.set_xlabel('Distance to Nearest RepeatMasker Hit (bp)')
    ax3.set_ylabel('Cumulative Fraction')
    ax3.set_title('Cumulative Distribution of Distances')
    ax3.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
    ax3.axhline(y=0.9, color='orange', linestyle='--', alpha=0.5)
    ax3.grid(True, alpha=0.3)

    # Find median and 90th percentile
    if sorted_dist:
        median_dist = sorted_dist[len(sorted_dist) // 2]
        p90_dist = sorted_dist[int(len(sorted_dist) * 0.9)]
        ax3.axvline(x=median_dist, color='red', linestyle='--', alpha=0.5, label=f'Median: {median_dist:,.0f}bp')
        ax3.axvline(x=p90_dist, color='orange', linestyle='--', alpha=0.5, label=f'90th pct: {p90_dist:,.0f}bp')
        ax3.legend()

    # Plot 4: Pie chart of overlap vs non-overlap
    ax4 = axes[1, 1]
    labels = ['Overlapping RM', 'Within 100bp', 'Within 1kb', 'Within 10kb', 'Beyond 10kb']
    within_100 = sum(1 for d in distances if d <= 100)
    within_1k = sum(1 for d in distances if 100 < d <= 1000)
    within_10k = sum(1 for d in distances if 1000 < d <= 10000)
    beyond_10k = sum(1 for r in results if r['rm_distance'] > 10000)

    sizes = [overlap_count, within_100, within_1k, within_10k, beyond_10k]
    colors = ['#ff6b6b', '#ffa94d', '#69db7c', '#74c0fc', '#e9ecef']
    explode = (0.05, 0, 0, 0, 0)

    ax4.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax4.set_title('BLAST Hits by Distance to RepeatMasker')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved distance plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--rm-file',
        type=Path,
        default=get_references_dir() / 'dm6.fa.out',
        help='RepeatMasker .out file'
    )
    parser.add_argument(
        '--exon-metadata',
        type=Path,
        default=Path('data/queries/genome_wide/exon_metadata.tsv'),
        help='Exon metadata with genomic coordinates'
    )
    parser.add_argument(
        '--blast-file',
        type=Path,
        default=get_results_dir() / 'exon_analysis' / 'genome_wide_exons.tsv',
        help='BLAST results TSV'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=get_results_dir() / 'exon_analysis' / 'repeatmasker_overlap',
        help='Output directory'
    )

    args = parser.parse_args()

    print("=" * 70)
    print("Exon BLAST Hits vs RepeatMasker Overlap Analysis")
    print("=" * 70)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load RepeatMasker annotations
    print(f"\nLoading RepeatMasker file: {args.rm_file}")
    repeats = parse_repeatmasker_out(args.rm_file)
    print(f"  Loaded {len(repeats):,} repeat annotations")

    # Build RM index
    print("\nBuilding RepeatMasker interval index...")
    rm_index = build_rm_index(repeats)
    rm_starts_index = build_rm_starts_index(rm_index)
    print(f"  Indexed {len(rm_index)} chromosomes")

    # Load exon metadata
    print(f"\nLoading exon metadata: {args.exon_metadata}")
    exon_meta = load_exon_metadata(args.exon_metadata)
    print(f"  Loaded {len(exon_meta):,} exons")

    # Load BLAST results
    print(f"\nLoading BLAST results: {args.blast_file}")
    blast_df = load_blast_results(args.blast_file)
    print(f"  Loaded {len(blast_df):,} BLAST hits")

    # Analyze hits
    print("\nAnalyzing BLAST hits for RepeatMasker overlap...")
    results = analyze_hits(blast_df, exon_meta, rm_index, rm_starts_index)
    print(f"  Analyzed {len(results):,} hits with valid exon coordinates")

    # Summary statistics
    overlap_count = sum(1 for r in results if r['rm_overlap'])
    non_overlap_count = len(results) - overlap_count

    print(f"\n{'='*70}")
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"\nTotal hits analyzed: {len(results):,}")
    print(f"  Overlapping RepeatMasker: {overlap_count:,} ({100*overlap_count/len(results):.1f}%)")
    print(f"  NOT overlapping RepeatMasker: {non_overlap_count:,} ({100*non_overlap_count/len(results):.1f}%)")

    # Distance statistics for non-overlapping hits
    non_overlap_distances = [r['rm_distance'] for r in results if not r['rm_overlap'] and r['rm_distance'] < float('inf')]
    if non_overlap_distances:
        sorted_dist = sorted(non_overlap_distances)
        median_dist = sorted_dist[len(sorted_dist) // 2]
        mean_dist = sum(sorted_dist) / len(sorted_dist)
        p10 = sorted_dist[int(len(sorted_dist) * 0.1)]
        p90 = sorted_dist[int(len(sorted_dist) * 0.9)]

        print(f"\nDistance to nearest RM (non-overlapping hits):")
        print(f"  Mean distance: {mean_dist:,.0f} bp")
        print(f"  Median distance: {median_dist:,.0f} bp")
        print(f"  10th percentile: {p10:,.0f} bp")
        print(f"  90th percentile: {p90:,.0f} bp")

        # Bin by distance
        within_100 = sum(1 for d in sorted_dist if d <= 100)
        within_500 = sum(1 for d in sorted_dist if d <= 500)
        within_1k = sum(1 for d in sorted_dist if d <= 1000)
        within_5k = sum(1 for d in sorted_dist if d <= 5000)
        within_10k = sum(1 for d in sorted_dist if d <= 10000)

        print(f"\n  Within 100bp:   {within_100:>10,} ({100*within_100/len(sorted_dist):>5.1f}%)")
        print(f"  Within 500bp:   {within_500:>10,} ({100*within_500/len(sorted_dist):>5.1f}%)")
        print(f"  Within 1kb:     {within_1k:>10,} ({100*within_1k/len(sorted_dist):>5.1f}%)")
        print(f"  Within 5kb:     {within_5k:>10,} ({100*within_5k/len(sorted_dist):>5.1f}%)")
        print(f"  Within 10kb:    {within_10k:>10,} ({100*within_10k/len(sorted_dist):>5.1f}%)")

    # Breakdown by RM class for overlapping hits
    rm_class_counts = defaultdict(int)
    for r in results:
        if r['rm_overlap'] and r['rm_class']:
            rm_class_counts[r['rm_class']] += 1

    print(f"\nOverlapping hits by RepeatMasker class:")
    for cls, count in sorted(rm_class_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {cls:<25} {count:>10,} ({100*count/overlap_count:>5.1f}%)")

    # Compare overlapping vs non-overlapping by exon position
    print(f"\nOverlap rate by exon position:")
    position_stats = defaultdict(lambda: {'total': 0, 'overlap': 0})
    for r in results:
        pos = r['position']
        position_stats[pos]['total'] += 1
        if r['rm_overlap']:
            position_stats[pos]['overlap'] += 1

    print(f"  {'Position':<15} {'Total':>12} {'Overlapping':>12} {'Rate':>8}")
    print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*8}")
    for pos in ['first_exon', 'internal_exon', 'last_exon', 'single_exon']:
        s = position_stats[pos]
        rate = 100 * s['overlap'] / s['total'] if s['total'] > 0 else 0
        print(f"  {pos:<15} {s['total']:>12,} {s['overlap']:>12,} {rate:>7.1f}%")

    # Compare by UTR overlap
    print(f"\nOverlap rate by UTR status:")
    utr_stats = defaultdict(lambda: {'total': 0, 'overlap': 0})
    for r in results:
        utr = r['utr_overlap']
        utr_stats[utr]['total'] += 1
        if r['rm_overlap']:
            utr_stats[utr]['overlap'] += 1

    print(f"  {'UTR Overlap':<12} {'Total':>12} {'Overlapping':>12} {'Rate':>8}")
    print(f"  {'-'*12} {'-'*12} {'-'*12} {'-'*8}")
    for utr in ['none', '5utr', '3utr', 'both']:
        s = utr_stats[utr]
        rate = 100 * s['overlap'] / s['total'] if s['total'] > 0 else 0
        print(f"  {utr:<12} {s['total']:>12,} {s['overlap']:>12,} {rate:>7.1f}%")

    # Save detailed results
    overlap_file = args.output_dir / 'exon_hits_rm_overlap.tsv'
    novel_file = args.output_dir / 'exon_hits_rm_novel.tsv'
    all_file = args.output_dir / 'exon_hits_rm_distances.tsv'

    print(f"\nSaving results...")

    # Save all results with distances
    with open(all_file, 'w') as f:
        f.write('qseqid\tqstart\tqend\tsseqid\tpident\tlength\tgene_symbol\tposition\tutr_overlap\t'
                'rm_overlap\trm_distance\trm_name\trm_class\trm_divergence\n')
        for r in sorted(results, key=lambda x: x['rm_distance']):
            div_str = f"{r['rm_divergence']:.1f}" if r['rm_divergence'] is not None else ''
            f.write(f"{r['qseqid']}\t{r['qstart']}\t{r['qend']}\t{r['sseqid']}\t"
                    f"{r['pident']:.1f}\t{r['length']}\t{r['gene_symbol']}\t{r['position']}\t"
                    f"{r['utr_overlap']}\t{r['rm_overlap']}\t{r['rm_distance']}\t"
                    f"{r['rm_name'] or ''}\t{r['rm_class'] or ''}\t{div_str}\n")
    print(f"  Saved all hits: {all_file}")

    # Save overlapping hits
    overlap_results = [r for r in results if r['rm_overlap']]
    with open(overlap_file, 'w') as f:
        f.write('qseqid\tqstart\tqend\tsseqid\tpident\tlength\tgene_symbol\t'
                'rm_name\trm_class\trm_divergence\n')
        for r in sorted(overlap_results, key=lambda x: -x['pident']):
            f.write(f"{r['qseqid']}\t{r['qstart']}\t{r['qend']}\t{r['sseqid']}\t"
                    f"{r['pident']:.1f}\t{r['length']}\t{r['gene_symbol']}\t"
                    f"{r['rm_name']}\t{r['rm_class']}\t{r['rm_divergence']:.1f}\n")
    print(f"  Saved overlapping hits: {overlap_file}")

    # Save novel (non-overlapping) hits
    novel_results = [r for r in results if not r['rm_overlap']]
    with open(novel_file, 'w') as f:
        f.write('qseqid\tqstart\tqend\tsseqid\tpident\tlength\tgene_symbol\t'
                'nearest_rm_distance\tnearest_rm_name\tnearest_rm_class\n')
        for r in sorted(novel_results, key=lambda x: -x['pident'])[:50000]:
            f.write(f"{r['qseqid']}\t{r['qstart']}\t{r['qend']}\t{r['sseqid']}\t"
                    f"{r['pident']:.1f}\t{r['length']}\t{r['gene_symbol']}\t"
                    f"{r['rm_distance']}\t{r['rm_name'] or ''}\t{r['rm_class'] or ''}\n")
    print(f"  Saved novel hits (top 50k by identity): {novel_file}")

    # Generate plots
    print(f"\nGenerating distance distribution plot...")
    plot_file = args.output_dir / 'distance_to_rm_distribution.png'
    generate_distance_histogram(results, plot_file)

    # Save summary statistics as JSON
    summary = {
        'total_hits': len(results),
        'overlapping_rm': overlap_count,
        'not_overlapping_rm': non_overlap_count,
        'overlap_pct': 100 * overlap_count / len(results),
        'distance_stats': {
            'mean': mean_dist if non_overlap_distances else None,
            'median': median_dist if non_overlap_distances else None,
            'p10': p10 if non_overlap_distances else None,
            'p90': p90 if non_overlap_distances else None,
        },
        'by_rm_class': dict(rm_class_counts),
        'by_position': {k: dict(v) for k, v in position_stats.items()},
        'by_utr_overlap': {k: dict(v) for k, v in utr_stats.items()},
    }

    summary_file = args.output_dir / 'overlap_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved summary: {summary_file}")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
