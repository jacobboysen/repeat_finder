#!/usr/bin/env python3
"""
Compare real vs shuffled 3'UTR TE hits - V2

For cases where:
1. Real UTR A and shuffled version B both hit the same TE C
2. The hits OVERLAP on the UTR coordinates (qstart-qend overlap)

Question: Within that overlapping UTR window, are the nucleotides matching
the same positions in the TE, or different positions?
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir


def ranges_overlap(s1, e1, s2, e2):
    """Check if two ranges overlap and return overlap amount."""
    # Normalize ranges
    start1, end1 = min(s1, e1), max(s1, e1)
    start2, end2 = min(s2, e2), max(s2, e2)

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def overlap_fraction(s1, e1, s2, e2):
    """Calculate overlap as fraction of smaller range."""
    start1, end1 = min(s1, e1), max(s1, e1)
    start2, end2 = min(s2, e2), max(s2, e2)

    overlap = ranges_overlap(s1, e1, s2, e2)
    if overlap == 0:
        return 0.0

    min_len = min(end1 - start1, end2 - start2)
    return overlap / min_len if min_len > 0 else 0.0


def main():
    results_dir = get_results_dir()

    # Use raw genome-wide data for fair comparison (not deduplicated)
    print("Loading real 3'UTR TE hits (genome-wide raw)...")
    real_path = results_dir / 'genome_wide_all_3utrs.tsv'

    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
               'qlen', 'slen', 'qseq', 'sseq', 'strand']

    real_df = pd.read_csv(real_path, sep='\t', header=None, names=columns, low_memory=False)
    real_df['transcript'] = real_df['qseqid']
    print(f"  Loaded {len(real_df):,} real hits")

    # Load shuffled from shuffled_full (most recent)
    print("\nLoading shuffled control hits (shuffled_full replicate 1)...")
    shuffled_path = results_dir / 'shuffled_full' / 'replicate_01_blast.tsv'

    # Shuffled has 16 columns (no strand)
    shuf_columns = columns[:16]
    shuffled_df = pd.read_csv(shuffled_path, sep='\t', header=None, names=shuf_columns, low_memory=False)
    shuffled_df['transcript'] = shuffled_df['qseqid'].str.replace(r'_shuf\d+$', '', regex=True)
    print(f"  Loaded {len(shuffled_df):,} shuffled hits")

    # Find transcript-TE pairs in both
    print("\nFinding transcript-TE pairs in both datasets...")
    real_df['pair_key'] = real_df['transcript'] + '|' + real_df['sseqid']
    shuffled_df['pair_key'] = shuffled_df['transcript'] + '|' + shuffled_df['sseqid']

    real_pairs = set(real_df['pair_key'].unique())
    shuffled_pairs = set(shuffled_df['pair_key'].unique())
    common_pairs = real_pairs & shuffled_pairs

    print(f"  Real pairs: {len(real_pairs):,}")
    print(f"  Shuffled pairs: {len(shuffled_pairs):,}")
    print(f"  Common pairs: {len(common_pairs):,}")

    # Filter to common pairs
    real_common = real_df[real_df['pair_key'].isin(common_pairs)].copy()
    shuf_common = shuffled_df[shuffled_df['pair_key'].isin(common_pairs)].copy()

    print(f"\n  Real hits in common pairs: {len(real_common):,}")
    print(f"  Shuffled hits in common pairs: {len(shuf_common):,}")

    # Now find hits that OVERLAP on UTR coordinates (qstart-qend)
    print("\nFinding hits with overlapping UTR coordinates...")

    # Group by pair_key for comparison
    real_grouped = real_common.groupby('pair_key')
    shuf_grouped = shuf_common.groupby('pair_key')

    results = []
    n_pairs_with_utr_overlap = 0
    n_total_overlapping_hit_pairs = 0

    for pair_key in common_pairs:
        if pair_key not in real_grouped.groups or pair_key not in shuf_grouped.groups:
            continue

        real_hits = real_grouped.get_group(pair_key)
        shuf_hits = shuf_grouped.get_group(pair_key)

        pair_has_overlap = False

        # Compare all combinations of real vs shuffled hits for this pair
        for _, real_hit in real_hits.iterrows():
            for _, shuf_hit in shuf_hits.iterrows():
                # Check UTR coordinate overlap
                utr_overlap = ranges_overlap(
                    real_hit['qstart'], real_hit['qend'],
                    shuf_hit['qstart'], shuf_hit['qend']
                )

                if utr_overlap > 0:
                    pair_has_overlap = True
                    n_total_overlapping_hit_pairs += 1

                    # Calculate UTR overlap fraction
                    utr_overlap_frac = overlap_fraction(
                        real_hit['qstart'], real_hit['qend'],
                        shuf_hit['qstart'], shuf_hit['qend']
                    )

                    # Now check: within this overlapping UTR window, where do they map on TE?
                    te_overlap = ranges_overlap(
                        real_hit['sstart'], real_hit['send'],
                        shuf_hit['sstart'], shuf_hit['send']
                    )
                    te_overlap_frac = overlap_fraction(
                        real_hit['sstart'], real_hit['send'],
                        shuf_hit['sstart'], shuf_hit['send']
                    )

                    results.append({
                        'pair_key': pair_key,
                        'transcript': pair_key.split('|')[0],
                        'te_id': pair_key.split('|')[1],
                        # UTR coordinates
                        'real_utr_start': real_hit['qstart'],
                        'real_utr_end': real_hit['qend'],
                        'shuf_utr_start': shuf_hit['qstart'],
                        'shuf_utr_end': shuf_hit['qend'],
                        'utr_overlap_bp': utr_overlap,
                        'utr_overlap_frac': utr_overlap_frac,
                        # TE coordinates
                        'real_te_start': min(real_hit['sstart'], real_hit['send']),
                        'real_te_end': max(real_hit['sstart'], real_hit['send']),
                        'shuf_te_start': min(shuf_hit['sstart'], shuf_hit['send']),
                        'shuf_te_end': max(shuf_hit['sstart'], shuf_hit['send']),
                        'te_overlap_bp': te_overlap,
                        'te_overlap_frac': te_overlap_frac,
                        # Quality
                        'real_pident': real_hit['pident'],
                        'shuf_pident': shuf_hit['pident'],
                        'real_length': real_hit['length'],
                        'shuf_length': shuf_hit['length'],
                        'real_bitscore': real_hit['bitscore'],
                        'shuf_bitscore': shuf_hit['bitscore'],
                    })

        if pair_has_overlap:
            n_pairs_with_utr_overlap += 1

    results_df = pd.DataFrame(results)

    print(f"\n  Pairs with at least one UTR-overlapping hit: {n_pairs_with_utr_overlap:,}")
    print(f"  Total overlapping hit pairs: {n_total_overlapping_hit_pairs:,}")

    # Analysis
    print("\n" + "="*70)
    print("RESULTS: When real and shuffled hits OVERLAP on the UTR,")
    print("         do they match the SAME positions on the TE?")
    print("="*70)

    print(f"\nTotal UTR-overlapping hit pairs analyzed: {len(results_df):,}")

    # Categorize by TE overlap
    print("\n--- TE Position Overlap for UTR-overlapping hits ---")

    no_te_overlap = (results_df['te_overlap_bp'] == 0).sum()
    some_te_overlap = (results_df['te_overlap_bp'] > 0).sum()
    high_te_overlap = (results_df['te_overlap_frac'] > 0.5).sum()
    very_high_te_overlap = (results_df['te_overlap_frac'] > 0.8).sum()

    print(f"  No TE overlap (different TE regions): {no_te_overlap:,} ({100*no_te_overlap/len(results_df):.1f}%)")
    print(f"  Some TE overlap (>0 bp): {some_te_overlap:,} ({100*some_te_overlap/len(results_df):.1f}%)")
    print(f"  High TE overlap (>50%): {high_te_overlap:,} ({100*high_te_overlap/len(results_df):.1f}%)")
    print(f"  Very high TE overlap (>80%): {very_high_te_overlap:,} ({100*very_high_te_overlap/len(results_df):.1f}%)")

    print("\n--- Statistics ---")
    print(f"  Mean TE overlap fraction: {results_df['te_overlap_frac'].mean():.3f}")
    print(f"  Median TE overlap fraction: {results_df['te_overlap_frac'].median():.3f}")
    print(f"  Mean UTR overlap: {results_df['utr_overlap_bp'].mean():.1f} bp")

    # Breakdown by UTR overlap amount
    print("\n--- TE overlap by UTR overlap amount ---")

    utr_bins = [(0, 10), (10, 25), (25, 50), (50, 100), (100, float('inf'))]
    for lo, hi in utr_bins:
        subset = results_df[(results_df['utr_overlap_bp'] > lo) & (results_df['utr_overlap_bp'] <= hi)]
        if len(subset) > 0:
            label = f"{lo+1}-{hi}" if hi != float('inf') else f">{lo}"
            no_te = (subset['te_overlap_bp'] == 0).sum()
            some_te = (subset['te_overlap_bp'] > 0).sum()
            print(f"\n  UTR overlap {label} bp: {len(subset):,} pairs")
            print(f"    No TE overlap: {no_te} ({100*no_te/len(subset):.1f}%)")
            print(f"    Some TE overlap: {some_te} ({100*some_te/len(subset):.1f}%)")
            print(f"    Mean TE overlap frac: {subset['te_overlap_frac'].mean():.3f}")

    # For high-quality real hits
    print("\n--- High-quality real hits (≥80% identity, ≥50bp) ---")
    hq = results_df[(results_df['real_pident'] >= 80) & (results_df['real_length'] >= 50)]
    if len(hq) > 0:
        print(f"  Total HQ real hits with UTR overlap to shuffled: {len(hq):,}")
        no_te = (hq['te_overlap_bp'] == 0).sum()
        some_te = (hq['te_overlap_bp'] > 0).sum()
        print(f"    No TE overlap (different TE region): {no_te} ({100*no_te/len(hq):.1f}%)")
        print(f"    Some TE overlap (same TE region): {some_te} ({100*some_te/len(hq):.1f}%)")
        print(f"    Mean TE overlap frac: {hq['te_overlap_frac'].mean():.3f}")

    # Examples
    print("\n" + "="*70)
    print("EXAMPLES")
    print("="*70)

    # High UTR overlap + high TE overlap (same region)
    print("\n--- Same TE region (high UTR AND TE overlap) ---")
    same_region = results_df[(results_df['utr_overlap_frac'] > 0.5) & (results_df['te_overlap_frac'] > 0.5)].head(5)
    for _, row in same_region.iterrows():
        print(f"\n  {row['transcript']} -> {row['te_id']}")
        print(f"    UTR: real {row['real_utr_start']}-{row['real_utr_end']}, "
              f"shuf {row['shuf_utr_start']}-{row['shuf_utr_end']} "
              f"(overlap {row['utr_overlap_frac']:.0%})")
        print(f"    TE:  real {row['real_te_start']}-{row['real_te_end']}, "
              f"shuf {row['shuf_te_start']}-{row['shuf_te_end']} "
              f"(overlap {row['te_overlap_frac']:.0%})")
        print(f"    Quality: real {row['real_pident']:.1f}%/{row['real_length']}bp, "
              f"shuf {row['shuf_pident']:.1f}%/{row['shuf_length']}bp")

    # High UTR overlap + NO TE overlap (different region)
    print("\n--- Different TE region (high UTR overlap but NO TE overlap) ---")
    diff_region = results_df[(results_df['utr_overlap_frac'] > 0.5) & (results_df['te_overlap_bp'] == 0)].head(5)
    for _, row in diff_region.iterrows():
        print(f"\n  {row['transcript']} -> {row['te_id']}")
        print(f"    UTR: real {row['real_utr_start']}-{row['real_utr_end']}, "
              f"shuf {row['shuf_utr_start']}-{row['shuf_utr_end']} "
              f"(overlap {row['utr_overlap_frac']:.0%})")
        print(f"    TE:  real {row['real_te_start']}-{row['real_te_end']}, "
              f"shuf {row['shuf_te_start']}-{row['shuf_te_end']} "
              f"(NO overlap)")
        print(f"    Quality: real {row['real_pident']:.1f}%/{row['real_length']}bp, "
              f"shuf {row['shuf_pident']:.1f}%/{row['shuf_length']}bp")

    # Save results
    output_path = results_dir / 'real_vs_shuffled_utr_overlap_comparison.tsv'
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n\nResults saved to: {output_path}")

    return results_df


if __name__ == '__main__':
    results_df = main()
