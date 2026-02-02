#!/usr/bin/env python3
"""
Analyze which TE families have better/worse real vs shuffled binding patterns.

For each TE family:
- What % of positions show "same pattern" (both match or both mismatch)?
- What % show "different pattern" (one matches, other doesn't)?
- Which families have more "real-only matches" vs "shuffled-only matches"?
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir


def get_te_range(sstart, send):
    """Return normalized (min, max) TE coordinates."""
    return min(sstart, send), max(sstart, send)


def ranges_overlap(s1, e1, s2, e2):
    """Return overlap start, end if ranges overlap, else None."""
    start1, end1 = min(s1, e1), max(s1, e1)
    start2, end2 = min(s2, e2), max(s2, e2)

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start < overlap_end:
        return overlap_start, overlap_end
    return None


def extract_alignment_at_te_positions(sstart, send, sseq, qseq, target_start, target_end):
    """Extract alignment portion covering target TE positions."""
    if pd.isna(sseq) or pd.isna(qseq):
        return None

    sseq = str(sseq)
    qseq = str(qseq)

    if sstart <= send:
        te_pos = sstart
        direction = 1
    else:
        te_pos = sstart
        direction = -1

    te_portion = []
    query_portion = []
    positions = []

    for i, (s_char, q_char) in enumerate(zip(sseq, qseq)):
        if target_start <= te_pos <= target_end:
            te_portion.append(s_char)
            query_portion.append(q_char)
            positions.append(te_pos)

        if s_char != '-':
            te_pos += direction

    if te_portion:
        return ''.join(te_portion), ''.join(query_portion), positions
    return None


def compare_alignments(te_seq1, query_seq1, te_seq2, query_seq2, positions1, positions2):
    """Compare two alignments to the same TE region."""
    pos_set1 = set(positions1)
    pos_set2 = set(positions2)
    common_positions = pos_set1 & pos_set2

    if not common_positions:
        return None

    map1 = {p: (te_seq1[i], query_seq1[i]) for i, p in enumerate(positions1) if p in common_positions}
    map2 = {p: (te_seq2[i], query_seq2[i]) for i, p in enumerate(positions2) if p in common_positions}

    same_match = 0
    same_mismatch = 0
    real_only_match = 0
    shuf_only_match = 0

    for pos in common_positions:
        te1, q1 = map1[pos]
        te2, q2 = map2[pos]

        if te1 == '-' or te2 == '-':
            continue

        match1 = (q1.upper() == te1.upper()) and q1 != '-'
        match2 = (q2.upper() == te2.upper()) and q2 != '-'

        if match1 and match2:
            same_match += 1
        elif not match1 and not match2:
            same_mismatch += 1
        elif match1:
            real_only_match += 1
        else:
            shuf_only_match += 1

    total = same_match + same_mismatch + real_only_match + shuf_only_match
    if total == 0:
        return None

    return {
        'same_match': same_match,
        'same_mismatch': same_mismatch,
        'real_only_match': real_only_match,
        'shuf_only_match': shuf_only_match,
        'total': total,
    }


def extract_te_family(te_id):
    """Extract TE family name from TE ID (e.g., FBti0019634 -> roo, etc.)."""
    # The TE ID format is FBtiXXXXXXX, we need to look up the family
    # For now, return the ID itself - we'll aggregate by the first part
    return te_id


def main():
    results_dir = get_results_dir()

    # Load data
    print("Loading real 3'UTR TE hits...")
    real_path = results_dir / 'genome_wide_all_3utrs.tsv'
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
               'qlen', 'slen', 'qseq', 'sseq', 'strand']
    real_df = pd.read_csv(real_path, sep='\t', header=None, names=columns, low_memory=False)
    real_df['transcript'] = real_df['qseqid']
    print(f"  Loaded {len(real_df):,} real hits")

    print("\nLoading shuffled control hits...")
    shuffled_path = results_dir / 'shuffled_full' / 'replicate_01_blast.tsv'
    shuf_columns = columns[:16]
    shuffled_df = pd.read_csv(shuffled_path, sep='\t', header=None, names=shuf_columns, low_memory=False)
    shuffled_df['transcript'] = shuffled_df['qseqid'].str.replace(r'_shuf\d+$', '', regex=True)
    print(f"  Loaded {len(shuffled_df):,} shuffled hits")

    # Find common transcript-TE pairs
    print("\nFinding transcript-TE pairs in both...")
    real_df['pair_key'] = real_df['transcript'] + '|' + real_df['sseqid']
    shuffled_df['pair_key'] = shuffled_df['transcript'] + '|' + shuffled_df['sseqid']

    common_pairs = set(real_df['pair_key'].unique()) & set(shuffled_df['pair_key'].unique())
    print(f"  Common pairs: {len(common_pairs):,}")

    # Filter to common pairs
    real_common = real_df[real_df['pair_key'].isin(common_pairs)].copy()
    shuf_common = shuffled_df[shuffled_df['pair_key'].isin(common_pairs)].copy()

    # Group by pair
    real_grouped = real_common.groupby('pair_key')
    shuf_grouped = shuf_common.groupby('pair_key')

    # Aggregate by TE ID
    te_stats = defaultdict(lambda: {
        'same_match': 0,
        'same_mismatch': 0,
        'real_only_match': 0,
        'shuf_only_match': 0,
        'total': 0,
        'n_comparisons': 0,
    })

    print("\nComparing binding patterns by TE...")

    n_checked = 0
    for pair_key in common_pairs:
        if pair_key not in real_grouped.groups or pair_key not in shuf_grouped.groups:
            continue

        te_id = pair_key.split('|')[1]

        real_hits = real_grouped.get_group(pair_key)
        shuf_hits = shuf_grouped.get_group(pair_key)

        for _, real_hit in real_hits.iterrows():
            real_te_start, real_te_end = get_te_range(real_hit['sstart'], real_hit['send'])

            for _, shuf_hit in shuf_hits.iterrows():
                shuf_te_start, shuf_te_end = get_te_range(shuf_hit['sstart'], shuf_hit['send'])

                overlap = ranges_overlap(real_te_start, real_te_end, shuf_te_start, shuf_te_end)

                if overlap:
                    overlap_start, overlap_end = overlap
                    overlap_len = overlap_end - overlap_start

                    if overlap_len < 10:
                        continue

                    real_extract = extract_alignment_at_te_positions(
                        real_hit['sstart'], real_hit['send'],
                        real_hit['sseq'], real_hit['qseq'],
                        overlap_start, overlap_end
                    )

                    shuf_extract = extract_alignment_at_te_positions(
                        shuf_hit['sstart'], shuf_hit['send'],
                        shuf_hit['sseq'], shuf_hit['qseq'],
                        overlap_start, overlap_end
                    )

                    if real_extract and shuf_extract:
                        comparison = compare_alignments(
                            real_extract[0], real_extract[1],
                            shuf_extract[0], shuf_extract[1],
                            real_extract[2], shuf_extract[2]
                        )

                        if comparison:
                            te_stats[te_id]['same_match'] += comparison['same_match']
                            te_stats[te_id]['same_mismatch'] += comparison['same_mismatch']
                            te_stats[te_id]['real_only_match'] += comparison['real_only_match']
                            te_stats[te_id]['shuf_only_match'] += comparison['shuf_only_match']
                            te_stats[te_id]['total'] += comparison['total']
                            te_stats[te_id]['n_comparisons'] += 1

        n_checked += 1
        if n_checked % 20000 == 0:
            print(f"  Checked {n_checked:,} pairs...")

    print(f"\n  Total TEs with comparisons: {len(te_stats):,}")

    # Convert to DataFrame and calculate metrics
    rows = []
    for te_id, stats in te_stats.items():
        if stats['total'] < 100:  # Require at least 100 compared positions
            continue

        total = stats['total']
        same_pattern = stats['same_match'] + stats['same_mismatch']
        diff_pattern = stats['real_only_match'] + stats['shuf_only_match']

        # Real advantage: does real match more often than shuffled when they differ?
        if stats['real_only_match'] + stats['shuf_only_match'] > 0:
            real_advantage = stats['real_only_match'] / (stats['real_only_match'] + stats['shuf_only_match'])
        else:
            real_advantage = 0.5

        rows.append({
            'te_id': te_id,
            'n_comparisons': stats['n_comparisons'],
            'total_positions': total,
            'same_match': stats['same_match'],
            'same_mismatch': stats['same_mismatch'],
            'real_only_match': stats['real_only_match'],
            'shuf_only_match': stats['shuf_only_match'],
            'pct_both_match': 100 * stats['same_match'] / total,
            'pct_same_pattern': 100 * same_pattern / total,
            'pct_different': 100 * diff_pattern / total,
            'real_advantage': real_advantage,
        })

    results_df = pd.DataFrame(rows)
    results_df = results_df.sort_values('total_positions', ascending=False)

    # Summary
    print("\n" + "="*80)
    print("RESULTS: TE Family Binding Pattern Analysis")
    print("="*80)

    print(f"\nTEs with ≥100 compared positions: {len(results_df):,}")

    # Overall stats
    print("\n--- Overall Statistics ---")
    print(f"  Mean % both match: {results_df['pct_both_match'].mean():.1f}%")
    print(f"  Mean % same pattern: {results_df['pct_same_pattern'].mean():.1f}%")
    print(f"  Mean % different: {results_df['pct_different'].mean():.1f}%")
    print(f"  Mean real advantage: {results_df['real_advantage'].mean():.3f}")

    # Top TEs by total positions (most data)
    print("\n--- Top 20 TEs by Data (most compared positions) ---")
    top_data = results_df.head(20)
    print(f"{'TE ID':<15} {'Positions':>10} {'%Match':>8} {'%Same':>8} {'%Diff':>8} {'RealAdv':>8}")
    print("-" * 65)
    for _, row in top_data.iterrows():
        print(f"{row['te_id']:<15} {row['total_positions']:>10,} {row['pct_both_match']:>7.1f}% "
              f"{row['pct_same_pattern']:>7.1f}% {row['pct_different']:>7.1f}% {row['real_advantage']:>7.2f}")

    # TEs where real has strongest advantage (different from shuffled)
    print("\n--- Top 20 TEs: Real Matches Better (real_advantage > 0.5) ---")
    real_better = results_df[results_df['total_positions'] >= 500].sort_values('real_advantage', ascending=False).head(20)
    print(f"{'TE ID':<15} {'Positions':>10} {'%Match':>8} {'%Diff':>8} {'RealAdv':>8}")
    print("-" * 55)
    for _, row in real_better.iterrows():
        print(f"{row['te_id']:<15} {row['total_positions']:>10,} {row['pct_both_match']:>7.1f}% "
              f"{row['pct_different']:>7.1f}% {row['real_advantage']:>7.2f}")

    # TEs where shuffled matches better
    print("\n--- Top 20 TEs: Shuffled Matches Better (real_advantage < 0.5) ---")
    shuf_better = results_df[results_df['total_positions'] >= 500].sort_values('real_advantage', ascending=True).head(20)
    print(f"{'TE ID':<15} {'Positions':>10} {'%Match':>8} {'%Diff':>8} {'RealAdv':>8}")
    print("-" * 55)
    for _, row in shuf_better.iterrows():
        print(f"{row['te_id']:<15} {row['total_positions']:>10,} {row['pct_both_match']:>7.1f}% "
              f"{row['pct_different']:>7.1f}% {row['real_advantage']:>7.2f}")

    # TEs with highest % both match (most compositionally similar)
    print("\n--- Top 20 TEs: Highest % Both Match (most compositional) ---")
    most_similar = results_df[results_df['total_positions'] >= 500].sort_values('pct_both_match', ascending=False).head(20)
    print(f"{'TE ID':<15} {'Positions':>10} {'%Match':>8} {'%Same':>8} {'RealAdv':>8}")
    print("-" * 55)
    for _, row in most_similar.iterrows():
        print(f"{row['te_id']:<15} {row['total_positions']:>10,} {row['pct_both_match']:>7.1f}% "
              f"{row['pct_same_pattern']:>7.1f}% {row['real_advantage']:>7.2f}")

    # TEs with lowest % both match (most sequence-specific)
    print("\n--- Top 20 TEs: Lowest % Both Match (most sequence-specific) ---")
    most_specific = results_df[results_df['total_positions'] >= 500].sort_values('pct_both_match', ascending=True).head(20)
    print(f"{'TE ID':<15} {'Positions':>10} {'%Match':>8} {'%Same':>8} {'RealAdv':>8}")
    print("-" * 55)
    for _, row in most_specific.iterrows():
        print(f"{row['te_id']:<15} {row['total_positions']:>10,} {row['pct_both_match']:>7.1f}% "
              f"{row['pct_same_pattern']:>7.1f}% {row['real_advantage']:>7.2f}")

    # Save results
    output_path = results_dir / 'te_family_binding_patterns.tsv'
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n\nResults saved to: {output_path}")

    return results_df


if __name__ == '__main__':
    main()
