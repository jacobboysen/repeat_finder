#!/usr/bin/env python3
"""
Compare real vs shuffled TE binding patterns.

When both real UTR and shuffled version hit the SAME REGION of a TE:
1. Are the match/mismatch patterns the same or different?
2. Which nucleotides are most commonly matched similarly?
3. Do they "bind" the same way or differently?
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
    """
    Extract the portion of the alignment that covers target TE positions.

    Returns: (te_seq_portion, query_seq_portion, te_positions) or None
    """
    if pd.isna(sseq) or pd.isna(qseq):
        return None

    sseq = str(sseq)
    qseq = str(qseq)

    # Determine direction
    if sstart <= send:
        # Forward: TE positions increase left to right
        te_pos = sstart
        direction = 1
    else:
        # Reverse: TE positions decrease left to right
        te_pos = sstart
        direction = -1

    # Walk through alignment, tracking TE positions
    te_portion = []
    query_portion = []
    positions = []

    for i, (s_char, q_char) in enumerate(zip(sseq, qseq)):
        # Check if this TE position is in our target range
        if target_start <= te_pos <= target_end:
            te_portion.append(s_char)
            query_portion.append(q_char)
            positions.append(te_pos)

        # Advance TE position (gaps in subject don't advance position)
        if s_char != '-':
            te_pos += direction

    if te_portion:
        return ''.join(te_portion), ''.join(query_portion), positions
    return None


def compare_alignments(te_seq1, query_seq1, te_seq2, query_seq2, positions1, positions2):
    """
    Compare two alignments to the same TE region.

    Returns dict with comparison metrics.
    """
    # Find common positions
    pos_set1 = set(positions1)
    pos_set2 = set(positions2)
    common_positions = pos_set1 & pos_set2

    if not common_positions:
        return None

    # Build position -> (te_char, query_char) maps
    map1 = {p: (te_seq1[i], query_seq1[i]) for i, p in enumerate(positions1) if p in common_positions}
    map2 = {p: (te_seq2[i], query_seq2[i]) for i, p in enumerate(positions2) if p in common_positions}

    # Compare at each common position
    same_match = 0  # Both match the TE
    same_mismatch = 0  # Both mismatch the TE
    diff_pattern = 0  # One matches, one mismatches

    # Track which TE nucleotides are commonly matched by both
    both_match_nt = Counter()  # TE nucleotide that both match
    both_mismatch_nt = Counter()  # TE nucleotide that both mismatch
    real_only_match_nt = Counter()  # TE nucleotide that only real matches
    shuf_only_match_nt = Counter()  # TE nucleotide that only shuffled matches

    # Track what query nucleotides are used
    real_query_when_both_match = Counter()
    shuf_query_when_both_match = Counter()

    for pos in common_positions:
        te1, q1 = map1[pos]
        te2, q2 = map2[pos]

        # Skip gaps in TE sequence
        if te1 == '-' or te2 == '-':
            continue

        match1 = (q1.upper() == te1.upper()) and q1 != '-'
        match2 = (q2.upper() == te2.upper()) and q2 != '-'

        te_nt = te1.upper()

        if match1 and match2:
            same_match += 1
            both_match_nt[te_nt] += 1
            real_query_when_both_match[q1.upper()] += 1
            shuf_query_when_both_match[q2.upper()] += 1
        elif not match1 and not match2:
            same_mismatch += 1
            both_mismatch_nt[te_nt] += 1
        else:
            diff_pattern += 1
            if match1:
                real_only_match_nt[te_nt] += 1
            else:
                shuf_only_match_nt[te_nt] += 1

    total = same_match + same_mismatch + diff_pattern
    if total == 0:
        return None

    return {
        'common_positions': len(common_positions),
        'same_match': same_match,
        'same_mismatch': same_mismatch,
        'diff_pattern': diff_pattern,
        'total_compared': total,
        'same_pattern_frac': (same_match + same_mismatch) / total,
        'both_match_frac': same_match / total,
        'both_match_nt': dict(both_match_nt),
        'both_mismatch_nt': dict(both_mismatch_nt),
        'real_only_match_nt': dict(real_only_match_nt),
        'shuf_only_match_nt': dict(shuf_only_match_nt),
    }


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

    # Find overlapping TE regions and compare
    print("\nComparing binding patterns for overlapping TE regions...")

    all_comparisons = []
    global_both_match_nt = Counter()
    global_both_mismatch_nt = Counter()
    global_real_only_match_nt = Counter()
    global_shuf_only_match_nt = Counter()

    n_pairs_checked = 0
    n_pairs_with_te_overlap = 0

    for pair_key in common_pairs:
        if pair_key not in real_grouped.groups or pair_key not in shuf_grouped.groups:
            continue

        n_pairs_checked += 1
        real_hits = real_grouped.get_group(pair_key)
        shuf_hits = shuf_grouped.get_group(pair_key)

        pair_has_overlap = False

        for _, real_hit in real_hits.iterrows():
            real_te_start, real_te_end = get_te_range(real_hit['sstart'], real_hit['send'])

            for _, shuf_hit in shuf_hits.iterrows():
                shuf_te_start, shuf_te_end = get_te_range(shuf_hit['sstart'], shuf_hit['send'])

                # Check for TE region overlap
                overlap = ranges_overlap(real_te_start, real_te_end, shuf_te_start, shuf_te_end)

                if overlap:
                    overlap_start, overlap_end = overlap
                    overlap_len = overlap_end - overlap_start

                    if overlap_len < 10:  # Skip tiny overlaps
                        continue

                    pair_has_overlap = True

                    # Extract alignments at overlapping positions
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
                            real_extract[0], real_extract[1],   # te_seq1, query_seq1
                            shuf_extract[0], shuf_extract[1],   # te_seq2, query_seq2
                            real_extract[2], shuf_extract[2]    # positions1, positions2
                        )

                        if comparison:
                            comparison['pair_key'] = pair_key
                            comparison['te_overlap_start'] = overlap_start
                            comparison['te_overlap_end'] = overlap_end
                            comparison['te_overlap_len'] = overlap_len
                            comparison['real_pident'] = real_hit['pident']
                            comparison['shuf_pident'] = shuf_hit['pident']
                            comparison['real_length'] = real_hit['length']
                            comparison['shuf_length'] = shuf_hit['length']
                            all_comparisons.append(comparison)

                            # Aggregate nucleotide counts
                            for nt, count in comparison['both_match_nt'].items():
                                global_both_match_nt[nt] += count
                            for nt, count in comparison['both_mismatch_nt'].items():
                                global_both_mismatch_nt[nt] += count
                            for nt, count in comparison['real_only_match_nt'].items():
                                global_real_only_match_nt[nt] += count
                            for nt, count in comparison['shuf_only_match_nt'].items():
                                global_shuf_only_match_nt[nt] += count

        if pair_has_overlap:
            n_pairs_with_te_overlap += 1

        if n_pairs_checked % 10000 == 0:
            print(f"  Checked {n_pairs_checked:,} pairs, found {len(all_comparisons):,} overlapping alignments...")

    print(f"\n  Total pairs checked: {n_pairs_checked:,}")
    print(f"  Pairs with TE region overlap (≥10bp): {n_pairs_with_te_overlap:,}")
    print(f"  Total overlapping alignment comparisons: {len(all_comparisons):,}")

    # Results
    print("\n" + "="*70)
    print("RESULTS: Binding Pattern Comparison")
    print("="*70)

    if not all_comparisons:
        print("No overlapping alignments found!")
        return

    comp_df = pd.DataFrame(all_comparisons)

    # Overall pattern similarity
    total_same_match = comp_df['same_match'].sum()
    total_same_mismatch = comp_df['same_mismatch'].sum()
    total_diff = comp_df['diff_pattern'].sum()
    total_positions = total_same_match + total_same_mismatch + total_diff

    print(f"\n--- Overall Match/Mismatch Pattern (across {total_positions:,} compared positions) ---")
    print(f"  Both match TE:      {total_same_match:,} ({100*total_same_match/total_positions:.1f}%)")
    print(f"  Both mismatch TE:   {total_same_mismatch:,} ({100*total_same_mismatch/total_positions:.1f}%)")
    print(f"  Different pattern:  {total_diff:,} ({100*total_diff/total_positions:.1f}%)")
    print(f"  Same pattern total: {total_same_match + total_same_mismatch:,} ({100*(total_same_match + total_same_mismatch)/total_positions:.1f}%)")

    # Nucleotide breakdown
    print("\n--- Which TE Nucleotides Are Most Commonly Matched by BOTH? ---")
    total_both_match = sum(global_both_match_nt.values())
    for nt in ['A', 'T', 'G', 'C']:
        count = global_both_match_nt.get(nt, 0)
        pct = 100 * count / total_both_match if total_both_match > 0 else 0
        print(f"  {nt}: {count:,} ({pct:.1f}%)")

    print("\n--- Which TE Nucleotides Are Most Commonly Mismatched by BOTH? ---")
    total_both_mismatch = sum(global_both_mismatch_nt.values())
    for nt in ['A', 'T', 'G', 'C']:
        count = global_both_mismatch_nt.get(nt, 0)
        pct = 100 * count / total_both_mismatch if total_both_mismatch > 0 else 0
        print(f"  {nt}: {count:,} ({pct:.1f}%)")

    print("\n--- Which TE Nucleotides Does ONLY Real Match (Shuffled Mismatches)? ---")
    total_real_only = sum(global_real_only_match_nt.values())
    for nt in ['A', 'T', 'G', 'C']:
        count = global_real_only_match_nt.get(nt, 0)
        pct = 100 * count / total_real_only if total_real_only > 0 else 0
        print(f"  {nt}: {count:,} ({pct:.1f}%)")

    print("\n--- Which TE Nucleotides Does ONLY Shuffled Match (Real Mismatches)? ---")
    total_shuf_only = sum(global_shuf_only_match_nt.values())
    for nt in ['A', 'T', 'G', 'C']:
        count = global_shuf_only_match_nt.get(nt, 0)
        pct = 100 * count / total_shuf_only if total_shuf_only > 0 else 0
        print(f"  {nt}: {count:,} ({pct:.1f}%)")

    # Comparison: real-only vs shuf-only
    print("\n--- Real vs Shuffled Unique Match Preference ---")
    print("  (Which nucleotides does real match that shuffled doesn't, and vice versa?)")
    for nt in ['A', 'T', 'G', 'C']:
        real_count = global_real_only_match_nt.get(nt, 0)
        shuf_count = global_shuf_only_match_nt.get(nt, 0)
        if real_count + shuf_count > 0:
            ratio = real_count / shuf_count if shuf_count > 0 else float('inf')
            print(f"  {nt}: Real-only {real_count:,} vs Shuf-only {shuf_count:,} (ratio: {ratio:.2f}x)")

    # By overlap size
    print("\n--- Pattern Similarity by TE Overlap Size ---")
    overlap_bins = [(10, 25), (25, 50), (50, 100), (100, float('inf'))]
    for lo, hi in overlap_bins:
        subset = comp_df[(comp_df['te_overlap_len'] >= lo) & (comp_df['te_overlap_len'] < hi)]
        if len(subset) > 0:
            label = f"{lo}-{int(hi)}" if hi != float('inf') else f"≥{lo}"
            sm = subset['same_match'].sum()
            smm = subset['same_mismatch'].sum()
            diff = subset['diff_pattern'].sum()
            tot = sm + smm + diff
            print(f"\n  TE overlap {label} bp: {len(subset):,} comparisons, {tot:,} positions")
            print(f"    Both match: {100*sm/tot:.1f}%, Both mismatch: {100*smm/tot:.1f}%, Different: {100*diff/tot:.1f}%")

    # Save detailed results
    # Remove dict columns for CSV
    save_df = comp_df.drop(columns=['both_match_nt', 'both_mismatch_nt', 'real_only_match_nt', 'shuf_only_match_nt'])
    output_path = results_dir / 'te_binding_pattern_comparison.tsv'
    save_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n\nDetailed results saved to: {output_path}")

    return comp_df, global_both_match_nt, global_real_only_match_nt, global_shuf_only_match_nt


if __name__ == '__main__':
    main()
