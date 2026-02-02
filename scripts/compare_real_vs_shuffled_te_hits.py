#!/usr/bin/env python3
"""
Compare real vs shuffled 3'UTR TE hits to assess whether overlapping hits
(where both hit the same TE sequence) are matching the same or different
regions of the TE.

This calculates edit distance between aligned sequences for:
1. Same transcript + same TE
2. Looking at whether the TE region (sstart-send) is similar
"""

import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
import sys

# Add scripts directory to path for utils
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir

def levenshtein_distance(s1: str, s2: str) -> int:
    """Calculate Levenshtein edit distance between two strings."""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def normalized_edit_distance(s1: str, s2: str) -> float:
    """Calculate normalized edit distance (0-1 scale)."""
    if not s1 and not s2:
        return 0.0
    max_len = max(len(s1), len(s2))
    return levenshtein_distance(s1, s2) / max_len


def overlap_fraction(start1: int, end1: int, start2: int, end2: int) -> float:
    """Calculate overlap fraction between two ranges."""
    # Normalize so start < end
    s1, e1 = min(start1, end1), max(start1, end1)
    s2, e2 = min(start2, end2), max(start2, end2)

    overlap_start = max(s1, s2)
    overlap_end = min(e1, e2)

    if overlap_start >= overlap_end:
        return 0.0

    overlap_len = overlap_end - overlap_start
    min_len = min(e1 - s1, e2 - s2)

    return overlap_len / min_len if min_len > 0 else 0.0


def load_blast_hits(filepath: Path, is_shuffled: bool = False) -> pd.DataFrame:
    """Load BLAST hits from TSV file."""
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
               'qlen', 'slen', 'qseq', 'sseq']

    # Shuffled files may have 16 columns (no strand), deduplicated may have more
    df = pd.read_csv(filepath, sep='\t', header=None, names=columns[:16],
                     usecols=range(16), low_memory=False)

    if is_shuffled:
        # Extract original transcript ID (remove _shufN suffix)
        df['transcript'] = df['qseqid'].str.replace(r'_shuf\d+$', '', regex=True)
    else:
        df['transcript'] = df['qseqid']

    return df


def main():
    results_dir = get_results_dir()

    # Load real hits (deduplicated for higher quality)
    print("Loading real 3'UTR TE hits...")
    real_path = results_dir / '3utr_deduplicated' / '3utr_deduplicated_hits.tsv'
    real_df = load_blast_hits(real_path, is_shuffled=False)
    print(f"  Loaded {len(real_df):,} real hits")

    # Load shuffled hits (use replicate 1 as representative)
    print("\nLoading shuffled control hits (replicate 1)...")
    shuffled_path = results_dir / 'shuffled_full' / 'replicate_01_blast.tsv'
    shuffled_df = load_blast_hits(shuffled_path, is_shuffled=True)
    print(f"  Loaded {len(shuffled_df):,} shuffled hits")

    # Find overlapping transcript-TE pairs
    print("\nFinding transcript-TE pairs present in both real and shuffled...")

    # Create keys for matching: (transcript, TE)
    real_df['pair_key'] = real_df['transcript'] + '|' + real_df['sseqid']
    shuffled_df['pair_key'] = shuffled_df['transcript'] + '|' + shuffled_df['sseqid']

    real_pairs = set(real_df['pair_key'].unique())
    shuffled_pairs = set(shuffled_df['pair_key'].unique())

    overlapping_pairs = real_pairs & shuffled_pairs
    print(f"  Real unique transcript-TE pairs: {len(real_pairs):,}")
    print(f"  Shuffled unique transcript-TE pairs: {len(shuffled_pairs):,}")
    print(f"  Overlapping pairs: {len(overlapping_pairs):,}")

    # Filter to overlapping pairs
    real_overlap = real_df[real_df['pair_key'].isin(overlapping_pairs)].copy()
    shuffled_overlap = shuffled_df[shuffled_df['pair_key'].isin(overlapping_pairs)].copy()

    print(f"\n  Real hits in overlapping pairs: {len(real_overlap):,}")
    print(f"  Shuffled hits in overlapping pairs: {len(shuffled_overlap):,}")

    # For each overlapping pair, compare the best hit from each
    print("\nComparing edit distances for overlapping pairs...")

    # Get best hit per pair (by bitscore)
    real_best = real_overlap.loc[real_overlap.groupby('pair_key')['bitscore'].idxmax()]
    shuffled_best = shuffled_overlap.loc[shuffled_overlap.groupby('pair_key')['bitscore'].idxmax()]

    # Merge for comparison
    comparison = real_best.merge(
        shuffled_best,
        on='pair_key',
        suffixes=('_real', '_shuf')
    )

    print(f"\n  Comparing {len(comparison):,} transcript-TE pairs...")

    # Calculate metrics
    results = []
    for idx, row in comparison.iterrows():
        # TE region overlap
        te_overlap = overlap_fraction(
            row['sstart_real'], row['send_real'],
            row['sstart_shuf'], row['send_shuf']
        )

        # Only calculate edit distance for sequence comparison if there's some overlap
        # or if sequences are similar length (meaningful comparison)
        seq_real = str(row['qseq_real']).replace('-', '') if pd.notna(row['qseq_real']) else ''
        seq_shuf = str(row['qseq_shuf']).replace('-', '') if pd.notna(row['qseq_shuf']) else ''

        # Aligned sequence comparison (with gaps)
        aligned_real = str(row['sseq_real']) if pd.notna(row['sseq_real']) else ''
        aligned_shuf = str(row['sseq_shuf']) if pd.notna(row['sseq_shuf']) else ''

        # Calculate edit distance on the TE-aligned sequences
        if len(aligned_real) > 0 and len(aligned_shuf) > 0:
            norm_edit_dist = normalized_edit_distance(aligned_real, aligned_shuf)
        else:
            norm_edit_dist = 1.0

        results.append({
            'pair_key': row['pair_key'],
            'transcript': row['pair_key'].split('|')[0],
            'te_id': row['pair_key'].split('|')[1],
            'real_pident': row['pident_real'],
            'shuf_pident': row['pident_shuf'],
            'real_length': row['length_real'],
            'shuf_length': row['length_shuf'],
            'real_bitscore': row['bitscore_real'],
            'shuf_bitscore': row['bitscore_shuf'],
            'real_te_start': min(row['sstart_real'], row['send_real']),
            'real_te_end': max(row['sstart_real'], row['send_real']),
            'shuf_te_start': min(row['sstart_shuf'], row['send_shuf']),
            'shuf_te_end': max(row['sstart_shuf'], row['send_shuf']),
            'te_region_overlap': te_overlap,
            'norm_edit_distance': norm_edit_dist,
        })

    results_df = pd.DataFrame(results)

    # Categorize the comparisons
    def categorize(row):
        if row['te_region_overlap'] > 0.5 and row['norm_edit_distance'] < 0.3:
            return 'same_match'
        elif row['te_region_overlap'] > 0.3:
            return 'close_match'
        else:
            return 'different_match'

    results_df['match_category'] = results_df.apply(categorize, axis=1)

    # Summary statistics
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)

    print(f"\nTotal overlapping transcript-TE pairs analyzed: {len(results_df):,}")

    print("\n--- Match Categories ---")
    category_counts = results_df['match_category'].value_counts()
    for cat, count in category_counts.items():
        pct = 100 * count / len(results_df)
        print(f"  {cat}: {count:,} ({pct:.1f}%)")

    print("\n--- TE Region Overlap Statistics ---")
    print(f"  Mean overlap: {results_df['te_region_overlap'].mean():.3f}")
    print(f"  Median overlap: {results_df['te_region_overlap'].median():.3f}")
    print(f"  Pairs with >50% overlap: {(results_df['te_region_overlap'] > 0.5).sum():,}")
    print(f"  Pairs with >80% overlap: {(results_df['te_region_overlap'] > 0.8).sum():,}")
    print(f"  Pairs with 0% overlap: {(results_df['te_region_overlap'] == 0).sum():,}")

    print("\n--- Edit Distance Statistics ---")
    print(f"  Mean normalized edit distance: {results_df['norm_edit_distance'].mean():.3f}")
    print(f"  Median normalized edit distance: {results_df['norm_edit_distance'].median():.3f}")
    print(f"  Pairs with <0.2 edit distance (very similar): {(results_df['norm_edit_distance'] < 0.2).sum():,}")
    print(f"  Pairs with <0.3 edit distance (similar): {(results_df['norm_edit_distance'] < 0.3).sum():,}")

    print("\n--- Quality Comparison (Real vs Shuffled) ---")
    print(f"  Mean real pident: {results_df['real_pident'].mean():.1f}%")
    print(f"  Mean shuffled pident: {results_df['shuf_pident'].mean():.1f}%")
    print(f"  Mean real length: {results_df['real_length'].mean():.1f} bp")
    print(f"  Mean shuffled length: {results_df['shuf_length'].mean():.1f} bp")
    print(f"  Mean real bitscore: {results_df['real_bitscore'].mean():.1f}")
    print(f"  Mean shuffled bitscore: {results_df['shuf_bitscore'].mean():.1f}")

    # Breakdown by overlap category
    print("\n--- Detailed Breakdown by TE Region Overlap ---")
    overlap_bins = [0, 0.01, 0.25, 0.5, 0.75, 1.0]
    overlap_labels = ['0%', '1-25%', '25-50%', '50-75%', '75-100%']
    results_df['overlap_bin'] = pd.cut(results_df['te_region_overlap'],
                                        bins=overlap_bins,
                                        labels=overlap_labels,
                                        include_lowest=True)

    for bin_label in overlap_labels:
        subset = results_df[results_df['overlap_bin'] == bin_label]
        if len(subset) > 0:
            print(f"\n  TE Region Overlap: {bin_label}")
            print(f"    Count: {len(subset):,} ({100*len(subset)/len(results_df):.1f}%)")
            print(f"    Mean edit distance: {subset['norm_edit_distance'].mean():.3f}")
            print(f"    Mean real pident: {subset['real_pident'].mean():.1f}%")
            print(f"    Mean shuffled pident: {subset['shuf_pident'].mean():.1f}%")

    # Save results
    output_path = results_dir / 'real_vs_shuffled_comparison.tsv'
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\nDetailed results saved to: {output_path}")

    # Show some examples
    print("\n" + "="*70)
    print("EXAMPLE COMPARISONS")
    print("="*70)

    print("\n--- Same Match Examples (high overlap, low edit distance) ---")
    same_matches = results_df[results_df['match_category'] == 'same_match'].head(5)
    for _, row in same_matches.iterrows():
        print(f"\n  {row['transcript']} -> {row['te_id']}")
        print(f"    TE region: real {row['real_te_start']}-{row['real_te_end']}, "
              f"shuf {row['shuf_te_start']}-{row['shuf_te_end']}")
        print(f"    Overlap: {row['te_region_overlap']:.1%}, Edit dist: {row['norm_edit_distance']:.3f}")
        print(f"    Real: {row['real_pident']:.1f}% identity, {row['real_length']}bp")
        print(f"    Shuf: {row['shuf_pident']:.1f}% identity, {row['shuf_length']}bp")

    print("\n--- Different Match Examples (low overlap) ---")
    diff_matches = results_df[results_df['match_category'] == 'different_match'].head(5)
    for _, row in diff_matches.iterrows():
        print(f"\n  {row['transcript']} -> {row['te_id']}")
        print(f"    TE region: real {row['real_te_start']}-{row['real_te_end']}, "
              f"shuf {row['shuf_te_start']}-{row['shuf_te_end']}")
        print(f"    Overlap: {row['te_region_overlap']:.1%}, Edit dist: {row['norm_edit_distance']:.3f}")
        print(f"    Real: {row['real_pident']:.1f}% identity, {row['real_length']}bp")
        print(f"    Shuf: {row['shuf_pident']:.1f}% identity, {row['shuf_length']}bp")

    return results_df


if __name__ == '__main__':
    results_df = main()
