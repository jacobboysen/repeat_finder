#!/usr/bin/env python3
"""
Compare conservation of TE hits that contain specific motifs vs those that don't.

Uses existing pre-computed conservation data to avoid bigWig queries.

Approach:
1. Load TE hit sequences (qseq from BLAST results)
2. Load pre-computed conservation scores for TE hits
3. For each motif, identify which TE hits contain it
4. Compare conservation: hits WITH motif vs hits WITHOUT motif
5. Also compare: motifs in TE hits vs overall UTR conservation
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def load_conservation_scores(tab_path: Path) -> dict:
    """
    Load conservation scores from bigWigAverageOverBed output.
    Returns: {hit_name: phylop_mean}
    """
    scores = {}
    with open(tab_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                name = parts[0]
                mean_score = float(parts[5])
                scores[name] = mean_score
    return scores


def sequence_contains_motif(sequence: str, motif: str) -> bool:
    """Check if sequence contains motif (after removing gaps)."""
    clean_seq = sequence.replace('-', '').upper()
    return motif.upper() in clean_seq


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--motifs', type=str,
                        default='AATAAA,ATTAAA,TGTATA,TGTAAA,TATTTA,CAGCAG,AGCAGC,ACACAC,CACACA,GGTAAG',
                        help='Comma-separated list of motifs to analyze')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_conservation_by_te_hit.tsv"

    # File paths
    blast_path = results_dir / "genome_wide_all_3utrs.tsv"
    conservation_path = results_dir / "repeatmasker_analysis" / "te_hits_all_conservation.tab"
    utr_conservation_path = results_dir / "repeatmasker_analysis" / "utrs_conservation.tab"

    motifs = [m.strip().upper() for m in args.motifs.split(',')]

    print("=" * 70)
    print("MOTIF CONSERVATION ANALYSIS (using pre-computed scores)")
    print("=" * 70)
    print(f"Motifs: {', '.join(motifs)}")

    # Load conservation scores
    print("\n1. Loading conservation scores...")
    conservation = load_conservation_scores(conservation_path)
    print(f"   Loaded scores for {len(conservation):,} TE hits")

    # Load UTR background conservation
    print("\n2. Loading UTR background conservation...")
    utr_conservation = load_conservation_scores(utr_conservation_path)
    utr_mean = np.mean(list(utr_conservation.values()))
    print(f"   Loaded {len(utr_conservation):,} UTRs, mean phyloP: {utr_mean:.3f}")

    # Build hit name to sequence mapping and check motif presence
    print("\n3. Loading BLAST results and checking motif presence...")

    # We need to build the same hit names as in conservation file
    # Format: novel_ID|transcript|TE|pident|length or known_ID|...

    hit_data = {}  # hit_name -> {sequence, phylop, has_motif_X, ...}
    n_processed = 0

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        sseqid = hit.get('sseqid', '')
        pident = hit.get('pident', 0)
        length = hit.get('length', 0)
        qseq = hit.get('qseq', '')

        if not qseq:
            continue

        n_processed += 1
        if n_processed % 500000 == 0:
            print(f"   Processed {n_processed:,} hits...")

        # Try to match conservation file naming
        # Conservation file has: novel_ID|transcript|TE|pident|length
        # We need to reconstruct this or use a simpler matching approach

        # For now, store by a composite key
        # Note: conservation file uses sequential IDs, so exact matching is tricky
        # Let's store what we can and match by transcript/position later

        clean_seq = qseq.replace('-', '').upper()

        record = {
            'transcript': qseqid,
            'te': sseqid,
            'pident': pident,
            'length': length,
            'sequence': clean_seq,
        }

        # Check each motif
        for motif in motifs:
            record[f'has_{motif}'] = motif in clean_seq

        hit_data[n_processed] = record

    print(f"   Total hits processed: {n_processed:,}")

    # Match conservation scores to hits
    # The conservation file uses format: category_idx|transcript|te|pident|length
    print("\n4. Matching conservation scores to hits...")

    # Parse conservation names to extract info
    cons_by_key = {}
    for name, score in conservation.items():
        # Parse: novel_1234|FBtr0123|FBti0456|80.5|50
        parts = name.split('|')
        if len(parts) >= 5:
            transcript = parts[1]
            te = parts[2]
            try:
                pident = float(parts[3])
                length = int(parts[4])
                key = (transcript, te, round(pident, 1), length)
                cons_by_key[key] = score
            except:
                continue

    print(f"   Parsed {len(cons_by_key):,} conservation entries")

    # Match hits to conservation
    matched = 0
    for idx, record in hit_data.items():
        key = (record['transcript'], record['te'],
               round(record['pident'], 1), record['length'])
        if key in cons_by_key:
            record['phylop'] = cons_by_key[key]
            matched += 1
        else:
            record['phylop'] = np.nan

    print(f"   Matched conservation for {matched:,} / {len(hit_data):,} hits ({100*matched/len(hit_data):.1f}%)")

    # Convert to DataFrame
    df = pd.DataFrame(list(hit_data.values()))
    df = df.dropna(subset=['phylop'])
    print(f"\n   Hits with conservation data: {len(df):,}")

    # Analyze conservation by motif presence
    print("\n" + "=" * 70)
    print("RESULTS: Conservation of TE Hits by Motif Presence")
    print("=" * 70)

    print(f"\nBackground: Mean UTR phyloP = {utr_mean:.3f}")
    print(f"            Mean TE hit phyloP = {df['phylop'].mean():.3f}")

    results = []

    for motif in motifs:
        col = f'has_{motif}'
        if col not in df.columns:
            continue

        with_motif = df[df[col] == True]
        without_motif = df[df[col] == False]

        if len(with_motif) < 10:
            continue

        with_mean = with_motif['phylop'].mean()
        without_mean = without_motif['phylop'].mean()
        diff = with_mean - without_mean

        with_conserved = (with_motif['phylop'] > 1).mean() * 100
        without_conserved = (without_motif['phylop'] > 1).mean() * 100

        results.append({
            'motif': motif,
            'n_with': len(with_motif),
            'n_without': len(without_motif),
            'pct_with': len(with_motif) / len(df) * 100,
            'mean_with': with_mean,
            'mean_without': without_mean,
            'diff': diff,
            'pct_cons_with': with_conserved,
            'pct_cons_without': without_conserved,
        })

    results_df = pd.DataFrame(results)

    # Print results
    print(f"\n{'Motif':<8} {'N with':>10} {'N w/o':>10} {'%hits':>7} {'Mean(w)':>9} {'Mean(w/o)':>10} {'Diff':>8} {'%Cons(w)':>9} {'%Cons(w/o)':>10}")
    print("-" * 100)

    for _, row in results_df.iterrows():
        diff_str = f"{row['diff']:+.4f}"
        print(f"{row['motif']:<8} {row['n_with']:>10,} {row['n_without']:>10,} "
              f"{row['pct_with']:>6.1f}% {row['mean_with']:>9.4f} {row['mean_without']:>10.4f} "
              f"{diff_str:>8} {row['pct_cons_with']:>8.1f}% {row['pct_cons_without']:>9.1f}%")

    # Save results
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\nResults saved to: {output_path}")

    # Statistical summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    more_conserved = results_df[results_df['diff'] > 0]
    less_conserved = results_df[results_df['diff'] < 0]

    print(f"\nMotifs where TE hits WITH motif are MORE conserved: {len(more_conserved)}")
    if len(more_conserved) > 0:
        for _, row in more_conserved.iterrows():
            print(f"  {row['motif']}: +{row['diff']:.4f} phyloP")

    print(f"\nMotifs where TE hits WITH motif are LESS conserved: {len(less_conserved)}")
    if len(less_conserved) > 0:
        for _, row in less_conserved.iterrows():
            print(f"  {row['motif']}: {row['diff']:.4f} phyloP")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
This analysis compares conservation of TE hits that CONTAIN a motif
vs TE hits that DON'T contain the motif.

'Diff' = Mean phyloP (with motif) - Mean phyloP (without motif)
  Positive = TE hits containing this motif are MORE conserved
  Negative = TE hits containing this motif are LESS conserved

If regulatory motifs in TE hits are functional, we expect:
  - Functional motifs to show POSITIVE diff (more conserved)
  - Non-functional/recent motifs to show ~0 or negative diff
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
