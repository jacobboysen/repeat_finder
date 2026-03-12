#!/usr/bin/env python3
"""
Analyze positional distribution of enriched motifs within 3'UTRs.

Maps where specific motifs occur within UTRs (5' end vs 3' end) and
compares real vs shuffled distributions.

Output: results/motif_analysis/motif_position_density.tsv
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple

import numpy as np
import pandas as pd

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def find_motif_positions(sequence: str, motif: str) -> List[int]:
    """
    Find all positions of a motif in a sequence (after removing gaps).

    Args:
        sequence: Aligned sequence (may contain '-' gaps)
        motif: Motif to search for

    Returns:
        List of start positions (0-indexed) in the ungapped sequence
    """
    clean_seq = sequence.replace('-', '').upper()
    motif = motif.upper()

    positions = []
    start = 0
    while True:
        pos = clean_seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1

    return positions


def compute_relative_positions(
    blast_file: Path,
    motifs: List[str],
    n_bins: int = 10
) -> Tuple[dict, dict]:
    """
    Compute relative position distributions for motifs.

    Args:
        blast_file: Path to BLAST TSV file
        motifs: List of motifs to analyze
        n_bins: Number of bins for distribution (default: 10 for deciles)

    Returns:
        Tuple of:
        - dict mapping motif -> array of position counts per bin
        - dict mapping motif -> total occurrences
    """
    motif_positions = {m: [] for m in motifs}
    motif_counts = {m: 0 for m in motifs}

    for hit in iter_blast_results(blast_file):
        qseq = hit.get('qseq', '')
        qstart = hit.get('qstart', 0)
        qlen = hit.get('qlen', 0)

        if not qseq or qlen == 0:
            continue

        # Remove gaps and get ungapped length
        clean_seq = qseq.replace('-', '').upper()

        for motif in motifs:
            positions = find_motif_positions(qseq, motif)
            for pos in positions:
                # Compute absolute position in UTR
                # Position in alignment + qstart gives UTR position
                # Note: pos is position in ungapped aligned region
                # We need to map this to the original UTR position

                # Approximate: fraction through aligned region × alignment span
                align_len = len(clean_seq)
                if align_len > 0:
                    frac_through_align = pos / align_len
                    utr_pos = qstart + frac_through_align * (hit.get('qend', qstart) - qstart)

                    # Relative position in UTR (0 = 5' end, 1 = 3' end)
                    rel_pos = utr_pos / qlen
                    rel_pos = max(0, min(1, rel_pos))  # Clamp to [0, 1]

                    motif_positions[motif].append(rel_pos)
                    motif_counts[motif] += 1

    # Bin the positions into deciles
    bin_edges = np.linspace(0, 1, n_bins + 1)
    position_hist = {}

    for motif in motifs:
        if motif_positions[motif]:
            counts, _ = np.histogram(motif_positions[motif], bins=bin_edges)
            position_hist[motif] = counts
        else:
            position_hist[motif] = np.zeros(n_bins, dtype=int)

    return position_hist, motif_counts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--motif-results',
        type=Path,
        default=None,
        help='Motif enrichment results (default: results/motif_analysis/motif_enrichment_6mer.tsv)'
    )
    parser.add_argument(
        '--real',
        type=Path,
        default=None,
        help='Real BLAST results (default: results/genome_wide_all_3utrs.tsv)'
    )
    parser.add_argument(
        '--shuffled-dir',
        type=Path,
        default=None,
        help='Directory with shuffled BLAST results (default: results/shuffled_full)'
    )
    parser.add_argument(
        '--n-replicates',
        type=int,
        default=10,
        help='Number of shuffled replicates to use (default: 10)'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=50,
        help='Number of top enriched/depleted motifs to analyze (default: 50)'
    )
    parser.add_argument(
        '--n-bins',
        type=int,
        default=10,
        help='Number of position bins/deciles (default: 10)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output file (default: results/motif_analysis/motif_position_density.tsv)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    motif_dir = results_dir / "motif_analysis"

    # Set paths
    motif_results_path = args.motif_results or motif_dir / "motif_enrichment_6mer.tsv"
    real_path = args.real or results_dir / "genome_wide_all_3utrs.tsv"
    shuffled_dir = args.shuffled_dir or results_dir / "shuffled_full"
    output_path = args.output or motif_dir / "motif_position_density.tsv"

    # Check inputs
    if not motif_results_path.exists():
        print(f"Error: Motif results not found: {motif_results_path}")
        print("Run analyze_te_motifs.py first.")
        return 1

    if not real_path.exists():
        print(f"Error: Real BLAST file not found: {real_path}")
        return 1

    print("=" * 70)
    print("MOTIF POSITIONAL DISTRIBUTION ANALYSIS")
    print("=" * 70)

    # Load motif enrichment results
    print(f"\nLoading motif results from: {motif_results_path}")
    motif_df = pd.read_csv(motif_results_path, sep='\t')
    print(f"  Total motifs: {len(motif_df)}")

    # Select top enriched and depleted motifs
    sig_motifs = motif_df[motif_df['q_value'] < 0.05]
    top_enriched = sig_motifs.nlargest(args.top_n, 'enrichment')['motif'].tolist()
    top_depleted = sig_motifs.nsmallest(args.top_n, 'enrichment')['motif'].tolist()

    # Also include known regulatory motifs
    known_motifs = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'ATTTTT', 'TTTTTA']
    all_motifs = list(set(top_enriched + top_depleted + known_motifs))
    print(f"  Analyzing {len(all_motifs)} motifs ({len(top_enriched)} enriched, {len(top_depleted)} depleted, known)")

    # Compute positions in real data
    print(f"\nAnalyzing positions in real data...")
    real_hist, real_counts = compute_relative_positions(real_path, all_motifs, args.n_bins)

    total_real = sum(real_counts.values())
    print(f"  Total motif occurrences: {total_real:,}")

    # Compute positions in shuffled data
    print(f"\nAnalyzing positions in shuffled controls...")
    shuf_hists = []

    for rep in range(1, args.n_replicates + 1):
        shuf_path = shuffled_dir / f"replicate_{rep:02d}_blast.tsv"
        if not shuf_path.exists():
            print(f"  Warning: Missing replicate {rep}")
            continue

        print(f"  Replicate {rep}...", end=" ", flush=True)
        shuf_hist, shuf_counts = compute_relative_positions(shuf_path, all_motifs, args.n_bins)
        shuf_hists.append(shuf_hist)
        print(f"{sum(shuf_counts.values()):,} occurrences")

    if not shuf_hists:
        print("Error: No shuffled replicates processed")
        return 1

    # Compile results
    print("\nCompiling results...")
    results = []

    for motif in all_motifs:
        # Get enrichment info
        motif_info = motif_df[motif_df['motif'] == motif]
        if len(motif_info) > 0:
            enrichment = motif_info.iloc[0]['enrichment']
            q_value = motif_info.iloc[0]['q_value']
        else:
            enrichment = 1.0
            q_value = 1.0

        # Real distribution
        real_dist = real_hist.get(motif, np.zeros(args.n_bins))
        real_total = real_counts.get(motif, 0)

        # Shuffled distributions (averaged)
        shuf_dists = [h.get(motif, np.zeros(args.n_bins)) for h in shuf_hists]
        shuf_mean_dist = np.mean(shuf_dists, axis=0)
        shuf_std_dist = np.std(shuf_dists, axis=0)
        shuf_total = np.mean([sum(h.get(motif, np.zeros(args.n_bins))) for h in shuf_hists])

        # Normalize to density
        if real_total > 0:
            real_density = real_dist / real_total
        else:
            real_density = np.zeros(args.n_bins)

        if shuf_total > 0:
            shuf_density = shuf_mean_dist / shuf_total
        else:
            shuf_density = np.zeros(args.n_bins)

        # Add row for each bin
        for i in range(args.n_bins):
            bin_start = i / args.n_bins
            bin_end = (i + 1) / args.n_bins

            results.append({
                'motif': motif,
                'enrichment': enrichment,
                'q_value': q_value,
                'bin': i + 1,
                'bin_start': bin_start,
                'bin_end': bin_end,
                'real_count': int(real_dist[i]),
                'real_density': real_density[i],
                'shuf_mean_count': shuf_mean_dist[i],
                'shuf_std_count': shuf_std_dist[i],
                'shuf_density': shuf_density[i],
                'density_ratio': real_density[i] / shuf_density[i] if shuf_density[i] > 0 else np.nan,
            })

    # Create DataFrame
    df = pd.DataFrame(results)

    # Save results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("POSITIONAL BIAS SUMMARY")
    print("=" * 70)

    # Analyze positional bias for top motifs
    print("\nMotifs with strong 5' bias (enriched in bins 1-3):")
    print("-" * 50)

    motif_bias = []
    for motif in all_motifs:
        motif_data = df[df['motif'] == motif]
        if motif_data['real_count'].sum() < 100:
            continue

        five_prime = motif_data[motif_data['bin'] <= 3]['real_density'].sum()
        three_prime = motif_data[motif_data['bin'] >= 8]['real_density'].sum()

        if five_prime + three_prime > 0:
            bias_ratio = five_prime / (five_prime + three_prime)
            motif_bias.append({
                'motif': motif,
                'five_prime_density': five_prime,
                'three_prime_density': three_prime,
                'bias_ratio': bias_ratio,  # >0.5 = 5' biased, <0.5 = 3' biased
                'total': motif_data['real_count'].sum()
            })

    bias_df = pd.DataFrame(motif_bias)
    if len(bias_df) > 0:
        # 5' biased motifs
        five_biased = bias_df[bias_df['bias_ratio'] > 0.6].nlargest(10, 'bias_ratio')
        if len(five_biased) > 0:
            hdr_5p = "5' density"
            hdr_3p = "3' density"
            print(f"{'Motif':<10} {hdr_5p:>12} {hdr_3p:>12} {'Bias ratio':>12}")
            print("-" * 50)
            for _, row in five_biased.iterrows():
                print(f"{row['motif']:<10} {row['five_prime_density']:>12.3f} "
                      f"{row['three_prime_density']:>12.3f} {row['bias_ratio']:>12.3f}")

        # 3' biased motifs
        print("\nMotifs with strong 3' bias (enriched in bins 8-10):")
        print("-" * 50)
        three_biased = bias_df[bias_df['bias_ratio'] < 0.4].nsmallest(10, 'bias_ratio')
        if len(three_biased) > 0:
            hdr_5p = "5' density"
            hdr_3p = "3' density"
            print(f"{'Motif':<10} {hdr_5p:>12} {hdr_3p:>12} {'Bias ratio':>12}")
            print("-" * 50)
            for _, row in three_biased.iterrows():
                print(f"{row['motif']:<10} {row['five_prime_density']:>12.3f} "
                      f"{row['three_prime_density']:>12.3f} {row['bias_ratio']:>12.3f}")

    # Check known regulatory elements
    print("\n" + "=" * 70)
    print("KNOWN REGULATORY ELEMENT POSITIONS")
    print("=" * 70)

    known_info = {
        'AATAAA': 'Poly(A) signal - expected near 3\' end',
        'ATTAAA': 'Alt poly(A) signal - expected near 3\' end',
        'ATTTTT': 'AU-rich element - variable position',
    }

    for motif, desc in known_info.items():
        motif_data = df[df['motif'] == motif]
        if len(motif_data) > 0 and motif_data['real_count'].sum() > 0:
            total = motif_data['real_count'].sum()
            peak_bin = motif_data.loc[motif_data['real_density'].idxmax()]['bin']
            five_prime = motif_data[motif_data['bin'] <= 3]['real_density'].sum()
            three_prime = motif_data[motif_data['bin'] >= 8]['real_density'].sum()

            print(f"\n{motif} - {desc}")
            print(f"  Total occurrences: {total:,}")
            print(f"  Peak position: bin {peak_bin} ({(peak_bin-1)*10}-{peak_bin*10}%)")
            print(f"  5' (bins 1-3): {five_prime*100:.1f}% of density")
            print(f"  3' (bins 8-10): {three_prime*100:.1f}% of density")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
