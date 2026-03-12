#!/usr/bin/env python3
"""
Core motif discovery: Extract k-mers from TE-aligned sequences and compute enrichment.

Compares k-mer frequencies in real TE hits vs shuffled controls to identify
enriched/depleted motifs that may represent functional sequence elements.

Output: results/motif_analysis/motif_enrichment_6mer.tsv
"""

import argparse
import sys
from pathlib import Path
from collections import Counter
from typing import Iterator, Tuple

import numpy as np
import pandas as pd
from scipy import stats

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def extract_kmers(sequence: str, k: int = 6) -> Iterator[str]:
    """
    Extract k-mers from a sequence, removing gaps.

    Args:
        sequence: Aligned sequence (may contain '-' gaps)
        k: K-mer length (default: 6)

    Yields:
        K-mers from the ungapped sequence
    """
    # Remove gaps from aligned sequence
    clean_seq = sequence.replace('-', '').upper()

    # Skip if too short
    if len(clean_seq) < k:
        return

    # Extract all k-mers
    for i in range(len(clean_seq) - k + 1):
        kmer = clean_seq[i:i + k]
        # Skip k-mers with ambiguous bases
        if all(base in 'ACGT' for base in kmer):
            yield kmer


def count_kmers_from_blast(
    blast_file: Path,
    k: int = 6,
    use_qseq: bool = True,
    verbose: bool = False
) -> Tuple[Counter, int]:
    """
    Count k-mers from BLAST results file using streaming.

    Args:
        blast_file: Path to BLAST TSV file
        k: K-mer length
        use_qseq: If True, use query sequence; if False, use subject sequence
        verbose: Print progress

    Returns:
        Tuple of (Counter of k-mer counts, total hits processed)
    """
    kmer_counts = Counter()
    n_hits = 0
    seq_col = 'qseq' if use_qseq else 'sseq'

    for hit in iter_blast_results(blast_file):
        seq = hit.get(seq_col, '')
        if seq:
            for kmer in extract_kmers(seq, k):
                kmer_counts[kmer] += 1
        n_hits += 1

        if verbose and n_hits % 500000 == 0:
            print(f"    Processed {n_hits:,} hits, {len(kmer_counts):,} unique k-mers")

    return kmer_counts, n_hits


def compute_enrichment(
    real_counts: Counter,
    shuffled_counts_list: list,
    min_count: int = 100,
    pseudocount: float = 1.0
) -> pd.DataFrame:
    """
    Compute k-mer enrichment statistics comparing real vs shuffled.

    Args:
        real_counts: K-mer counts from real data
        shuffled_counts_list: List of Counters from shuffled replicates
        min_count: Minimum total count (real + shuffled) for inclusion
        pseudocount: Pseudocount to add for log calculations

    Returns:
        DataFrame with enrichment statistics
    """
    # Get all k-mers seen in any dataset
    all_kmers = set(real_counts.keys())
    for shuf_counts in shuffled_counts_list:
        all_kmers.update(shuf_counts.keys())

    # Compute statistics for each k-mer
    results = []
    n_reps = len(shuffled_counts_list)

    for kmer in all_kmers:
        real_count = real_counts.get(kmer, 0)

        # Get shuffled counts across replicates
        shuf_counts = [sc.get(kmer, 0) for sc in shuffled_counts_list]
        shuf_mean = np.mean(shuf_counts)
        shuf_std = np.std(shuf_counts, ddof=1) if n_reps > 1 else 0

        # Skip low-count k-mers
        total_count = real_count + sum(shuf_counts)
        if total_count < min_count:
            continue

        # Enrichment ratio (with pseudocount)
        enrichment = (real_count + pseudocount) / (shuf_mean + pseudocount)
        log2_enrichment = np.log2(enrichment)

        # Z-score (how many std devs from shuffled mean)
        if shuf_std > 0:
            z_score = (real_count - shuf_mean) / shuf_std
        else:
            # If no variance, use a simple comparison
            z_score = np.sign(real_count - shuf_mean) * 10 if real_count != shuf_mean else 0

        # P-value from z-score (two-tailed)
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))

        results.append({
            'motif': kmer,
            'real_count': real_count,
            'shuf_mean': shuf_mean,
            'shuf_std': shuf_std,
            'enrichment': enrichment,
            'log2_enrichment': log2_enrichment,
            'z_score': z_score,
            'p_value': p_value,
        })

    df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    if len(df) > 0:
        try:
            from statsmodels.stats.multitest import multipletests
            _, q_values, _, _ = multipletests(df['p_value'].values, method='fdr_bh')
            df['q_value'] = q_values
        except ImportError:
            # Fallback: Bonferroni
            df['q_value'] = np.minimum(df['p_value'] * len(df), 1.0)

    # Sort by significance
    df = df.sort_values('p_value')

    return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
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
        '-k', '--kmer-size',
        type=int,
        default=6,
        help='K-mer size (default: 6)'
    )
    parser.add_argument(
        '--min-count',
        type=int,
        default=100,
        help='Minimum total count for inclusion (default: 100)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output file (default: results/motif_analysis/motif_enrichment_6mer.tsv)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()

    # Set paths
    real_path = args.real or results_dir / "genome_wide_all_3utrs.tsv"
    shuffled_dir = args.shuffled_dir or results_dir / "shuffled_full"
    output_dir = results_dir / "motif_analysis"
    output_path = args.output or output_dir / f"motif_enrichment_{args.kmer_size}mer.tsv"

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check inputs
    if not real_path.exists():
        print(f"Error: Real BLAST file not found: {real_path}")
        return 1

    if not shuffled_dir.exists():
        print(f"Error: Shuffled directory not found: {shuffled_dir}")
        return 1

    print("=" * 70)
    print(f"TE MOTIF ENRICHMENT ANALYSIS ({args.kmer_size}-mer)")
    print("=" * 70)

    # Count k-mers in real data
    print(f"\nCounting {args.kmer_size}-mers in real data...")
    print(f"  File: {real_path}")
    real_counts, real_hits = count_kmers_from_blast(
        real_path, k=args.kmer_size, verbose=args.verbose
    )
    print(f"  Total hits: {real_hits:,}")
    print(f"  Unique {args.kmer_size}-mers: {len(real_counts):,}")
    print(f"  Total {args.kmer_size}-mers: {sum(real_counts.values()):,}")

    # Count k-mers in shuffled replicates
    print(f"\nCounting {args.kmer_size}-mers in shuffled controls...")
    shuffled_counts_list = []
    total_shuf_hits = 0

    for rep in range(1, args.n_replicates + 1):
        shuf_path = shuffled_dir / f"replicate_{rep:02d}_blast.tsv"
        if not shuf_path.exists():
            print(f"  Warning: Missing replicate {rep}")
            continue

        print(f"  Replicate {rep}...", end=" ", flush=True)
        shuf_counts, shuf_hits = count_kmers_from_blast(
            shuf_path, k=args.kmer_size, verbose=False
        )
        shuffled_counts_list.append(shuf_counts)
        total_shuf_hits += shuf_hits
        print(f"{shuf_hits:,} hits, {len(shuf_counts):,} unique k-mers")

    if not shuffled_counts_list:
        print("Error: No shuffled replicates found")
        return 1

    print(f"\n  Total shuffled hits: {total_shuf_hits:,} across {len(shuffled_counts_list)} replicates")

    # Compute enrichment
    print("\nComputing enrichment statistics...")
    df = compute_enrichment(
        real_counts,
        shuffled_counts_list,
        min_count=args.min_count
    )
    print(f"  K-mers passing filter (count >= {args.min_count}): {len(df):,}")

    # Save results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_path}")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    sig_enriched = df[(df['q_value'] < 0.05) & (df['log2_enrichment'] > 1)]
    sig_depleted = df[(df['q_value'] < 0.05) & (df['log2_enrichment'] < -1)]

    print(f"\nSignificant results (q < 0.05, |log2 enrichment| > 1):")
    print(f"  Enriched: {len(sig_enriched)} motifs")
    print(f"  Depleted: {len(sig_depleted)} motifs")

    # Top enriched motifs
    print(f"\nTop 20 enriched motifs (enrichment > 2x, q < 0.05):")
    print("-" * 70)
    top_enriched = df[(df['enrichment'] > 2) & (df['q_value'] < 0.05)].nlargest(20, 'enrichment')
    if len(top_enriched) > 0:
        print(f"{'Motif':<10} {'Real':>10} {'Shuf mean':>12} {'Enrich':>8} {'Z-score':>10} {'Q-value':>12}")
        print("-" * 70)
        for _, row in top_enriched.iterrows():
            print(f"{row['motif']:<10} {row['real_count']:>10,.0f} {row['shuf_mean']:>12,.1f} "
                  f"{row['enrichment']:>8.2f}x {row['z_score']:>10.2f} {row['q_value']:>12.2e}")
    else:
        print("  No significant enriched motifs found")

    # Top depleted motifs
    print(f"\nTop 20 depleted motifs (enrichment < 0.5x, q < 0.05):")
    print("-" * 70)
    top_depleted = df[(df['enrichment'] < 0.5) & (df['q_value'] < 0.05)].nsmallest(20, 'enrichment')
    if len(top_depleted) > 0:
        print(f"{'Motif':<10} {'Real':>10} {'Shuf mean':>12} {'Enrich':>8} {'Z-score':>10} {'Q-value':>12}")
        print("-" * 70)
        for _, row in top_depleted.iterrows():
            print(f"{row['motif']:<10} {row['real_count']:>10,.0f} {row['shuf_mean']:>12,.1f} "
                  f"{row['enrichment']:>8.2f}x {row['z_score']:>10.2f} {row['q_value']:>12.2e}")
    else:
        print("  No significant depleted motifs found")

    # Known regulatory elements check
    print("\n" + "=" * 70)
    print("KNOWN REGULATORY ELEMENTS")
    print("=" * 70)

    known_elements = {
        'AATAAA': 'Canonical poly(A) signal',
        'ATTAAA': 'Alternative poly(A) signal',
        'AGTAAA': 'Alternative poly(A) signal',
        'TATAAA': 'TATA-like motif',
        'ATTTTT': 'AU-rich element (ARE)',
        'TTTTTA': 'AU-rich element (ARE)',
        'GCGCGC': 'GC-rich element',
    }

    print(f"\n{'Motif':<10} {'Description':<30} {'Real':>10} {'Shuf':>10} {'Enrich':>8} {'Q-val':>10}")
    print("-" * 80)

    for motif, desc in known_elements.items():
        row = df[df['motif'] == motif]
        if len(row) > 0:
            row = row.iloc[0]
            print(f"{motif:<10} {desc:<30} {row['real_count']:>10,.0f} {row['shuf_mean']:>10,.1f} "
                  f"{row['enrichment']:>8.2f}x {row['q_value']:>10.2e}")
        else:
            real = real_counts.get(motif, 0)
            shuf = np.mean([sc.get(motif, 0) for sc in shuffled_counts_list])
            print(f"{motif:<10} {desc:<30} {real:>10,.0f} {shuf:>10,.1f} {'N/A':>8} {'N/A':>10}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
