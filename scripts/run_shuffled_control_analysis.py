#!/usr/bin/env python3
"""
Run shuffled control analysis to validate TE hit statistics.

Creates dinucleotide-shuffled versions of 3'UTRs, runs BLAST, and compares
hit statistics between real and shuffled sequences.

This establishes baseline expectations for random sequence composition effects.
"""

import argparse
import subprocess
import sys
import random
import tempfile
from pathlib import Path
from collections import defaultdict
import os

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))

def count_dinucs(seq):
    """Count dinucleotide frequencies."""
    counts = defaultdict(int)
    seq = str(seq).upper()
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if 'N' not in dinuc and all(c in 'ACGT' for c in dinuc):
            counts[dinuc] += 1
    return dict(counts)


def dinucleotide_shuffle(sequence, rng):
    """Shuffle sequence preserving dinucleotide frequencies."""
    seq = str(sequence).upper()

    if len(seq) < 2:
        return seq

    # Build edge list
    edges = []
    for i in range(len(seq) - 1):
        if seq[i] in 'ACGT' and seq[i+1] in 'ACGT':
            edges.append((seq[i], seq[i+1]))

    if not edges:
        return seq

    rng.shuffle(edges)

    # Build adjacency list
    graph = defaultdict(list)
    for src, dst in edges:
        graph[src].append(dst)

    # Find Eulerian path
    start = seq[0] if seq[0] in graph and graph[seq[0]] else next((n for n in graph if graph[n]), seq[0])

    stack = [start]
    path = []

    while stack:
        current = stack[-1]
        if graph[current]:
            stack.append(graph[current].pop())
        else:
            path.append(stack.pop())

    return ''.join(reversed(path))


def parse_fasta(filepath):
    """Simple FASTA parser."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def write_fasta(sequences, filepath):
    """Write sequences to FASTA file."""
    with open(filepath, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Write in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def run_blast(query_fasta, db_path, output_file, threads=4):
    """Run BLAST with optimized parameters."""
    cmd = [
        'blastn',
        '-query', str(query_fasta),
        '-db', str(db_path),
        '-out', str(output_file),
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-word_size', '7',
        '-gapopen', '2',
        '-gapextend', '1',
        '-penalty', '-1',
        '-reward', '1',
        '-dust', 'yes',
        '-evalue', '10',
        '-num_threads', str(threads),
        '-max_target_seqs', '500'
    ]

    subprocess.run(cmd, check=True, capture_output=True)


def parse_blast_results(filepath):
    """Parse BLAST results and return hit statistics."""
    hits = []

    if not filepath.exists() or filepath.stat().st_size == 0:
        return hits

    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                hits.append({
                    'qseqid': parts[0],
                    'sseqid': parts[1],
                    'pident': float(parts[2]),
                    'length': int(parts[3]),
                    'evalue': float(parts[10]),
                    'bitscore': float(parts[11])
                })

    return hits


def calculate_statistics(hits):
    """Calculate summary statistics for hits."""
    if not hits:
        return {
            'total_hits': 0,
            'unique_queries': 0,
            'mean_pident': 0,
            'mean_length': 0,
            'median_pident': 0,
            'median_length': 0,
            'hits_80plus': 0,
            'hits_100bp_plus': 0,
            'hits_hq': 0  # >=80% and >=50bp
        }

    pidents = [h['pident'] for h in hits]
    lengths = [h['length'] for h in hits]
    queries = set(h['qseqid'] for h in hits)

    pidents_sorted = sorted(pidents)
    lengths_sorted = sorted(lengths)
    n = len(hits)

    return {
        'total_hits': n,
        'unique_queries': len(queries),
        'mean_pident': sum(pidents) / n,
        'mean_length': sum(lengths) / n,
        'median_pident': pidents_sorted[n//2],
        'median_length': lengths_sorted[n//2],
        'hits_80plus': sum(1 for p in pidents if p >= 80),
        'hits_100bp_plus': sum(1 for l in lengths if l >= 100),
        'hits_hq': sum(1 for h in hits if h['pident'] >= 80 and h['length'] >= 50)
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--utr-fasta', type=Path,
                        default=Path('data/references/dmel_3utr.fasta'),
                        help='Input 3\'UTR FASTA file')
    parser.add_argument('--te-db', type=Path,
                        default=Path('data/blastdb/dmel_te_flybase'),
                        help='TE BLAST database')
    parser.add_argument('--output-dir', type=Path,
                        default=Path('results/shuffled_controls'),
                        help='Output directory')
    parser.add_argument('--n-shuffles', type=int, default=10,
                        help='Number of shuffle replicates')
    parser.add_argument('--sample-frac', type=float, default=0.1,
                        help='Fraction of sequences to sample (default: 0.1 = 10%%)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--threads', type=int, default=4,
                        help='BLAST threads')
    parser.add_argument('--full', action='store_true',
                        help='Use all sequences (no sampling)')

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70)
    print("SHUFFLED CONTROL ANALYSIS")
    print("="*70)
    print(f"\nInput: {args.utr_fasta}")
    print(f"TE database: {args.te_db}")
    print(f"Shuffles: {args.n_shuffles}")
    print(f"Sample fraction: {args.sample_frac if not args.full else 1.0}")
    print(f"Random seed: {args.seed}")

    # Load sequences
    print("\nLoading sequences...")
    all_sequences = parse_fasta(args.utr_fasta)
    print(f"  Loaded {len(all_sequences):,} sequences")

    # Sample if requested
    rng = random.Random(args.seed)
    if args.full:
        sample_ids = list(all_sequences.keys())
    else:
        sample_ids = rng.sample(list(all_sequences.keys()),
                                int(len(all_sequences) * args.sample_frac))

    sample_sequences = {k: all_sequences[k] for k in sample_ids}
    print(f"  Using {len(sample_sequences):,} sequences for analysis")

    # Write sampled real sequences
    real_fasta = args.output_dir / 'real_sample.fasta'
    write_fasta(sample_sequences, real_fasta)

    # Run BLAST on real sequences
    print("\n" + "-"*50)
    print("Running BLAST on real sequences...")
    real_blast = args.output_dir / 'real_blast.tsv'
    run_blast(real_fasta, args.te_db, real_blast, args.threads)
    real_hits = parse_blast_results(real_blast)
    real_stats = calculate_statistics(real_hits)
    print(f"  Real hits: {real_stats['total_hits']:,}")

    # Run shuffled replicates
    print("\n" + "-"*50)
    print(f"Running {args.n_shuffles} shuffled replicates...")

    shuffle_stats_list = []

    for i in range(args.n_shuffles):
        print(f"  Replicate {i+1}/{args.n_shuffles}...", end=' ', flush=True)

        # Create shuffled sequences with unique seed for each replicate
        rep_rng = random.Random(args.seed + i + 1)
        shuffled_seqs = {}
        for seq_id, seq in sample_sequences.items():
            shuffled_seqs[f"{seq_id}_shuf{i+1}"] = dinucleotide_shuffle(seq, rep_rng)

        # Write and BLAST
        shuf_fasta = args.output_dir / f'shuffled_rep{i+1}.fasta'
        write_fasta(shuffled_seqs, shuf_fasta)

        shuf_blast = args.output_dir / f'shuffled_rep{i+1}_blast.tsv'
        run_blast(shuf_fasta, args.te_db, shuf_blast, args.threads)

        shuf_hits = parse_blast_results(shuf_blast)
        shuf_stats = calculate_statistics(shuf_hits)
        shuffle_stats_list.append(shuf_stats)

        print(f"{shuf_stats['total_hits']:,} hits")

        # Clean up intermediate files to save space
        shuf_fasta.unlink()

    # Calculate shuffled statistics
    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)

    def mean_std(values):
        if not values:
            return 0, 0
        n = len(values)
        mean = sum(values) / n
        if n > 1:
            variance = sum((x - mean)**2 for x in values) / (n - 1)
            std = variance ** 0.5
        else:
            std = 0
        return mean, std

    # Aggregate shuffled stats
    shuf_total_hits = [s['total_hits'] for s in shuffle_stats_list]
    shuf_mean_pident = [s['mean_pident'] for s in shuffle_stats_list]
    shuf_mean_length = [s['mean_length'] for s in shuffle_stats_list]
    shuf_hq = [s['hits_hq'] for s in shuffle_stats_list]

    hits_mean, hits_std = mean_std(shuf_total_hits)
    pident_mean, pident_std = mean_std(shuf_mean_pident)
    length_mean, length_std = mean_std(shuf_mean_length)
    hq_mean, hq_std = mean_std(shuf_hq)

    print(f"\n{'Metric':<25} {'Real':>15} {'Shuffled (mean±SD)':>25} {'Fold Change':>15}")
    print("-"*80)

    fold_hits = real_stats['total_hits'] / hits_mean if hits_mean > 0 else float('inf')
    print(f"{'Total hits':<25} {real_stats['total_hits']:>15,} {f'{hits_mean:,.0f} ± {hits_std:,.0f}':>25} {fold_hits:>15.1f}x")

    fold_hq = real_stats['hits_hq'] / hq_mean if hq_mean > 0 else float('inf')
    print(f"{'HQ hits (≥80%, ≥50bp)':<25} {real_stats['hits_hq']:>15,} {f'{hq_mean:,.0f} ± {hq_std:,.0f}':>25} {fold_hq:>15.1f}x")

    print(f"{'Mean identity (%)':<25} {real_stats['mean_pident']:>15.1f} {f'{pident_mean:.1f} ± {pident_std:.1f}':>25} {'-':>15}")
    print(f"{'Mean length (bp)':<25} {real_stats['mean_length']:>15.1f} {f'{length_mean:.1f} ± {length_std:.1f}':>25} {'-':>15}")

    # Calculate z-score for total hits
    if hits_std > 0:
        z_score = (real_stats['total_hits'] - hits_mean) / hits_std
        print(f"\nZ-score (total hits): {z_score:.1f}")
        if z_score > 3:
            print("  >> Highly significant: real hits far exceed shuffled expectation")

    # Save detailed results
    results_file = args.output_dir / 'control_comparison.tsv'
    with open(results_file, 'w') as f:
        f.write("metric\treal\tshuffled_mean\tshuffled_std\tfold_change\n")
        f.write(f"total_hits\t{real_stats['total_hits']}\t{hits_mean:.0f}\t{hits_std:.0f}\t{fold_hits:.2f}\n")
        f.write(f"hits_hq\t{real_stats['hits_hq']}\t{hq_mean:.0f}\t{hq_std:.0f}\t{fold_hq:.2f}\n")
        f.write(f"mean_pident\t{real_stats['mean_pident']:.2f}\t{pident_mean:.2f}\t{pident_std:.2f}\t-\n")
        f.write(f"mean_length\t{real_stats['mean_length']:.2f}\t{length_mean:.2f}\t{length_std:.2f}\t-\n")
        f.write(f"unique_queries\t{real_stats['unique_queries']}\t{mean_std([s['unique_queries'] for s in shuffle_stats_list])[0]:.0f}\t{mean_std([s['unique_queries'] for s in shuffle_stats_list])[1]:.0f}\t-\n")

    print(f"\nResults saved to: {results_file}")

    # Save per-replicate data
    replicate_file = args.output_dir / 'shuffled_replicates.tsv'
    with open(replicate_file, 'w') as f:
        f.write("replicate\ttotal_hits\tunique_queries\tmean_pident\tmean_length\thits_hq\n")
        for i, stats in enumerate(shuffle_stats_list):
            f.write(f"{i+1}\t{stats['total_hits']}\t{stats['unique_queries']}\t{stats['mean_pident']:.2f}\t{stats['mean_length']:.2f}\t{stats['hits_hq']}\n")

    print(f"Per-replicate data: {replicate_file}")

    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)

    if fold_hits > 2:
        print(f"\nReal sequences have {fold_hits:.1f}x more TE hits than shuffled controls.")
        print("This indicates genuine TE-derived content beyond sequence composition effects.")
    else:
        print(f"\nReal sequences have only {fold_hits:.1f}x more TE hits than shuffled.")
        print("Much of the signal may be due to sequence composition rather than true TE content.")

    if fold_hq > 5:
        print(f"\nHigh-quality hits are {fold_hq:.1f}x enriched in real sequences.")
        print("Strong signal for genuine TE-derived sequences.")

    return 0


if __name__ == '__main__':
    sys.exit(main())
