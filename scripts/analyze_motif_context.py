#!/usr/bin/env python3
"""
Compare motif frequencies across four contexts:
1. Full 3'UTR sequences (background)
2. TE-aligned regions within 3'UTRs (real hits)
3. Full shuffled 3'UTR sequences (shuffled background)
4. TE-aligned regions within shuffled 3'UTRs (shuffled hits)

This reveals whether motifs are enriched in TEs specifically vs just common in 3'UTRs.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter
from itertools import product

import pandas as pd
from Bio import SeqIO

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def count_kmers_in_sequence(sequence: str, k: int = 6) -> Counter:
    """Count all k-mers in a sequence."""
    sequence = sequence.upper().replace('-', '')
    counts = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:  # Skip k-mers with N
            counts[kmer] += 1
    return counts


def count_kmers_in_fasta(fasta_path: Path, k: int = 6) -> Counter:
    """Count all k-mers across all sequences in a FASTA file."""
    total_counts = Counter()
    n_seqs = 0
    total_bp = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        total_bp += len(seq)
        n_seqs += 1
        counts = count_kmers_in_sequence(seq, k)
        total_counts.update(counts)

    return total_counts, n_seqs, total_bp


def count_kmers_in_blast_hits(blast_path: Path, k: int = 6) -> Counter:
    """Count all k-mers in the qseq column of BLAST results."""
    total_counts = Counter()
    n_hits = 0
    total_bp = 0

    for hit in iter_blast_results(blast_path):
        qseq = hit.get('qseq', '')
        if qseq:
            clean_seq = qseq.replace('-', '')
            total_bp += len(clean_seq)
            n_hits += 1
            counts = count_kmers_in_sequence(qseq, k)
            total_counts.update(counts)

    return total_counts, n_hits, total_bp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k', type=int, default=6, help='K-mer size (default: 6)')
    parser.add_argument('--top-n', type=int, default=100, help='Number of top motifs to show (default: 100)')
    parser.add_argument('--output', type=Path, default=None, help='Output file')
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()

    # File paths
    real_utr_fasta = project_root / "data" / "references" / "dmel_3utr.fasta"
    real_blast = results_dir / "genome_wide_all_3utrs.tsv"
    shuffled_dir = results_dir / "shuffled_full"
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_context_comparison.tsv"

    print("=" * 70)
    print("MOTIF CONTEXT COMPARISON ANALYSIS")
    print("=" * 70)

    # 1. Count k-mers in full 3'UTR sequences
    print(f"\n1. Counting {args.k}-mers in full 3'UTR sequences...")
    utr_counts, n_utrs, utr_bp = count_kmers_in_fasta(real_utr_fasta, args.k)
    print(f"   Sequences: {n_utrs:,}")
    print(f"   Total bp: {utr_bp:,}")
    print(f"   Total {args.k}-mers: {sum(utr_counts.values()):,}")

    # 2. Count k-mers in TE hits (real)
    print(f"\n2. Counting {args.k}-mers in real TE hits...")
    te_counts, n_hits, te_bp = count_kmers_in_blast_hits(real_blast, args.k)
    print(f"   TE hits: {n_hits:,}")
    print(f"   Total aligned bp: {te_bp:,}")
    print(f"   Total {args.k}-mers: {sum(te_counts.values()):,}")

    # 3. Count k-mers in shuffled UTRs and their TE hits
    print(f"\n3. Counting {args.k}-mers in shuffled sequences...")

    # Shuffled FASTA files are in the same directory as shuffled BLAST results
    shuffled_fasta_dir = shuffled_dir

    shuf_utr_counts = Counter()
    shuf_te_counts = Counter()
    n_shuf_utrs = 0
    shuf_utr_bp = 0
    n_shuf_hits = 0
    shuf_te_bp = 0
    n_reps = 0

    # Process shuffled replicates
    for rep in range(1, 11):
        # Shuffled FASTA (same naming as BLAST but .fasta extension)
        shuf_fasta = shuffled_fasta_dir / f"replicate_{rep:02d}.fasta"
        if not shuf_fasta.exists():
            # Try alternate naming
            shuf_fasta = shuffled_fasta_dir / f"replicate_{rep}.fasta"
        if shuf_fasta.exists():
            counts, n_seqs, bp = count_kmers_in_fasta(shuf_fasta, args.k)
            shuf_utr_counts.update(counts)
            n_shuf_utrs += n_seqs
            shuf_utr_bp += bp
            print(f"   Replicate {rep} FASTA: {n_seqs:,} seqs, {bp:,} bp")

        # Shuffled BLAST
        shuf_blast = shuffled_dir / f"replicate_{rep:02d}_blast.tsv"
        if shuf_blast.exists():
            counts, n_hits_rep, bp = count_kmers_in_blast_hits(shuf_blast, args.k)
            shuf_te_counts.update(counts)
            n_shuf_hits += n_hits_rep
            shuf_te_bp += bp
            n_reps += 1

    print(f"\n   Total shuffled UTRs: {n_shuf_utrs:,}")
    print(f"   Total shuffled UTR bp: {shuf_utr_bp:,}")
    print(f"   Total shuffled TE hits: {n_shuf_hits:,}")
    print(f"   Total shuffled TE bp: {shuf_te_bp:,}")

    # Generate all possible k-mers
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=args.k)]

    # Compile results
    print("\n4. Compiling results...")
    results = []

    for kmer in all_kmers:
        # Raw counts
        utr_count = utr_counts.get(kmer, 0)
        te_count = te_counts.get(kmer, 0)
        shuf_utr_count = shuf_utr_counts.get(kmer, 0)
        shuf_te_count = shuf_te_counts.get(kmer, 0)

        # Frequencies (per million k-mers)
        utr_total = sum(utr_counts.values())
        te_total = sum(te_counts.values())
        shuf_utr_total = sum(shuf_utr_counts.values())
        shuf_te_total = sum(shuf_te_counts.values())

        utr_freq = (utr_count / utr_total * 1e6) if utr_total > 0 else 0
        te_freq = (te_count / te_total * 1e6) if te_total > 0 else 0
        shuf_utr_freq = (shuf_utr_count / shuf_utr_total * 1e6) if shuf_utr_total > 0 else 0
        shuf_te_freq = (shuf_te_count / shuf_te_total * 1e6) if shuf_te_total > 0 else 0

        # Enrichment ratios
        # TE vs UTR (real): Are motifs enriched in TE hits vs background UTR?
        te_vs_utr = (te_freq / utr_freq) if utr_freq > 0 else float('inf')

        # TE vs shuffled TE: Are motifs enriched in real TE hits vs shuffled TE hits?
        te_vs_shuf_te = (te_freq / shuf_te_freq) if shuf_te_freq > 0 else float('inf')

        # Shuffled TE vs shuffled UTR: Compositional baseline
        shuf_te_vs_shuf_utr = (shuf_te_freq / shuf_utr_freq) if shuf_utr_freq > 0 else float('inf')

        results.append({
            'motif': kmer,
            # Raw counts
            'utr_count': utr_count,
            'te_hit_count': te_count,
            'shuf_utr_count': shuf_utr_count,
            'shuf_te_hit_count': shuf_te_count,
            # Frequencies (per million)
            'utr_freq_ppm': utr_freq,
            'te_hit_freq_ppm': te_freq,
            'shuf_utr_freq_ppm': shuf_utr_freq,
            'shuf_te_hit_freq_ppm': shuf_te_freq,
            # Enrichment ratios
            'te_vs_utr_ratio': te_vs_utr,
            'te_vs_shuf_te_ratio': te_vs_shuf_te,
            'shuf_te_vs_shuf_utr_ratio': shuf_te_vs_shuf_utr,
        })

    df = pd.DataFrame(results)

    # Save full results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.4g')
    print(f"\nFull results saved to: {output_path}")

    # Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    hdr_context = "Context"
    hdr_total = "Total k-mers"
    hdr_unique = "Unique motifs"
    label_utr = "Full 3'UTRs"
    label_te = "Real TE hits"
    label_shuf_utr = "Shuffled UTRs (10 reps)"
    label_shuf_te = "Shuffled TE hits (10 reps)"

    print(f"\n{hdr_context:<25} {hdr_total:>15} {hdr_unique:>15}")
    print("-" * 60)
    print(f"{label_utr:<25} {sum(utr_counts.values()):>15,} {len(utr_counts):>15,}")
    print(f"{label_te:<25} {sum(te_counts.values()):>15,} {len(te_counts):>15,}")
    print(f"{label_shuf_utr:<25} {sum(shuf_utr_counts.values()):>15,} {len(shuf_utr_counts):>15,}")
    print(f"{label_shuf_te:<25} {sum(shuf_te_counts.values()):>15,} {len(shuf_te_counts):>15,}")

    # Top motifs by different enrichment measures
    print("\n" + "=" * 70)
    print("TOP MOTIFS: Enriched in TE hits vs full UTR background")
    print("(te_vs_utr_ratio - motifs specifically in TEs, not just common in UTRs)")
    print("=" * 70)

    top_te_vs_utr = df[df['te_hit_count'] >= 1000].nlargest(20, 'te_vs_utr_ratio')
    print(f"\n{'Motif':<10} {'UTR count':>12} {'TE count':>12} {'TE/UTR ratio':>12}")
    print("-" * 50)
    for _, row in top_te_vs_utr.iterrows():
        print(f"{row['motif']:<10} {row['utr_count']:>12,} {row['te_hit_count']:>12,} {row['te_vs_utr_ratio']:>12.2f}x")

    print("\n" + "=" * 70)
    print("TOP MOTIFS: Enriched in real TE hits vs shuffled TE hits")
    print("(te_vs_shuf_te_ratio - the original analysis metric)")
    print("=" * 70)

    top_te_vs_shuf = df[df['te_hit_count'] >= 1000].nlargest(20, 'te_vs_shuf_te_ratio')
    print(f"\n{'Motif':<10} {'Real TE':>12} {'Shuf TE':>12} {'Real/Shuf':>12}")
    print("-" * 50)
    for _, row in top_te_vs_shuf.iterrows():
        shuf_per_rep = row['shuf_te_hit_count'] / 10 if n_reps == 10 else row['shuf_te_hit_count']
        print(f"{row['motif']:<10} {row['te_hit_count']:>12,} {shuf_per_rep:>12,.0f} {row['te_vs_shuf_te_ratio']:>12.2f}x")

    print("\n" + "=" * 70)
    print("KNOWN REGULATORY MOTIFS - FULL CONTEXT")
    print("=" * 70)

    known_motifs = ['AATAAA', 'ATTAAA', 'TGTATA', 'TGTAAA', 'TATTTA', 'CAGCAG', 'ACACAC', 'GGTAAG']

    print(f"\n{'Motif':<10} {'UTR':>10} {'TE hit':>10} {'Shuf UTR':>12} {'Shuf TE':>10} {'TE/UTR':>8} {'TE/ShufTE':>10}")
    print("-" * 80)

    for motif in known_motifs:
        row = df[df['motif'] == motif].iloc[0]
        shuf_utr_per_rep = row['shuf_utr_count'] / 10
        shuf_te_per_rep = row['shuf_te_hit_count'] / 10
        print(f"{row['motif']:<10} {row['utr_count']:>10,} {row['te_hit_count']:>10,} "
              f"{shuf_utr_per_rep:>12,.0f} {shuf_te_per_rep:>10,.0f} "
              f"{row['te_vs_utr_ratio']:>8.2f}x {row['te_vs_shuf_te_ratio']:>9.2f}x")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    print("""
Two key ratios tell different stories:

1. TE/UTR ratio: Is this motif MORE common in TE-aligned regions than in
   3'UTRs overall? High ratio = motif is specifically associated with TEs,
   not just common background.

2. TE/ShufTE ratio: Is this motif MORE common in real TE hits than in
   shuffled TE hits? High ratio = cannot be explained by composition alone.

Motifs with BOTH high ratios are the strongest candidates for functional
TE-derived regulatory elements.
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
