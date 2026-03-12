#!/usr/bin/env python3
"""
For each motif occurrence in UTR sequences, determine if it overlaps
with a TE hit region. Compare real UTRs vs shuffled UTRs.

This answers: "Of all AATAAA motifs in 3'UTRs, what fraction fall within
TE-aligned regions? And how does this compare to shuffled controls?"

Approach:
1. Real analysis: Find motifs in real UTRs, check overlap with real TE hits
2. Shuffled analysis: Find motifs in shuffled UTRs, check overlap with shuffled TE hits
3. Compare the fractions
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter
from itertools import product

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def find_all_motif_positions(sequence: str, motif: str) -> list:
    """Find all start positions of motif in sequence (0-indexed)."""
    sequence = sequence.upper()
    motif = motif.upper()
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def build_te_coverage(blast_path: Path) -> dict:
    """
    Build dictionary mapping transcript_id -> list of merged (start, end) intervals.
    """
    hits_by_transcript = defaultdict(list)

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        qstart = hit.get('qstart', 0)
        qend = hit.get('qend', 0)

        if qseqid and qstart and qend:
            start = min(qstart, qend) - 1  # Convert to 0-based
            end = max(qstart, qend)
            hits_by_transcript[qseqid].append((start, end))

    # Merge overlapping intervals
    for transcript_id in hits_by_transcript:
        intervals = sorted(hits_by_transcript[transcript_id])
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        hits_by_transcript[transcript_id] = merged

    return dict(hits_by_transcript)


def position_overlaps_intervals(pos: int, length: int, intervals: list) -> bool:
    """Check if position [pos, pos+length) overlaps any interval."""
    motif_end = pos + length
    for start, end in intervals:
        if pos < end and motif_end > start:
            return True
    return False


def analyze_motif_te_overlap(sequences: dict, te_coverage: dict, k: int) -> tuple:
    """
    For each k-mer, count total occurrences and how many fall in TE regions.

    Returns: (motif_total, motif_in_te, total_bp, te_bp)
    """
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=k)]

    motif_total = Counter()
    motif_in_te = Counter()

    total_bp = 0
    te_bp = 0

    for seq_id, sequence in sequences.items():
        total_bp += len(sequence)

        # Calculate TE coverage for this sequence
        intervals = te_coverage.get(seq_id, [])
        for start, end in intervals:
            te_bp += (end - start)

        # Find all motif positions and check TE overlap
        for motif in all_kmers:
            positions = find_all_motif_positions(sequence, motif)
            for pos in positions:
                motif_total[motif] += 1
                if position_overlaps_intervals(pos, k, intervals):
                    motif_in_te[motif] += 1

    return motif_total, motif_in_te, total_bp, te_bp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k', type=int, default=6, help='K-mer size')
    parser.add_argument('--n-reps', type=int, default=3, help='Number of shuffled reps to analyze')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    shuffled_dir = results_dir / "shuffled_full"
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_te_overlap_v2.tsv"

    print("=" * 70)
    print("MOTIF TE OVERLAP ANALYSIS (v2)")
    print("=" * 70)

    # === REAL DATA ===
    print("\n" + "=" * 70)
    print("ANALYZING REAL DATA")
    print("=" * 70)

    # Load real UTR sequences
    print("\n1. Loading real UTR sequences...")
    utr_fasta = project_root / "data" / "references" / "dmel_3utr.fasta"
    real_utrs = {}
    for record in SeqIO.parse(utr_fasta, "fasta"):
        real_utrs[record.id] = str(record.seq).upper()
    print(f"   Loaded {len(real_utrs):,} sequences")

    # Build real TE coverage
    print("\n2. Building TE coverage from real BLAST hits...")
    real_blast = results_dir / "genome_wide_all_3utrs.tsv"
    real_te_cov = build_te_coverage(real_blast)
    print(f"   UTRs with TE hits: {len(real_te_cov):,}")

    # Analyze real motifs
    print("\n3. Analyzing motif-TE overlaps in real data...")
    real_total, real_in_te, real_bp, real_te_bp = analyze_motif_te_overlap(
        real_utrs, real_te_cov, args.k
    )
    real_coverage = 100 * real_te_bp / real_bp
    print(f"   Total bp: {real_bp:,}")
    print(f"   TE-covered bp: {real_te_bp:,}")
    print(f"   TE coverage: {real_coverage:.1f}%")

    # === SHUFFLED DATA ===
    print("\n" + "=" * 70)
    print(f"ANALYZING SHUFFLED DATA ({args.n_reps} replicates)")
    print("=" * 70)

    shuf_totals = []
    shuf_in_tes = []
    shuf_coverages = []

    for rep in range(1, args.n_reps + 1):
        print(f"\n--- Replicate {rep} ---")

        # Load shuffled sequences
        shuf_fasta = shuffled_dir / f"replicate_{rep:02d}.fasta"
        print(f"  Loading {shuf_fasta.name}...")
        shuf_seqs = {}
        for record in SeqIO.parse(shuf_fasta, "fasta"):
            shuf_seqs[record.id] = str(record.seq).upper()

        # Load shuffled TE coverage
        shuf_blast = shuffled_dir / f"replicate_{rep:02d}_blast.tsv"
        print(f"  Building TE coverage from {shuf_blast.name}...")
        shuf_te_cov = build_te_coverage(shuf_blast)

        # Analyze
        print(f"  Analyzing motif-TE overlaps...")
        s_total, s_in_te, s_bp, s_te_bp = analyze_motif_te_overlap(
            shuf_seqs, shuf_te_cov, args.k
        )

        coverage = 100 * s_te_bp / s_bp
        print(f"  TE coverage: {coverage:.1f}%")

        shuf_totals.append(s_total)
        shuf_in_tes.append(s_in_te)
        shuf_coverages.append(coverage)

    mean_shuf_coverage = sum(shuf_coverages) / len(shuf_coverages)

    # === COMPILE RESULTS ===
    print("\n" + "=" * 70)
    print("COMPILING RESULTS")
    print("=" * 70)

    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=args.k)]

    results = []
    for motif in all_kmers:
        # Real data
        r_total = real_total[motif]
        r_in_te = real_in_te[motif]
        r_frac = r_in_te / r_total if r_total > 0 else 0

        # Shuffled data (average across replicates)
        s_total = sum(st[motif] for st in shuf_totals) / len(shuf_totals)
        s_in_te = sum(sit[motif] for sit in shuf_in_tes) / len(shuf_in_tes)
        s_frac = s_in_te / s_total if s_total > 0 else 0

        # Enrichment vs coverage baseline
        enrich_vs_coverage = r_frac / (real_coverage / 100) if real_coverage > 0 else 0

        # Enrichment real vs shuffled
        enrich_vs_shuf = r_frac / s_frac if s_frac > 0 else float('inf')

        results.append({
            'motif': motif,
            'real_total': r_total,
            'real_in_te': r_in_te,
            'real_not_in_te': r_total - r_in_te,
            'real_pct_in_te': r_frac * 100,
            'shuf_total': s_total,
            'shuf_in_te': s_in_te,
            'shuf_pct_in_te': s_frac * 100,
            'enrich_vs_coverage': enrich_vs_coverage,
            'enrich_real_vs_shuf': enrich_vs_shuf,
        })

    df = pd.DataFrame(results)
    df.to_csv(output_path, sep='\t', index=False, float_format='%.4g')
    print(f"\nResults saved to: {output_path}")

    # === SUMMARY ===
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\nTE Coverage:")
    print(f"  Real UTRs: {real_coverage:.1f}%")
    print(f"  Shuffled UTRs: {mean_shuf_coverage:.1f}% (mean of {args.n_reps} reps)")

    print(f"\nExpected % in TE if random: {real_coverage:.1f}%")

    # Top motifs enriched in TE
    print("\n" + "-" * 70)
    print("TOP 20: Highest fraction in real TE hits (min 1000 occurrences)")
    print("-" * 70)

    hdr = f"{'Motif':<8} {'Total':>8} {'In TE':>8} {'%inTE':>7} {'Shuf%':>7} {'vsBase':>7} {'vsShuf':>7}"
    print(hdr)
    print("-" * len(hdr))

    top = df[df['real_total'] >= 1000].nlargest(20, 'real_pct_in_te')
    for _, row in top.iterrows():
        vs_shuf = f"{row['enrich_real_vs_shuf']:.2f}x" if row['enrich_real_vs_shuf'] < 100 else ">100x"
        print(f"{row['motif']:<8} {row['real_total']:>8,} {row['real_in_te']:>8,} "
              f"{row['real_pct_in_te']:>6.1f}% {row['shuf_pct_in_te']:>6.1f}% "
              f"{row['enrich_vs_coverage']:>6.2f}x {vs_shuf:>7}")

    # Bottom motifs (depleted in TE)
    print("\n" + "-" * 70)
    print("BOTTOM 20: Lowest fraction in real TE hits (min 1000 occurrences)")
    print("-" * 70)

    print(hdr)
    print("-" * len(hdr))

    bottom = df[df['real_total'] >= 1000].nsmallest(20, 'real_pct_in_te')
    for _, row in bottom.iterrows():
        vs_shuf = f"{row['enrich_real_vs_shuf']:.2f}x" if row['enrich_real_vs_shuf'] < 100 else ">100x"
        print(f"{row['motif']:<8} {row['real_total']:>8,} {row['real_in_te']:>8,} "
              f"{row['real_pct_in_te']:>6.1f}% {row['shuf_pct_in_te']:>6.1f}% "
              f"{row['enrich_vs_coverage']:>6.2f}x {vs_shuf:>7}")

    # Known regulatory motifs
    print("\n" + "=" * 70)
    print("KNOWN REGULATORY MOTIFS")
    print("=" * 70)

    known = ['AATAAA', 'ATTAAA', 'TGTATA', 'TGTAAA', 'TATTTA', 'CAGCAG', 'ACACAC', 'GGTAAG']

    print(f"\n{'Motif':<8} {'Total':>8} {'In TE':>8} {'Not TE':>8} {'%inTE':>7} {'Shuf%':>7} {'vsShuf':>7}")
    print("-" * 65)

    for motif in known:
        row = df[df['motif'] == motif].iloc[0]
        vs_shuf = f"{row['enrich_real_vs_shuf']:.2f}x" if row['enrich_real_vs_shuf'] < 100 else ">100x"
        print(f"{row['motif']:<8} {row['real_total']:>8,} {row['real_in_te']:>8,} "
              f"{row['real_not_in_te']:>8,} {row['real_pct_in_te']:>6.1f}% "
              f"{row['shuf_pct_in_te']:>6.1f}% {vs_shuf:>7}")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    te_cov_str = f"{real_coverage:.1f}"
    shuf_cov_str = f"{mean_shuf_coverage:.1f}"

    print(f"""
TE coverage: {te_cov_str}% of real UTR bp, {shuf_cov_str}% of shuffled UTR bp

For any motif:
  - If %inTE = {te_cov_str}% => motif is randomly distributed (neutral)
  - If %inTE > {te_cov_str}% => motif is ENRICHED in TE regions
  - If %inTE < {te_cov_str}% => motif is DEPLETED in TE regions

'vsShuf' compares real vs shuffled:
  - >1x = more of this motif falls in TE regions in real data than shuffled
  - <1x = less in TE regions compared to shuffled
  - =1x = same as shuffled (purely compositional)

Motifs with HIGH %inTE AND HIGH vsShuf = genuinely TE-associated
Motifs with HIGH %inTE but vsShuf~1 = common in UTRs that happen to have TEs
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
