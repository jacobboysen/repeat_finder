#!/usr/bin/env python3
"""
For each motif occurrence in a 3'UTR sequence, determine if it overlaps
with a TE hit region or not.

This answers: "Of all AATAAA motifs in 3'UTRs, what fraction fall within
TE-aligned regions?"

Approach:
1. Load UTR sequences and find all motif positions
2. Load deduplicated TE hits and build interval coverage per UTR
3. For each motif occurrence, check if it overlaps a TE hit
4. Compute: motifs_in_te / total_motifs for each motif type
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


def build_te_coverage(blast_path: Path, deduplicate: bool = True) -> dict:
    """
    Build a dictionary mapping transcript_id -> list of (start, end) intervals.

    If deduplicate=True, collapse overlapping intervals from same gene.
    """
    # First pass: collect all hits per transcript
    hits_by_transcript = defaultdict(list)

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        qstart = hit.get('qstart', 0)
        qend = hit.get('qend', 0)

        if qseqid and qstart and qend:
            # BLAST coordinates are 1-based, convert to 0-based
            start = min(qstart, qend) - 1
            end = max(qstart, qend)
            hits_by_transcript[qseqid].append((start, end))

    # Merge overlapping intervals if deduplicating
    if deduplicate:
        for transcript_id in hits_by_transcript:
            intervals = sorted(hits_by_transcript[transcript_id])
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1]:
                    # Overlaps with previous, extend it
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            hits_by_transcript[transcript_id] = merged

    return dict(hits_by_transcript)


def position_overlaps_intervals(pos: int, length: int, intervals: list) -> bool:
    """Check if position [pos, pos+length) overlaps any interval."""
    motif_end = pos + length
    for start, end in intervals:
        # Overlap if motif_start < interval_end AND motif_end > interval_start
        if pos < end and motif_end > start:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k', type=int, default=6, help='K-mer size (default: 6)')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()

    # File paths
    utr_fasta = project_root / "data" / "references" / "dmel_3utr.fasta"
    real_blast = results_dir / "genome_wide_all_3utrs.tsv"
    shuffled_dir = results_dir / "shuffled_full"
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_te_overlap_analysis.tsv"

    print("=" * 70)
    print("MOTIF TE OVERLAP ANALYSIS")
    print("=" * 70)

    # Load UTR sequences
    print("\n1. Loading UTR sequences...")
    utrs = {}
    for record in SeqIO.parse(utr_fasta, "fasta"):
        utrs[record.id] = str(record.seq).upper()
    print(f"   Loaded {len(utrs):,} UTR sequences")

    # Build TE coverage from real hits (deduplicated)
    print("\n2. Building TE coverage from real hits (deduplicated)...")
    real_te_coverage = build_te_coverage(real_blast, deduplicate=True)
    total_te_bp = sum(sum(end - start for start, end in intervals)
                      for intervals in real_te_coverage.values())
    print(f"   UTRs with TE hits: {len(real_te_coverage):,}")
    print(f"   Total TE-covered bp: {total_te_bp:,}")

    # Build TE coverage from shuffled hits (use first replicate as representative)
    print("\n3. Building TE coverage from shuffled hits (rep 1)...")
    shuf_blast = shuffled_dir / "replicate_01_blast.tsv"
    shuf_te_coverage = build_te_coverage(shuf_blast, deduplicate=True)
    shuf_te_bp = sum(sum(end - start for start, end in intervals)
                     for intervals in shuf_te_coverage.values())
    print(f"   UTRs with shuffled TE hits: {len(shuf_te_coverage):,}")
    print(f"   Total shuffled TE-covered bp: {shuf_te_bp:,}")

    # Generate all k-mers to analyze
    bases = ['A', 'C', 'G', 'T']
    all_kmers = [''.join(p) for p in product(bases, repeat=args.k)]

    # Count motif occurrences and TE overlaps
    print(f"\n4. Analyzing {len(all_kmers)} motifs across all UTRs...")

    # Initialize counters
    motif_total = Counter()           # Total occurrences in UTRs
    motif_in_real_te = Counter()      # Occurrences overlapping real TE hits
    motif_in_shuf_te = Counter()      # Occurrences overlapping shuffled TE hits

    total_utr_bp = 0
    n_processed = 0

    for transcript_id, sequence in utrs.items():
        total_utr_bp += len(sequence)
        n_processed += 1

        if n_processed % 5000 == 0:
            print(f"   Processed {n_processed:,} / {len(utrs):,} UTRs...")

        # Get TE intervals for this transcript
        real_intervals = real_te_coverage.get(transcript_id, [])
        shuf_intervals = shuf_te_coverage.get(transcript_id, [])

        # Find all motif occurrences
        for motif in all_kmers:
            positions = find_all_motif_positions(sequence, motif)

            for pos in positions:
                motif_total[motif] += 1

                if position_overlaps_intervals(pos, args.k, real_intervals):
                    motif_in_real_te[motif] += 1

                if position_overlaps_intervals(pos, args.k, shuf_intervals):
                    motif_in_shuf_te[motif] += 1

    print(f"   Total UTR bp: {total_utr_bp:,}")

    # Compile results
    print("\n5. Compiling results...")

    results = []
    for motif in all_kmers:
        total = motif_total[motif]
        in_real = motif_in_real_te[motif]
        in_shuf = motif_in_shuf_te[motif]
        not_in_real = total - in_real
        not_in_shuf = total - in_shuf

        # Fractions
        frac_in_real = in_real / total if total > 0 else 0
        frac_in_shuf = in_shuf / total if total > 0 else 0

        # Enrichment: observed / expected based on coverage
        # Expected = total * (te_bp / utr_bp)
        expected_real = total * (total_te_bp / total_utr_bp) if total_utr_bp > 0 else 0
        expected_shuf = total * (shuf_te_bp / total_utr_bp) if total_utr_bp > 0 else 0

        enrich_real = in_real / expected_real if expected_real > 0 else float('inf')
        enrich_shuf = in_shuf / expected_shuf if expected_shuf > 0 else float('inf')

        # Enrichment of real vs shuffled
        enrich_real_vs_shuf = frac_in_real / frac_in_shuf if frac_in_shuf > 0 else float('inf')

        results.append({
            'motif': motif,
            'total_in_utrs': total,
            'in_real_te': in_real,
            'not_in_real_te': not_in_real,
            'frac_in_real_te': frac_in_real,
            'in_shuf_te': in_shuf,
            'not_in_shuf_te': not_in_shuf,
            'frac_in_shuf_te': frac_in_shuf,
            'enrich_vs_coverage': enrich_real,
            'enrich_real_vs_shuf': enrich_real_vs_shuf,
        })

    df = pd.DataFrame(results)
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY: What fraction of each motif falls within TE hits?")
    print("=" * 70)

    # Overall coverage
    coverage_pct = 100 * total_te_bp / total_utr_bp
    shuf_coverage_pct = 100 * shuf_te_bp / total_utr_bp
    print(f"\nTE coverage of UTRs: {coverage_pct:.1f}% (real), {shuf_coverage_pct:.1f}% (shuffled)")
    print(f"Expected fraction in TE (by coverage): {coverage_pct:.1f}%")

    # Top motifs by fraction in TE
    print("\n" + "-" * 70)
    print("TOP 20 MOTIFS: Highest fraction found in TE-hit regions")
    print("-" * 70)

    top_frac = df[df['total_in_utrs'] >= 1000].nlargest(20, 'frac_in_real_te')
    print(f"\n{'Motif':<10} {'Total':>10} {'In TE':>10} {'% in TE':>10} {'% in Shuf':>10} {'Enrich':>10}")
    print("-" * 65)
    for _, row in top_frac.iterrows():
        print(f"{row['motif']:<10} {row['total_in_utrs']:>10,} {row['in_real_te']:>10,} "
              f"{row['frac_in_real_te']*100:>9.1f}% {row['frac_in_shuf_te']*100:>9.1f}% "
              f"{row['enrich_real_vs_shuf']:>9.2f}x")

    # Bottom motifs (lowest fraction in TE)
    print("\n" + "-" * 70)
    print("BOTTOM 20 MOTIFS: Lowest fraction found in TE-hit regions")
    print("-" * 70)

    bottom_frac = df[df['total_in_utrs'] >= 1000].nsmallest(20, 'frac_in_real_te')
    print(f"\n{'Motif':<10} {'Total':>10} {'In TE':>10} {'% in TE':>10} {'% in Shuf':>10} {'Enrich':>10}")
    print("-" * 65)
    for _, row in bottom_frac.iterrows():
        print(f"{row['motif']:<10} {row['total_in_utrs']:>10,} {row['in_real_te']:>10,} "
              f"{row['frac_in_real_te']*100:>9.1f}% {row['frac_in_shuf_te']*100:>9.1f}% "
              f"{row['enrich_real_vs_shuf']:>9.2f}x")

    # Known regulatory motifs
    print("\n" + "=" * 70)
    print("KNOWN REGULATORY MOTIFS")
    print("=" * 70)

    known = ['AATAAA', 'ATTAAA', 'TGTATA', 'TGTAAA', 'TATTTA', 'CAGCAG', 'ACACAC', 'GGTAAG']

    print(f"\n{'Motif':<10} {'Total':>10} {'In TE':>10} {'Not in TE':>10} {'% in TE':>10} {'vs Shuf':>10}")
    print("-" * 65)
    for motif in known:
        row = df[df['motif'] == motif].iloc[0]
        print(f"{row['motif']:<10} {row['total_in_utrs']:>10,} {row['in_real_te']:>10,} "
              f"{row['not_in_real_te']:>10,} {row['frac_in_real_te']*100:>9.1f}% "
              f"{row['enrich_real_vs_shuf']:>9.2f}x")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"""
TE coverage = {coverage_pct:.1f}% of UTR sequence

If motifs were randomly distributed:
  - Expected ~{coverage_pct:.1f}% of any motif would fall in TE regions

Motifs with HIGHER % in TE = enriched in TE-hit regions
Motifs with LOWER % in TE = depleted/avoided in TE-hit regions

'Enrich' = (% in real TE) / (% in shuffled TE)
  - >1 = more likely in real TE hits than shuffled
  - <1 = less likely in real TE hits than shuffled
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
