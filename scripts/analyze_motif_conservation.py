#!/usr/bin/env python3
"""
Compare phyloP conservation scores for motifs INSIDE vs OUTSIDE TE-hit regions.

This answers: Are regulatory motifs that fall within TE-aligned regions more
conserved than the same motifs outside TE regions?

Workflow:
1. Find all positions of target motifs in UTRs
2. Convert to genomic coordinates
3. Query phyloP bigWig for conservation scores
4. Classify each motif as in-TE or not-in-TE
5. Compare conservation distributions
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from collections import defaultdict
import re

import pandas as pd
import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def parse_gff_utrs(gff_path: Path) -> dict:
    """
    Parse GFF to get 3'UTR genomic coordinates.
    Returns: {transcript_id: (chrom, start, end, strand)}
    """
    utrs = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature = parts[2]
            if feature != 'three_prime_UTR':
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Parse attributes for parent transcript
            attrs = parts[8]
            parent_match = re.search(r'Parent=([^;,]+)', attrs)
            if parent_match:
                transcript_id = parent_match.group(1)
                # Handle multiple UTR segments per transcript - take union
                if transcript_id in utrs:
                    old = utrs[transcript_id]
                    utrs[transcript_id] = (
                        chrom,
                        min(old[1], start),
                        max(old[2], end),
                        strand
                    )
                else:
                    utrs[transcript_id] = (chrom, start, end, strand)

    return utrs


def find_motif_positions(sequence: str, motif: str) -> list:
    """Find all start positions of motif in sequence."""
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


def utr_to_genomic(utr_pos: int, motif_len: int, utr_info: tuple) -> tuple:
    """
    Convert UTR-relative position to genomic coordinates.

    Args:
        utr_pos: 0-based position in UTR sequence
        motif_len: length of motif
        utr_info: (chrom, start, end, strand) from GFF (1-based)

    Returns: (chrom, genomic_start, genomic_end) in 0-based BED format
    """
    chrom, utr_start, utr_end, strand = utr_info

    if strand == '+':
        # Forward strand: UTR sequence goes 5'->3' from utr_start to utr_end
        genomic_start = (utr_start - 1) + utr_pos  # Convert to 0-based
        genomic_end = genomic_start + motif_len
    else:
        # Reverse strand: UTR sequence is reverse complement
        # Position 0 in UTR corresponds to utr_end in genome
        genomic_end = utr_end - utr_pos
        genomic_start = genomic_end - motif_len

    return (chrom, genomic_start, genomic_end)


def build_te_intervals(blast_path: Path, utr_coords: dict) -> dict:
    """
    Build merged TE intervals in genomic coordinates.
    Returns: {chrom: [(start, end), ...]}
    """
    intervals_by_chrom = defaultdict(list)

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        qstart = hit.get('qstart', 0)
        qend = hit.get('qend', 0)

        if qseqid not in utr_coords:
            continue

        utr_info = utr_coords[qseqid]
        chrom, utr_g_start, utr_g_end, strand = utr_info

        # Convert hit coordinates to genomic
        hit_utr_start = min(qstart, qend) - 1  # 0-based
        hit_utr_end = max(qstart, qend)
        hit_len = hit_utr_end - hit_utr_start

        if strand == '+':
            g_start = (utr_g_start - 1) + hit_utr_start
            g_end = g_start + hit_len
        else:
            g_end = utr_g_end - hit_utr_start
            g_start = g_end - hit_len

        # Add chr prefix if not present
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom

        intervals_by_chrom[chrom].append((g_start, g_end))

    # Merge overlapping intervals
    merged = {}
    for chrom, intervals in intervals_by_chrom.items():
        intervals.sort()
        result = []
        for start, end in intervals:
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            else:
                result.append((start, end))
        merged[chrom] = result

    return merged


def position_in_intervals(chrom: str, start: int, end: int, intervals: dict) -> bool:
    """Check if position overlaps any interval."""
    if chrom not in intervals:
        return False
    for iv_start, iv_end in intervals[chrom]:
        if start < iv_end and end > iv_start:
            return True
    return False


def run_bigwig_average(bed_path: Path, bigwig_path: Path, output_path: Path):
    """Run bigWigAverageOverBed to get conservation scores."""
    cmd = ['bigWigAverageOverBed', str(bigwig_path), str(bed_path), str(output_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running bigWigAverageOverBed: {result.stderr}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--motifs', type=str, default='AATAAA,ATTAAA,TGTATA,CAGCAG,ACACAC,GGTAAG',
                        help='Comma-separated list of motifs to analyze')
    parser.add_argument('--sample-rate', type=float, default=0.1,
                        help='Fraction of motifs to sample (default: 0.1 = 10%%)')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_conservation_comparison.tsv"

    # File paths
    utr_fasta = project_root / "data" / "references" / "dmel_3utr.fasta"
    gff_path = project_root / "data" / "references" / "dmel-all-r6.66.gff"
    blast_path = results_dir / "genome_wide_all_3utrs.tsv"
    phylop_bw = project_root / "data" / "references" / "dm6.phyloP27way.bw"

    motifs = [m.strip().upper() for m in args.motifs.split(',')]

    print("=" * 70)
    print("MOTIF CONSERVATION ANALYSIS: TE vs Non-TE")
    print("=" * 70)
    print(f"Motifs: {', '.join(motifs)}")
    print(f"Sample rate: {args.sample_rate * 100:.0f}%")

    # Check for bigWigAverageOverBed
    try:
        subprocess.run(['bigWigAverageOverBed'], capture_output=True)
    except FileNotFoundError:
        print("\nError: bigWigAverageOverBed not found in PATH")
        print("Install with: conda install -c bioconda ucsc-bigwigaverageoverbed")
        return 1

    # Load UTR coordinates from GFF
    print("\n1. Loading UTR coordinates from GFF...")
    utr_coords = parse_gff_utrs(gff_path)
    print(f"   Found {len(utr_coords):,} UTRs with coordinates")

    # Load UTR sequences
    print("\n2. Loading UTR sequences...")
    utr_seqs = {}
    for record in SeqIO.parse(utr_fasta, "fasta"):
        utr_seqs[record.id] = str(record.seq).upper()
    print(f"   Loaded {len(utr_seqs):,} sequences")

    # Build TE intervals
    print("\n3. Building TE hit intervals...")
    te_intervals = build_te_intervals(blast_path, utr_coords)
    total_intervals = sum(len(v) for v in te_intervals.values())
    print(f"   Built {total_intervals:,} merged intervals across {len(te_intervals)} chromosomes")

    # Find all motif positions and classify
    print("\n4. Finding motif positions and classifying...")

    motif_records = []
    np.random.seed(42)

    for motif in motifs:
        k = len(motif)
        n_found = 0
        n_sampled = 0

        for transcript_id, sequence in utr_seqs.items():
            if transcript_id not in utr_coords:
                continue

            utr_info = utr_coords[transcript_id]
            positions = find_motif_positions(sequence, motif)

            for pos in positions:
                n_found += 1

                # Sample
                if np.random.random() > args.sample_rate:
                    continue
                n_sampled += 1

                # Convert to genomic
                chrom, g_start, g_end = utr_to_genomic(pos, k, utr_info)

                # Add chr prefix
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom

                # Check TE overlap
                in_te = position_in_intervals(chrom, g_start, g_end, te_intervals)

                motif_records.append({
                    'motif': motif,
                    'transcript': transcript_id,
                    'utr_pos': pos,
                    'chrom': chrom,
                    'start': g_start,
                    'end': g_end,
                    'in_te': in_te,
                })

        print(f"   {motif}: {n_found:,} found, {n_sampled:,} sampled")

    print(f"\n   Total sampled motif positions: {len(motif_records):,}")

    if len(motif_records) == 0:
        print("No motif positions found!")
        return 1

    # Create BED file for conservation query
    print("\n5. Querying phyloP conservation scores...")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_f:
        bed_path = Path(bed_f.name)
        for i, rec in enumerate(motif_records):
            # BED format: chrom, start, end, name
            bed_f.write(f"{rec['chrom']}\t{rec['start']}\t{rec['end']}\t{i}\n")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.tab', delete=False) as out_f:
        tab_path = Path(out_f.name)

    if not run_bigwig_average(bed_path, phylop_bw, tab_path):
        return 1

    # Parse conservation results
    print("   Parsing conservation scores...")
    conservation = {}
    with open(tab_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                idx = int(parts[0])
                mean_score = float(parts[5])  # mean column
                conservation[idx] = mean_score

    # Add conservation to records
    for i, rec in enumerate(motif_records):
        rec['phylop'] = conservation.get(i, np.nan)

    # Clean up temp files
    bed_path.unlink()
    tab_path.unlink()

    # Create DataFrame and analyze
    print("\n6. Analyzing results...")
    df = pd.DataFrame(motif_records)

    # Remove records with missing conservation
    df = df.dropna(subset=['phylop'])
    print(f"   Records with conservation scores: {len(df):,}")

    # Summary by motif and TE status
    print("\n" + "=" * 70)
    print("RESULTS: Conservation by Motif and TE Status")
    print("=" * 70)

    summary_rows = []

    for motif in motifs:
        motif_df = df[df['motif'] == motif]

        in_te = motif_df[motif_df['in_te'] == True]
        not_in_te = motif_df[motif_df['in_te'] == False]

        if len(in_te) > 0 and len(not_in_te) > 0:
            in_te_mean = in_te['phylop'].mean()
            not_in_te_mean = not_in_te['phylop'].mean()

            in_te_conserved = (in_te['phylop'] > 1).sum() / len(in_te) * 100
            not_in_te_conserved = (not_in_te['phylop'] > 1).sum() / len(not_in_te) * 100

            summary_rows.append({
                'motif': motif,
                'n_in_te': len(in_te),
                'n_not_in_te': len(not_in_te),
                'mean_phylop_in_te': in_te_mean,
                'mean_phylop_not_in_te': not_in_te_mean,
                'pct_conserved_in_te': in_te_conserved,
                'pct_conserved_not_in_te': not_in_te_conserved,
                'phylop_diff': in_te_mean - not_in_te_mean,
            })

    summary_df = pd.DataFrame(summary_rows)

    # Print results
    print(f"\n{'Motif':<8} {'N in TE':>10} {'N not TE':>10} {'Mean(TE)':>10} {'Mean(noTE)':>10} {'Diff':>8} {'%Cons(TE)':>10} {'%Cons(no)':>10}")
    print("-" * 90)

    for _, row in summary_df.iterrows():
        diff_str = f"{row['phylop_diff']:+.3f}"
        print(f"{row['motif']:<8} {row['n_in_te']:>10,} {row['n_not_in_te']:>10,} "
              f"{row['mean_phylop_in_te']:>10.3f} {row['mean_phylop_not_in_te']:>10.3f} "
              f"{diff_str:>8} {row['pct_conserved_in_te']:>9.1f}% {row['pct_conserved_not_in_te']:>9.1f}%")

    # Save detailed results
    df.to_csv(output_path, sep='\t', index=False)
    summary_path = output_path.with_name('motif_conservation_summary.tsv')
    summary_df.to_csv(summary_path, sep='\t', index=False)

    print(f"\nDetailed results: {output_path}")
    print(f"Summary: {summary_path}")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
phyloP > 0 = conserved (purifying selection)
phyloP < 0 = fast-evolving (positive selection or relaxed constraint)
phyloP > 1 = typically considered "conserved"
phyloP > 2 = highly conserved

'Diff' = Mean(in TE) - Mean(not in TE)
  Positive = motifs in TE regions are MORE conserved
  Negative = motifs in TE regions are LESS conserved
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
