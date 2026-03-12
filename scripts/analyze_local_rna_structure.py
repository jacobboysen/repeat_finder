#!/usr/bin/env python3
"""
Local RNA Structure Analysis using RNAplfold

Computes position-specific accessibility (probability of being unpaired) across
full UTR sequences, then compares TE-hit positions vs non-TE positions.

Uses literature-recommended parameters from Lange et al. 2012 (PMC3384308):
- L = 150 (max base pair span)
- W = 200 (window size)

Unlike global folding (RNAfold), RNAplfold:
- Uses sliding window - no manual tiling needed
- Handles long sequences efficiently O(n*L*L)
- Provides position-specific accessibility scores
- Better for UTR analysis (local structures dominate)

Author: Claude Code
Date: 2026-01-18
"""

import argparse
import subprocess
import tempfile
import os
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.utils.blast_io import iter_blast_results, BLAST_COLUMNS
from scripts.utils.paths import get_project_root


def parse_fasta(fasta_path):
    """Parse FASTA file into dict of {seqid: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract just the transcript ID
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def load_te_intervals(blast_path, max_hits=None):
    """
    Load TE hit intervals per transcript.
    Returns dict: {transcript_id: [(start, end), ...]}
    Coordinates are 0-based, half-open.
    """
    intervals = defaultdict(list)
    count = 0

    for hit in iter_blast_results(blast_path):
        qseqid = hit['qseqid']
        # Convert to 0-based
        start = min(hit['qstart'], hit['qend']) - 1
        end = max(hit['qstart'], hit['qend'])
        intervals[qseqid].append((start, end))

        count += 1
        if max_hits and count >= max_hits:
            break

    return intervals


def merge_intervals(intervals):
    """Merge overlapping intervals."""
    if not intervals:
        return []

    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]

    for start, end in sorted_intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    return merged


def run_rnaplfold(sequence, seq_id, span=150, winsize=200, ulength=31,
                  rnaplfold_path='RNAplfold', workdir=None):
    """
    Run RNAplfold on a sequence and return accessibility array.

    Returns numpy array of shape (seq_len,) with accessibility (unpaired probability)
    for each position.
    """
    if workdir is None:
        workdir = tempfile.mkdtemp()

    # Write sequence to temp file
    fasta_path = os.path.join(workdir, 'input.fa')
    with open(fasta_path, 'w') as f:
        f.write(f">{seq_id}\n{sequence}\n")

    # Run RNAplfold
    cmd = [
        rnaplfold_path,
        '-L', str(span),
        '-W', str(winsize),
        '-u', str(ulength),
    ]

    try:
        result = subprocess.run(
            cmd,
            stdin=open(fasta_path),
            capture_output=True,
            text=True,
            cwd=workdir,
            timeout=300  # 5 min timeout per sequence
        )
    except subprocess.TimeoutExpired:
        return None

    if result.returncode != 0:
        return None

    # Parse _lunp output file
    lunp_file = os.path.join(workdir, f"{seq_id}_lunp")
    if not os.path.exists(lunp_file):
        # Try alternative naming
        lunp_files = [f for f in os.listdir(workdir) if f.endswith('_lunp')]
        if lunp_files:
            lunp_file = os.path.join(workdir, lunp_files[0])
        else:
            return None

    # Parse lunp file
    # Format: header lines start with #, then:
    # position<tab>l=1<tab>l=2<tab>...
    # We use l=1 (single nucleotide unpaired probability)
    accessibility = np.zeros(len(sequence))

    with open(lunp_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            try:
                pos = int(parts[0]) - 1  # Convert to 0-based
                # Use l=1 (single nucleotide unpaired probability)
                val = parts[1]
                if val != 'NA' and pos < len(accessibility):
                    accessibility[pos] = float(val)
            except (ValueError, IndexError):
                continue

    return accessibility


def analyze_accessibility(sequences, te_intervals, args):
    """
    Analyze accessibility at TE vs non-TE positions.
    """
    rnaplfold_path = args.rnaplfold_path

    te_accessibilities = []
    non_te_accessibilities = []

    results_by_seq = []

    processed = 0
    for seq_id, sequence in sequences.items():
        if len(sequence) < 50:  # Skip very short sequences
            continue

        if len(sequence) > args.max_length:
            continue

        if processed >= args.max_seqs:
            break

        # Get TE intervals for this sequence
        intervals = te_intervals.get(seq_id, [])
        merged = merge_intervals(intervals)

        # Skip if no TE hits (need both categories)
        if not merged:
            continue

        # Create position mask: True = TE position
        te_mask = np.zeros(len(sequence), dtype=bool)
        for start, end in merged:
            te_mask[max(0, start):min(len(sequence), end)] = True

        # Skip if all TE or all non-TE
        te_count = np.sum(te_mask)
        non_te_count = len(sequence) - te_count
        if te_count < 10 or non_te_count < 10:
            continue

        # Run RNAplfold
        with tempfile.TemporaryDirectory() as workdir:
            accessibility = run_rnaplfold(
                sequence, seq_id,
                span=args.span,
                winsize=args.winsize,
                ulength=args.ulength,
                rnaplfold_path=rnaplfold_path,
                workdir=workdir
            )

        if accessibility is None:
            continue

        # Collect accessibilities
        te_acc = accessibility[te_mask]
        non_te_acc = accessibility[~te_mask]

        te_accessibilities.extend(te_acc)
        non_te_accessibilities.extend(non_te_acc)

        results_by_seq.append({
            'seq_id': seq_id,
            'seq_len': len(sequence),
            'te_positions': int(te_count),
            'non_te_positions': int(non_te_count),
            'te_mean_acc': float(np.mean(te_acc)),
            'non_te_mean_acc': float(np.mean(non_te_acc)),
            'te_median_acc': float(np.median(te_acc)),
            'non_te_median_acc': float(np.median(non_te_acc)),
        })

        processed += 1
        if processed % 10 == 0:
            print(f"  Processed {processed} sequences...", flush=True)

    return te_accessibilities, non_te_accessibilities, results_by_seq


def main():
    parser = argparse.ArgumentParser(
        description='Local RNA structure analysis using RNAplfold'
    )
    parser.add_argument('--utr-fasta', type=Path, required=True,
                        help='Path to UTR FASTA file')
    parser.add_argument('--blast-hits', type=Path, required=True,
                        help='Path to BLAST results TSV')
    parser.add_argument('--span', '-L', type=int, default=150,
                        help='Max base pair span (default: 150, per literature)')
    parser.add_argument('--winsize', '-W', type=int, default=200,
                        help='Window size (default: 200, per literature)')
    parser.add_argument('--ulength', '-u', type=int, default=31,
                        help='Unpaired length parameter (default: 31)')
    parser.add_argument('--max-length', type=int, default=2000,
                        help='Maximum sequence length to analyze (default: 2000)')
    parser.add_argument('--max-seqs', type=int, default=500,
                        help='Maximum sequences to analyze (default: 500)')
    parser.add_argument('--rnaplfold-path', type=str,
                        default='/Users/jacobboysen/miniconda3/envs/bioinformatics-program/bin/RNAplfold',
                        help='Path to RNAplfold binary')
    parser.add_argument('--output-dir', type=Path, default=None,
                        help='Output directory')

    args = parser.parse_args()

    # Setup output directory
    if args.output_dir is None:
        args.output_dir = get_project_root() / 'results' / 'structure_analysis' / 'local_fold'
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Local RNA Structure Analysis (RNAplfold)")
    print("=" * 70)
    print(f"Parameters (per Lange et al. 2012, PMC3384308):")
    print(f"  Max base pair span (L): {args.span}")
    print(f"  Window size (W): {args.winsize}")
    print(f"  Unpaired length (u): {args.ulength}")
    print(f"  Max sequence length: {args.max_length}")
    print(f"  Max sequences: {args.max_seqs}")
    print()

    # Load data
    print("1. Loading UTR sequences...", flush=True)
    sequences = parse_fasta(args.utr_fasta)
    print(f"   Loaded {len(sequences)} sequences")

    print("\n2. Loading TE hit intervals...", flush=True)
    te_intervals = load_te_intervals(args.blast_hits)
    print(f"   Loaded intervals for {len(te_intervals)} transcripts")

    # Filter sequences that have TE hits
    seqs_with_hits = {k: v for k, v in sequences.items() if k in te_intervals}
    print(f"   {len(seqs_with_hits)} sequences have TE hits")

    # Analyze
    print(f"\n3. Running RNAplfold analysis (up to {args.max_seqs} sequences)...", flush=True)
    te_acc, non_te_acc, results = analyze_accessibility(seqs_with_hits, te_intervals, args)

    if not results:
        print("ERROR: No sequences were successfully analyzed")
        return

    # Summary statistics
    print(f"\n4. Computing summary statistics...")
    print(f"   Analyzed {len(results)} sequences")
    print(f"   TE positions: {len(te_acc)}")
    print(f"   Non-TE positions: {len(non_te_acc)}")

    te_acc = np.array(te_acc)
    non_te_acc = np.array(non_te_acc)

    print("\n" + "=" * 70)
    print("RESULTS: Accessibility (probability unpaired)")
    print("=" * 70)
    print(f"\n{'Category':<20} {'Mean':<10} {'Median':<10} {'Std':<10} {'N positions':<15}")
    print("-" * 65)
    print(f"{'TE positions':<20} {np.mean(te_acc):<10.4f} {np.median(te_acc):<10.4f} {np.std(te_acc):<10.4f} {len(te_acc):<15}")
    print(f"{'Non-TE positions':<20} {np.mean(non_te_acc):<10.4f} {np.median(non_te_acc):<10.4f} {np.std(non_te_acc):<10.4f} {len(non_te_acc):<15}")

    # Statistical test
    from scipy import stats
    t_stat, p_value = stats.ttest_ind(te_acc, non_te_acc)
    mw_stat, mw_p = stats.mannwhitneyu(te_acc, non_te_acc, alternative='two-sided')

    print(f"\nStatistical tests:")
    print(f"  t-test: t={t_stat:.2f}, p={p_value:.2e}")
    print(f"  Mann-Whitney U: U={mw_stat:.0f}, p={mw_p:.2e}")

    diff = np.mean(te_acc) - np.mean(non_te_acc)
    print(f"\nDifference (TE - non-TE): {diff:+.4f}")
    if diff > 0:
        print("  → TE positions are MORE accessible (less structured)")
    else:
        print("  → TE positions are LESS accessible (more structured)")

    # Save results
    print(f"\n5. Saving results to {args.output_dir}")

    # Per-sequence results
    with open(args.output_dir / 'local_fold_by_sequence.tsv', 'w') as f:
        headers = ['seq_id', 'seq_len', 'te_positions', 'non_te_positions',
                   'te_mean_acc', 'non_te_mean_acc', 'te_median_acc', 'non_te_median_acc']
        f.write('\t'.join(headers) + '\n')
        for r in results:
            f.write('\t'.join(str(r[h]) for h in headers) + '\n')

    # Summary results
    with open(args.output_dir / 'local_fold_summary.tsv', 'w') as f:
        f.write("category\tn_positions\tmean_accessibility\tmedian_accessibility\tstd_accessibility\n")
        f.write(f"te_positions\t{len(te_acc)}\t{np.mean(te_acc)}\t{np.median(te_acc)}\t{np.std(te_acc)}\n")
        f.write(f"non_te_positions\t{len(non_te_acc)}\t{np.mean(non_te_acc)}\t{np.median(non_te_acc)}\t{np.std(non_te_acc)}\n")

    # Write summary markdown
    with open(args.output_dir / 'LOCAL_FOLD_SUMMARY.md', 'w') as f:
        f.write("# Local RNA Structure Analysis (RNAplfold)\n\n")
        f.write(f"**Date**: 2026-01-18\n")
        f.write(f"**Tool**: ViennaRNA RNAplfold\n")
        f.write(f"**Reference**: Lange et al. 2012 (PMC3384308)\n\n")
        f.write("---\n\n")
        f.write("## Parameters\n\n")
        f.write("| Parameter | Value | Rationale |\n")
        f.write("|-----------|-------|----------|\n")
        f.write(f"| Max bp span (L) | {args.span} | 83% of regulatory bp have span <100nt |\n")
        f.write(f"| Window size (W) | {args.winsize} | L + 50, per literature |\n")
        f.write(f"| Unpaired length (u) | {args.ulength} | Single-nt accessibility |\n")
        f.write(f"| Max seq length | {args.max_length} | Computational limit |\n")
        f.write(f"| Sequences analyzed | {len(results)} | |\n\n")
        f.write("---\n\n")
        f.write("## Key Finding\n\n")
        f.write("| Category | Mean Accessibility | Median | N positions |\n")
        f.write("|----------|-------------------|--------|-------------|\n")
        f.write(f"| TE positions | **{np.mean(te_acc):.4f}** | {np.median(te_acc):.4f} | {len(te_acc):,} |\n")
        f.write(f"| Non-TE positions | **{np.mean(non_te_acc):.4f}** | {np.median(non_te_acc):.4f} | {len(non_te_acc):,} |\n\n")
        f.write(f"**Difference**: {diff:+.4f} (TE - non-TE)\n\n")
        if diff > 0:
            f.write("**Interpretation**: TE positions are MORE accessible (less structured)\n\n")
        else:
            f.write("**Interpretation**: TE positions are LESS accessible (more structured)\n\n")
        f.write("## Statistical Tests\n\n")
        f.write(f"- t-test: t={t_stat:.2f}, p={p_value:.2e}\n")
        f.write(f"- Mann-Whitney U: p={mw_p:.2e}\n\n")
        f.write("---\n\n")
        f.write("## Method\n\n")
        f.write("Unlike global folding (RNAfold on extracted segments), this analysis:\n\n")
        f.write("1. Folds **full UTR sequences** using sliding window approach\n")
        f.write("2. Computes **position-specific accessibility** (P(unpaired) for each nt)\n")
        f.write("3. Maps TE hit coordinates to positions\n")
        f.write("4. Compares accessibility at TE vs non-TE positions\n\n")
        f.write("This avoids arbitrary segment boundaries and border effects.\n")

    print("\nDone!")


if __name__ == '__main__':
    main()
