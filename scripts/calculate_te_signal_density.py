#!/usr/bin/env python3
"""
Calculate position-wise TE signal density from BLAST results.

Computes multiple density metrics for each 3'UTR:
1. Hit count density: Number of overlapping HSPs at each position
2. E-value weighted: sum(1/e-value) at each position
3. Bitscore weighted: sum(bitscore/length) at each position
4. Combined score: sum(pident/100 * 1/e-value * 1/length) at each position

Results are smoothed with a Gaussian kernel.
"""

import argparse
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


# BLAST output columns (matching te_parameter_sweep.py extended format)
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'slen', 'qseq', 'sseq'
]


def load_blast_results(results_file):
    """Load BLAST results from TSV file."""
    if not results_file.exists():
        raise FileNotFoundError(f"Results file not found: {results_file}")

    if results_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS)

    # Ensure qstart < qend (BLAST can report reverse)
    swap_mask = df['qstart'] > df['qend']
    df.loc[swap_mask, ['qstart', 'qend']] = df.loc[swap_mask, ['qend', 'qstart']].values

    return df


def calculate_density_arrays(df, query_lengths):
    """
    Calculate position-wise density arrays for each query sequence.

    Returns:
        Dictionary mapping query_id to density_data dict containing:
        - length: sequence length
        - hit_count: raw overlap count at each position
        - evalue_weighted: sum(1/evalue) at each position
        - bitscore_weighted: sum(bitscore/length) at each position
        - combined: sum(pident/100 * 1/evalue * 1/length) at each position
        - hits: list of HSP details
    """
    results = {}

    # Group by query
    for qseqid, group in df.groupby('qseqid'):
        # Get query length
        qlen = group['qlen'].iloc[0] if 'qlen' in group.columns else query_lengths.get(qseqid, 0)

        if qlen == 0:
            continue

        # Initialize arrays
        hit_count = np.zeros(qlen, dtype=np.float64)
        evalue_weighted = np.zeros(qlen, dtype=np.float64)
        bitscore_weighted = np.zeros(qlen, dtype=np.float64)
        combined = np.zeros(qlen, dtype=np.float64)

        # Store hit details
        hits = []

        for _, row in group.iterrows():
            start = int(row['qstart']) - 1  # Convert to 0-indexed
            end = int(row['qend'])  # Keep as exclusive end

            # Clamp to valid range
            start = max(0, start)
            end = min(qlen, end)

            if start >= end:
                continue

            evalue = max(row['evalue'], 1e-180)  # Avoid division by zero
            pident = row['pident']
            bitscore = row['bitscore']
            hsp_length = end - start

            # Accumulate at each position
            hit_count[start:end] += 1
            evalue_weighted[start:end] += 1.0 / evalue
            bitscore_weighted[start:end] += bitscore / hsp_length
            combined[start:end] += (pident / 100.0) * (1.0 / evalue) * (1.0 / hsp_length)

            # Store hit details
            hits.append({
                'sseqid': row['sseqid'],
                'qstart': int(row['qstart']),
                'qend': int(row['qend']),
                'pident': pident,
                'evalue': row['evalue'],
                'bitscore': bitscore,
                'length': hsp_length
            })

        results[qseqid] = {
            'length': qlen,
            'hit_count': hit_count,
            'evalue_weighted': evalue_weighted,
            'bitscore_weighted': bitscore_weighted,
            'combined': combined,
            'hits': hits,
            'total_hits': len(hits)
        }

    return results


def smooth_density(density_dict, sigma=25):
    """Apply Gaussian smoothing to density arrays."""
    metrics = ['hit_count', 'evalue_weighted', 'bitscore_weighted', 'combined']

    for qseqid, data in density_dict.items():
        for metric in metrics:
            if metric in data and len(data[metric]) > 0:
                # Store raw and smoothed versions
                data[f'{metric}_raw'] = data[metric].copy()
                data[metric] = gaussian_filter1d(data[metric], sigma=sigma)

    return density_dict


def find_peaks(density_array, threshold_pct=0.1, min_distance=10):
    """
    Find peaks in density array.

    Args:
        density_array: 1D numpy array of density values
        threshold_pct: Minimum peak height as percentile of max
        min_distance: Minimum distance between peaks

    Returns:
        List of peak positions (0-indexed)
    """
    if len(density_array) == 0 or np.max(density_array) == 0:
        return []

    threshold = np.max(density_array) * threshold_pct

    # Find local maxima
    peaks = []
    for i in range(1, len(density_array) - 1):
        if density_array[i] > threshold:
            if density_array[i] > density_array[i-1] and density_array[i] >= density_array[i+1]:
                peaks.append(i)

    # Filter by minimum distance
    if len(peaks) <= 1:
        return peaks

    filtered = [peaks[0]]
    for p in peaks[1:]:
        if p - filtered[-1] >= min_distance:
            filtered.append(p)

    return filtered


def annotate_peaks(density_dict, metric='combined', threshold_pct=0.1):
    """Add peak annotations to density dictionary."""
    for qseqid, data in density_dict.items():
        if metric not in data:
            data['peaks'] = []
            continue

        peaks = find_peaks(data[metric], threshold_pct=threshold_pct)

        peak_annotations = []
        for peak_pos in peaks:
            # Find contributing TEs at this peak
            contributing_tes = defaultdict(float)
            for hit in data['hits']:
                if hit['qstart'] - 1 <= peak_pos < hit['qend']:
                    te_family = hit['sseqid'].split(':')[0] if ':' in hit['sseqid'] else hit['sseqid']
                    contributing_tes[te_family] += 1.0 / max(hit['evalue'], 1e-180)

            # Sort by contribution
            sorted_tes = sorted(contributing_tes.items(), key=lambda x: -x[1])

            peak_annotations.append({
                'position': peak_pos + 1,  # 1-indexed for output
                'value': float(data[metric][peak_pos]),
                'contributing_tes': dict(sorted_tes[:5])  # Top 5
            })

        data['peaks'] = peak_annotations

    return density_dict


def convert_to_json_serializable(obj):
    """Convert numpy types to Python native types for JSON serialization."""
    if isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(i) for i in obj]
    return obj


def save_results(density_dict, output_dir):
    """Save density results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save full data as pickle (for downstream analysis)
    pickle_file = output_dir / 'density_data.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(density_dict, f)

    # Save summary TSV
    summary_rows = []
    for qseqid, data in density_dict.items():
        row = {
            'qseqid': qseqid,
            'length': int(data['length']),
            'total_hits': int(data['total_hits']),
            'num_peaks': len(data.get('peaks', [])),
            'max_hit_count': float(np.max(data['hit_count'])) if len(data['hit_count']) > 0 else 0,
            'max_evalue_weighted': float(np.max(data['evalue_weighted'])) if len(data['evalue_weighted']) > 0 else 0,
            'max_combined': float(np.max(data['combined'])) if len(data['combined']) > 0 else 0
        }
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_file = output_dir / 'density_summary.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)

    # Save peak details as JSON (convert numpy types)
    peaks_data = {}
    for qseqid, data in density_dict.items():
        peaks_data[qseqid] = {
            'length': int(data['length']),
            'peaks': convert_to_json_serializable(data.get('peaks', []))
        }

    peaks_file = output_dir / 'peaks.json'
    with open(peaks_file, 'w') as f:
        json.dump(peaks_data, f, indent=2)

    return pickle_file, summary_file, peaks_file


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'blast_results',
        type=Path,
        nargs='?',
        help='BLAST results TSV file (or directory containing blast_results.tsv)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Output directory (default: same as input with _density suffix)'
    )
    parser.add_argument(
        '--sigma',
        type=float,
        default=25.0,
        help='Gaussian smoothing sigma (default: 25)'
    )
    parser.add_argument(
        '--peak-threshold',
        type=float,
        default=0.1,
        help='Peak threshold as fraction of max (default: 0.1)'
    )
    parser.add_argument(
        '--query-fasta',
        type=Path,
        default=Path('data/queries/germ_plasm/3UTR_sense_tier1.fasta'),
        help='Query FASTA for sequence lengths (default: data/queries/germ_plasm/3UTR_sense_tier1.fasta)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find BLAST results file
    if args.blast_results is None:
        # Try to find most recent sweep results
        sweep_dir = Path('results/parameter_sweep')
        if sweep_dir.exists():
            subdirs = sorted([d for d in sweep_dir.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                # Use first combo's results from most recent sweep
                combo_dirs = sorted([d for d in subdirs[0].iterdir() if d.is_dir() and d.name.startswith('combo_')])
                if combo_dirs:
                    args.blast_results = combo_dirs[0] / 'blast_results.tsv'

        if args.blast_results is None:
            print("Error: No BLAST results specified and none found in results/parameter_sweep/",
                  file=sys.stderr)
            return 1

    # Handle directory input
    if args.blast_results.is_dir():
        args.blast_results = args.blast_results / 'blast_results.tsv'

    if not args.blast_results.exists():
        print(f"Error: BLAST results not found: {args.blast_results}", file=sys.stderr)
        return 1

    # Set output directory
    if args.output_dir is None:
        args.output_dir = args.blast_results.parent / 'density'

    print("TE Signal Density Calculator")
    print("=" * 60)
    print(f"Input: {args.blast_results}")
    print(f"Output: {args.output_dir}")
    print(f"Smoothing sigma: {args.sigma}")
    print()

    # Load query lengths from FASTA if available
    query_lengths = {}
    if args.query_fasta.exists():
        from Bio import SeqIO
        for record in SeqIO.parse(args.query_fasta, 'fasta'):
            query_lengths[record.id] = len(record.seq)
        if args.verbose:
            print(f"Loaded {len(query_lengths)} query lengths from FASTA")

    # Load BLAST results
    print("Loading BLAST results...")
    df = load_blast_results(args.blast_results)
    print(f"  Loaded {len(df)} hits for {df['qseqid'].nunique()} queries")

    if len(df) == 0:
        print("No hits found. Exiting.")
        return 0

    # Calculate density arrays
    print("\nCalculating density arrays...")
    density_dict = calculate_density_arrays(df, query_lengths)
    print(f"  Processed {len(density_dict)} sequences")

    # Apply smoothing
    print(f"\nApplying Gaussian smoothing (sigma={args.sigma})...")
    density_dict = smooth_density(density_dict, sigma=args.sigma)

    # Find and annotate peaks
    print(f"\nFinding peaks (threshold={args.peak_threshold})...")
    density_dict = annotate_peaks(density_dict, threshold_pct=args.peak_threshold)

    total_peaks = sum(len(d.get('peaks', [])) for d in density_dict.values())
    print(f"  Found {total_peaks} peaks across all sequences")

    # Save results
    print("\nSaving results...")
    pickle_file, summary_file, peaks_file = save_results(density_dict, args.output_dir)
    print(f"  Saved: {pickle_file.name}")
    print(f"  Saved: {summary_file.name}")
    print(f"  Saved: {peaks_file.name}")

    # Print summary
    print("\n" + "=" * 60)
    print("Summary by Sequence:")
    print("-" * 60)
    print(f"{'Query':<35} {'Length':>8} {'Hits':>6} {'Peaks':>6} {'Max Signal':>12}")
    print("-" * 60)

    for qseqid in sorted(density_dict.keys()):
        data = density_dict[qseqid]
        max_signal = np.max(data['combined']) if len(data['combined']) > 0 else 0
        print(f"{qseqid:<35} {data['length']:>8} {data['total_hits']:>6} "
              f"{len(data.get('peaks', [])):>6} {max_signal:>12.2e}")

    print("\n" + "=" * 60)
    print("Density calculation complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
