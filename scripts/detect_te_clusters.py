#!/usr/bin/env python3
"""
Detect TE signal clusters from density data.

Identifies regions with concentrated TE signal by:
1. Finding peaks in the smoothed density signal
2. Merging nearby peaks into clusters
3. Computing cluster statistics and significance
4. Extracting cluster sequences for further analysis
"""

import argparse
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_density_data(density_file):
    """Load density data from pickle file."""
    with open(density_file, 'rb') as f:
        return pickle.load(f)


def find_cluster_regions(density_array, threshold_pct=0.1, merge_distance=50):
    """
    Find contiguous regions of high TE signal.

    Args:
        density_array: 1D numpy array of density values
        threshold_pct: Threshold as fraction of max signal
        merge_distance: Merge regions within this distance

    Returns:
        List of (start, end) tuples for cluster regions
    """
    if len(density_array) == 0 or np.max(density_array) == 0:
        return []

    threshold = np.max(density_array) * threshold_pct

    # Find positions above threshold
    above_threshold = density_array > threshold

    # Find contiguous regions
    regions = []
    in_region = False
    start = 0

    for i, val in enumerate(above_threshold):
        if val and not in_region:
            start = i
            in_region = True
        elif not val and in_region:
            regions.append((start, i))
            in_region = False

    if in_region:
        regions.append((start, len(density_array)))

    # Merge nearby regions
    if len(regions) <= 1:
        return regions

    merged = [regions[0]]
    for region in regions[1:]:
        if region[0] - merged[-1][1] <= merge_distance:
            # Merge
            merged[-1] = (merged[-1][0], region[1])
        else:
            merged.append(region)

    return merged


def annotate_cluster(cluster_start, cluster_end, density_data, sequence=None):
    """
    Annotate a cluster with statistics.

    Args:
        cluster_start: Start position (0-indexed)
        cluster_end: End position (exclusive)
        density_data: Dictionary containing density arrays and hits
        sequence: Optional sequence string for cluster extraction

    Returns:
        Dictionary with cluster annotations
    """
    # Get signal values in cluster region
    combined = density_data['combined'][cluster_start:cluster_end]
    hit_count = density_data['hit_count'][cluster_start:cluster_end]

    # Find contributing TEs
    contributing_tes = defaultdict(lambda: {'count': 0, 'score': 0})

    for hit in density_data.get('hits', []):
        hit_start = hit['qstart'] - 1  # Convert to 0-indexed
        hit_end = hit['qend']

        # Check overlap with cluster
        if hit_end > cluster_start and hit_start < cluster_end:
            te_id = hit['sseqid']
            # Extract TE family (before first underscore or colon)
            if ':' in te_id:
                te_family = te_id.split(':')[0]
            elif '_' in te_id:
                te_family = te_id.split('_')[0]
            else:
                te_family = te_id

            contributing_tes[te_family]['count'] += 1
            contributing_tes[te_family]['score'] += 1.0 / max(hit['evalue'], 1e-180)

    # Sort by score
    sorted_tes = sorted(contributing_tes.items(), key=lambda x: -x[1]['score'])

    annotation = {
        'start': cluster_start + 1,  # 1-indexed for output
        'end': cluster_end,
        'length': cluster_end - cluster_start,
        'max_signal': float(np.max(combined)) if len(combined) > 0 else 0,
        'mean_signal': float(np.mean(combined)) if len(combined) > 0 else 0,
        'total_signal': float(np.sum(combined)) if len(combined) > 0 else 0,
        'max_hit_count': int(np.max(hit_count)) if len(hit_count) > 0 else 0,
        'dominant_te': sorted_tes[0][0] if sorted_tes else None,
        'te_families': {k: {'count': int(v['count']), 'score': float(v['score'])}
                       for k, v in sorted_tes[:5]}
    }

    if sequence:
        annotation['sequence'] = sequence[cluster_start:cluster_end]

    return annotation


def detect_clusters(density_dict, sequences=None, threshold_pct=0.1, merge_distance=50,
                    min_cluster_length=20):
    """
    Detect TE clusters across all sequences.

    Args:
        density_dict: Dictionary of density data per sequence
        sequences: Optional dictionary mapping seq_id to sequence strings
        threshold_pct: Threshold for cluster detection
        merge_distance: Distance for merging nearby regions
        min_cluster_length: Minimum cluster length to report

    Returns:
        Dictionary mapping seq_id to list of cluster annotations
    """
    clusters = {}

    for qseqid, data in density_dict.items():
        # Find cluster regions
        regions = find_cluster_regions(
            data['combined'],
            threshold_pct=threshold_pct,
            merge_distance=merge_distance
        )

        # Annotate each cluster
        seq = sequences.get(qseqid) if sequences else None
        seq_clusters = []

        for start, end in regions:
            if end - start < min_cluster_length:
                continue

            annotation = annotate_cluster(start, end, data, sequence=seq)
            seq_clusters.append(annotation)

        if seq_clusters:
            clusters[qseqid] = seq_clusters

    return clusters


def save_clusters(clusters, output_dir, prefix=''):
    """Save cluster results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save summary TSV
    summary_rows = []
    for qseqid, seq_clusters in clusters.items():
        for i, cluster in enumerate(seq_clusters):
            row = {
                'qseqid': qseqid,
                'cluster_id': i + 1,
                'start': cluster['start'],
                'end': cluster['end'],
                'length': cluster['length'],
                'max_signal': cluster['max_signal'],
                'mean_signal': cluster['mean_signal'],
                'total_signal': cluster['total_signal'],
                'max_hit_count': cluster['max_hit_count'],
                'dominant_te': cluster['dominant_te']
            }
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_file = output_dir / f'{prefix}clusters_summary.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)

    # Save detailed JSON
    json_file = output_dir / f'{prefix}clusters_detail.json'
    # Remove sequences from JSON (too large)
    clusters_for_json = {}
    for qseqid, seq_clusters in clusters.items():
        clusters_for_json[qseqid] = [
            {k: v for k, v in c.items() if k != 'sequence'}
            for c in seq_clusters
        ]
    with open(json_file, 'w') as f:
        json.dump(clusters_for_json, f, indent=2)

    # Save cluster sequences as FASTA
    fasta_records = []
    for qseqid, seq_clusters in clusters.items():
        for i, cluster in enumerate(seq_clusters):
            if 'sequence' in cluster and cluster['sequence']:
                record = SeqRecord(
                    Seq(cluster['sequence']),
                    id=f"{qseqid}_cluster{i+1}",
                    description=f"pos={cluster['start']}-{cluster['end']} "
                                f"signal={cluster['max_signal']:.2e} "
                                f"te={cluster['dominant_te']}"
                )
                fasta_records.append(record)

    if fasta_records:
        fasta_file = output_dir / f'{prefix}cluster_sequences.fasta'
        SeqIO.write(fasta_records, fasta_file, 'fasta')
    else:
        fasta_file = None

    return summary_file, json_file, fasta_file


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'density_dir',
        type=Path,
        nargs='?',
        help='Directory containing density_data.pkl'
    )
    parser.add_argument(
        '--query-fasta',
        type=Path,
        default=Path('data/queries/germ_plasm/3UTR_sense_tier1.fasta'),
        help='Query FASTA for extracting cluster sequences'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/te_clusters'),
        help='Output directory (default: data/te_clusters)'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.1,
        help='Cluster detection threshold as fraction of max (default: 0.1)'
    )
    parser.add_argument(
        '--merge-distance',
        type=int,
        default=50,
        help='Merge clusters within this distance (default: 50bp)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=20,
        help='Minimum cluster length (default: 20bp)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find density data
    if args.density_dir is None:
        # Try to find from most recent parameter sweep
        sweep_dir = Path('results/parameter_sweep')
        if sweep_dir.exists():
            subdirs = sorted([d for d in sweep_dir.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                combo_dirs = sorted([d for d in subdirs[0].iterdir()
                                    if d.is_dir() and d.name.startswith('combo_')])
                if combo_dirs:
                    density_dir = combo_dirs[0] / 'density'
                    if density_dir.exists():
                        args.density_dir = density_dir

        if args.density_dir is None:
            print("Error: No density directory specified and none found", file=sys.stderr)
            return 1

    density_file = args.density_dir / 'density_data.pkl'
    if not density_file.exists():
        print(f"Error: Density file not found: {density_file}", file=sys.stderr)
        return 1

    print("TE Cluster Detection")
    print("=" * 60)
    print(f"Density data: {args.density_dir}")
    print(f"Output: {args.output_dir}")
    print(f"Threshold: {args.threshold}")
    print(f"Merge distance: {args.merge_distance}bp")
    print()

    # Load density data
    print("Loading density data...")
    density_dict = load_density_data(density_file)
    print(f"  Loaded data for {len(density_dict)} sequences")

    # Load sequences if available
    sequences = {}
    if args.query_fasta.exists():
        print(f"\nLoading sequences from: {args.query_fasta}")
        for record in SeqIO.parse(args.query_fasta, 'fasta'):
            sequences[record.id] = str(record.seq)
        print(f"  Loaded {len(sequences)} sequences")

    # Detect clusters
    print("\nDetecting clusters...")
    clusters = detect_clusters(
        density_dict,
        sequences=sequences,
        threshold_pct=args.threshold,
        merge_distance=args.merge_distance,
        min_cluster_length=args.min_length
    )

    total_clusters = sum(len(c) for c in clusters.values())
    print(f"  Found {total_clusters} clusters in {len(clusters)} sequences")

    if total_clusters == 0:
        print("\nNo clusters found with current parameters.")
        return 0

    # Save results
    print("\nSaving results...")
    summary_file, json_file, fasta_file = save_clusters(clusters, args.output_dir)
    print(f"  Saved: {summary_file.name}")
    print(f"  Saved: {json_file.name}")
    if fasta_file:
        print(f"  Saved: {fasta_file.name}")

    # Print summary
    print("\n" + "=" * 60)
    print("Cluster Summary:")
    print("-" * 60)
    print(f"{'Query':<30} {'Clusters':>10} {'Max Signal':>12} {'Dominant TE':<20}")
    print("-" * 60)

    for qseqid in sorted(clusters.keys()):
        seq_clusters = clusters[qseqid]
        max_signal = max(c['max_signal'] for c in seq_clusters)
        dominant_tes = [c['dominant_te'] for c in seq_clusters if c['dominant_te']]
        dominant = dominant_tes[0] if dominant_tes else 'N/A'
        print(f"{qseqid:<30} {len(seq_clusters):>10} {max_signal:>12.2e} {dominant:<20}")

    print("\n" + "=" * 60)
    print("Cluster detection complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
