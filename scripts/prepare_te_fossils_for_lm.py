#!/usr/bin/env python3
"""
Prepare TE fossil dataset for genomic language model scoring.

Reads annotated BLAST hits from the regulatory TE analysis pipeline and
produces a clean dataset bundling:
  1. Full regulatory region sequences (FASTA)
  2. Per-nucleotide TE masks (NumPy compressed)
  3. Per-hit annotations with TE family/class metadata
  4. TE-free regions as matched negative controls

Quality tiers (cumulative):
  strict:   e-value <= 1e-5  (mask=3)
  moderate: e-value <= 1e-3  (mask=2)
  relaxed:  e-value <= 0.01  (mask=1)

Usage:
    python scripts/prepare_te_fossils_for_lm.py -v
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir, get_te_database_path
from utils.data_loaders import parse_fasta, load_te_database, parse_te_class


def load_promoter_metadata(path):
    """Load promoter metadata into standardized dict keyed by region_id."""
    df = pd.read_csv(path, sep='\t')
    metadata = {}
    for _, row in df.iterrows():
        rid = row['promoter_id']
        metadata[rid] = {
            'region_type': 'promoter',
            'chrom': row['chrom'],
            'start': int(row['region_start']),
            'end': int(row['region_end']),
            'strand': row['strand'],
            'length': int(row['length']),
            'fbgn': row['fbgn'] if pd.notna(row.get('fbgn')) else '',
            'symbol': row['gene_symbol'] if pd.notna(row.get('gene_symbol')) else '',
            'source': row['source'] if pd.notna(row.get('source')) else '',
            'is_starr_seq': False,
        }
    return metadata


def load_enhancer_metadata(path):
    """Load enhancer metadata into standardized dict keyed by region_id."""
    df = pd.read_csv(path, sep='\t')
    metadata = {}
    for _, row in df.iterrows():
        rid = row['enhancer_id']
        # Prefer assoc_fbgn/assoc_symbol, fall back to nearest
        fbgn = ''
        symbol = ''
        if pd.notna(row.get('assoc_fbgn')) and row['assoc_fbgn']:
            fbgn = row['assoc_fbgn']
            symbol = row['assoc_symbol'] if pd.notna(row.get('assoc_symbol')) else ''
        elif pd.notna(row.get('nearest_fbgn')) and row['nearest_fbgn']:
            fbgn = row['nearest_fbgn']
            symbol = row['nearest_symbol'] if pd.notna(row.get('nearest_symbol')) else ''

        metadata[rid] = {
            'region_type': 'enhancer',
            'chrom': row['chrom'],
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': '.',
            'length': int(row['length']),
            'fbgn': fbgn,
            'symbol': symbol,
            'source': row['source'] if pd.notna(row.get('source')) else '',
            'is_starr_seq': bool(row.get('is_starr_seq', False))
                           if pd.notna(row.get('is_starr_seq')) else False,
        }
    return metadata


def classify_tier(evalue, e_strict, e_moderate, e_relaxed):
    """Return (tier_name, tier_value) for an e-value, or None if above relaxed."""
    if evalue <= e_strict:
        return 'strict', 3
    elif evalue <= e_moderate:
        return 'moderate', 2
    elif evalue <= e_relaxed:
        return 'relaxed', 1
    return None, 0


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--promoter-hits',
        type=Path,
        default=get_results_dir() / 'regulatory_analysis' / 'promoter' / 'promoter_te_hits_annotated.tsv',
    )
    parser.add_argument(
        '--enhancer-hits',
        type=Path,
        default=get_results_dir() / 'regulatory_analysis' / 'enhancer' / 'enhancer_te_hits_annotated.tsv',
    )
    parser.add_argument(
        '--promoter-meta',
        type=Path,
        default=get_project_root() / 'data' / 'queries' / 'regulatory' / 'promoters' / 'promoter_metadata.tsv',
    )
    parser.add_argument(
        '--enhancer-meta',
        type=Path,
        default=get_project_root() / 'data' / 'queries' / 'regulatory' / 'enhancers' / 'enhancer_metadata.tsv',
    )
    parser.add_argument(
        '--promoter-fasta',
        type=Path,
        default=get_project_root() / 'data' / 'queries' / 'regulatory' / 'promoters' / 'promoters_sense.fasta',
    )
    parser.add_argument(
        '--enhancer-fasta',
        type=Path,
        default=get_project_root() / 'data' / 'queries' / 'regulatory' / 'enhancers' / 'enhancers_sense.fasta',
    )
    parser.add_argument(
        '--te-database',
        type=Path,
        default=get_te_database_path(),
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=get_results_dir() / 'te_fossil_lm_dataset',
    )
    parser.add_argument('--evalue-strict', type=float, default=1e-5)
    parser.add_argument('--evalue-moderate', type=float, default=1e-3)
    parser.add_argument('--evalue-relaxed', type=float, default=0.01)
    parser.add_argument('--chunk-size', type=int, default=500000)
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    print("Prepare TE Fossil Dataset for Language Models")
    print("=" * 60)

    # ── Step 1: Load small files ─────────────────────────────────
    if args.verbose:
        print("\nStep 1: Loading metadata, sequences, and TE database...")

    te_db = load_te_database(args.te_database)
    if args.verbose:
        print(f"  TE database: {len(te_db):,} entries")

    # Load metadata
    promoter_meta = load_promoter_metadata(args.promoter_meta)
    enhancer_meta = load_enhancer_metadata(args.enhancer_meta)
    if args.verbose:
        print(f"  Promoter metadata: {len(promoter_meta):,} regions")
        print(f"  Enhancer metadata: {len(enhancer_meta):,} regions")

    # Load sequences
    promoter_seqs = parse_fasta(args.promoter_fasta)
    enhancer_seqs = parse_fasta(args.enhancer_fasta)
    if args.verbose:
        print(f"  Promoter sequences: {len(promoter_seqs):,}")
        print(f"  Enhancer sequences: {len(enhancer_seqs):,}")

    # Validate: keep only IDs present in both metadata and FASTA
    all_metadata = {}
    all_sequences = {}

    for rid in promoter_meta:
        if rid in promoter_seqs:
            all_metadata[rid] = promoter_meta[rid]
            all_sequences[rid] = promoter_seqs[rid]

    for rid in enhancer_meta:
        if rid in enhancer_seqs:
            all_metadata[rid] = enhancer_meta[rid]
            all_sequences[rid] = enhancer_seqs[rid]

    n_promoters = sum(1 for m in all_metadata.values() if m['region_type'] == 'promoter')
    n_enhancers = sum(1 for m in all_metadata.values() if m['region_type'] == 'enhancer')
    print(f"\nValid regions: {len(all_metadata):,} ({n_promoters:,} promoters, {n_enhancers:,} enhancers)")

    # Free raw dicts
    del promoter_meta, enhancer_meta, promoter_seqs, enhancer_seqs

    # ── Step 2: Allocate masks and counters ──────────────────────
    if args.verbose:
        print("\nStep 2: Allocating masks and hit counters...")

    masks = {rid: np.zeros(len(all_sequences[rid]), dtype=np.uint8) for rid in all_metadata}
    hit_counts = {rid: {'strict': 0, 'moderate': 0, 'relaxed': 0} for rid in all_metadata}

    total_mask_bytes = sum(m.nbytes for m in masks.values())
    if args.verbose:
        print(f"  Mask memory: {total_mask_bytes / 1e6:.1f} MB")

    # Prepare output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    hits_path = args.output_dir / 'hits.tsv'

    hit_columns = [
        'region_id', 'region_type', 'hit_start', 'hit_end',
        'genomic_start', 'genomic_end', 'te_id', 'te_family', 'te_class',
        'pident', 'length', 'evalue', 'bitscore',
        'in_repeatmasker', 'rm_name', 'rm_class', 'quality_tier',
    ]

    # ── Step 3: Stream annotated hits ────────────────────────────
    if args.verbose:
        print("\nStep 3: Processing annotated hits...")

    total_hits_written = 0

    hit_sources = [
        ('promoter', args.promoter_hits),
        ('enhancer', args.enhancer_hits),
    ]

    with open(hits_path, 'w') as hits_fh:
        hits_fh.write('\t'.join(hit_columns) + '\n')

        for region_type, hits_file in hit_sources:
            if not hits_file.exists():
                print(f"  WARNING: {hits_file} not found, skipping {region_type}")
                continue

            if args.verbose:
                print(f"  Reading {region_type} hits from {hits_file}...")

            chunk_iter = pd.read_csv(
                hits_file, sep='\t', chunksize=args.chunk_size,
                dtype={'chrom': str, 'fbgn': str, 'symbol': str,
                       'rm_name': str, 'rm_class': str},
            )
            type_hits = 0

            for chunk in chunk_iter:
                for _, row in chunk.iterrows():
                    rid = row['region_id']

                    # Skip if not in valid set
                    if rid not in all_metadata:
                        continue

                    evalue = float(row['evalue'])

                    # Skip if above relaxed threshold
                    if evalue > args.evalue_relaxed:
                        continue

                    tier_name, tier_value = classify_tier(
                        evalue, args.evalue_strict, args.evalue_moderate, args.evalue_relaxed
                    )
                    if tier_value == 0:
                        continue

                    # Convert query coordinates to 0-based half-open
                    qstart = int(row['query_start'])
                    qend = int(row['query_end'])
                    hit_start = qstart - 1  # 0-based
                    hit_end = qend          # exclusive

                    # Clamp to sequence bounds
                    seq_len = len(all_sequences[rid])
                    hit_start = max(0, hit_start)
                    hit_end = min(seq_len, hit_end)

                    if hit_start >= hit_end:
                        continue

                    # Update mask (take max = best tier)
                    mask = masks[rid]
                    mask[hit_start:hit_end] = np.maximum(
                        mask[hit_start:hit_end], tier_value
                    )

                    # Increment hit counters (cumulative)
                    if tier_name == 'strict':
                        hit_counts[rid]['strict'] += 1
                        hit_counts[rid]['moderate'] += 1
                        hit_counts[rid]['relaxed'] += 1
                    elif tier_name == 'moderate':
                        hit_counts[rid]['moderate'] += 1
                        hit_counts[rid]['relaxed'] += 1
                    else:  # relaxed
                        hit_counts[rid]['relaxed'] += 1

                    # Look up TE family/class
                    te_id = row['te_id']
                    te_family = ''
                    te_class = ''
                    if te_id in te_db:
                        te_family = te_db[te_id]['name']
                        te_class = te_db[te_id]['class']
                    elif pd.notna(row.get('rm_name')) and row['rm_name']:
                        te_family = str(row['rm_name'])
                        te_class = str(row['rm_class']) if pd.notna(row.get('rm_class')) else ''
                    if not te_class and te_family:
                        te_class = parse_te_class(te_family)

                    # Write hit row
                    in_rm = row['in_repeatmasker']
                    rm_name = str(row['rm_name']) if pd.notna(row.get('rm_name')) else ''
                    rm_class = str(row['rm_class']) if pd.notna(row.get('rm_class')) else ''

                    hit_row = [
                        rid,
                        all_metadata[rid]['region_type'],
                        str(hit_start),
                        str(hit_end),
                        str(row['genomic_start']),
                        str(row['genomic_end']),
                        te_id,
                        te_family,
                        te_class,
                        str(row['pident']),
                        str(row['length']),
                        str(evalue),
                        str(row['bitscore']),
                        str(in_rm),
                        rm_name,
                        rm_class,
                        tier_name,
                    ]
                    hits_fh.write('\t'.join(hit_row) + '\n')
                    type_hits += 1
                    total_hits_written += 1

            if args.verbose:
                print(f"    {region_type}: {type_hits:,} hits written")

    print(f"\nTotal hits written: {total_hits_written:,}")

    # ── Step 4: Write regions.tsv ────────────────────────────────
    if args.verbose:
        print("\nStep 4: Writing regions.tsv...")

    regions_path = args.output_dir / 'regions.tsv'
    region_columns = [
        'region_id', 'region_type', 'chrom', 'start', 'end', 'strand', 'length',
        'fbgn', 'symbol', 'source', 'is_starr_seq',
        'n_te_hits_strict', 'n_te_hits_moderate', 'n_te_hits_relaxed',
        'te_bp_strict', 'te_bp_moderate', 'te_bp_relaxed',
        'te_fraction_strict', 'te_fraction_moderate', 'te_fraction_relaxed',
        'has_te_strict', 'has_te_moderate', 'has_te_relaxed',
    ]

    # Collect region rows in sorted order for deterministic output
    sorted_rids = sorted(all_metadata.keys())

    with open(regions_path, 'w') as f:
        f.write('\t'.join(region_columns) + '\n')

        for rid in sorted_rids:
            meta = all_metadata[rid]
            mask = masks[rid]
            counts = hit_counts[rid]
            length = meta['length']

            te_bp_strict = int(np.sum(mask >= 3))
            te_bp_moderate = int(np.sum(mask >= 2))
            te_bp_relaxed = int(np.sum(mask >= 1))

            te_frac_strict = te_bp_strict / length if length > 0 else 0.0
            te_frac_moderate = te_bp_moderate / length if length > 0 else 0.0
            te_frac_relaxed = te_bp_relaxed / length if length > 0 else 0.0

            row = [
                rid,
                meta['region_type'],
                meta['chrom'],
                str(meta['start']),
                str(meta['end']),
                meta['strand'],
                str(length),
                meta['fbgn'],
                meta['symbol'],
                meta['source'],
                str(meta['is_starr_seq']),
                str(counts['strict']),
                str(counts['moderate']),
                str(counts['relaxed']),
                str(te_bp_strict),
                str(te_bp_moderate),
                str(te_bp_relaxed),
                f"{te_frac_strict:.6f}",
                f"{te_frac_moderate:.6f}",
                f"{te_frac_relaxed:.6f}",
                str(counts['strict'] > 0),
                str(counts['moderate'] > 0),
                str(counts['relaxed'] > 0),
            ]
            f.write('\t'.join(row) + '\n')

    print(f"Regions written: {len(sorted_rids):,}")

    # ── Step 5: Write sequences.fasta ────────────────────────────
    if args.verbose:
        print("\nStep 5: Writing sequences.fasta...")

    fasta_path = args.output_dir / 'sequences.fasta'
    with open(fasta_path, 'w') as f:
        for rid in sorted_rids:
            meta = all_metadata[rid]
            mask = masks[rid]
            length = meta['length']
            te_frac = int(np.sum(mask >= 3)) / length if length > 0 else 0.0

            header = (
                f">{rid} "
                f"region_type={meta['region_type']} "
                f"chrom={meta['chrom']} "
                f"start={meta['start']} "
                f"end={meta['end']} "
                f"strand={meta['strand']} "
                f"gene={meta['symbol'] or meta['fbgn'] or 'NA'} "
                f"te_fraction_strict={te_frac:.4f}"
            )
            f.write(header + '\n')

            # Write sequence in 80-char lines
            seq = all_sequences[rid]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f"FASTA entries written: {len(sorted_rids):,}")

    # ── Step 6: Save masks.npz ───────────────────────────────────
    if args.verbose:
        print("\nStep 6: Saving masks.npz...")

    masks_path = args.output_dir / 'masks.npz'
    # Only store regions with at least one TE-derived position
    nonzero_masks = {rid: mask for rid, mask in masks.items() if mask.any()}
    np.savez_compressed(str(masks_path), **nonzero_masks)
    print(f"Masks saved: {len(nonzero_masks):,} regions with TE hits (of {len(masks):,} total)")

    # ── Step 7: Write dataset_stats.json ─────────────────────────
    if args.verbose:
        print("\nStep 7: Computing dataset statistics...")

    n_te_bearing = {
        'strict': sum(1 for c in hit_counts.values() if c['strict'] > 0),
        'moderate': sum(1 for c in hit_counts.values() if c['moderate'] > 0),
        'relaxed': sum(1 for c in hit_counts.values() if c['relaxed'] > 0),
    }
    n_te_free = len(all_metadata) - n_te_bearing['relaxed']

    total_bp = sum(len(all_sequences[rid]) for rid in all_metadata)
    te_bp_strict_total = sum(int(np.sum(masks[rid] >= 3)) for rid in all_metadata)
    te_bp_moderate_total = sum(int(np.sum(masks[rid] >= 2)) for rid in all_metadata)
    te_bp_relaxed_total = sum(int(np.sum(masks[rid] >= 1)) for rid in all_metadata)

    total_hits_by_tier = {
        'strict': sum(c['strict'] for c in hit_counts.values()),
        'moderate': sum(c['moderate'] for c in hit_counts.values()),
        'relaxed': sum(c['relaxed'] for c in hit_counts.values()),
    }

    stats = {
        'n_regions': len(all_metadata),
        'n_promoters': n_promoters,
        'n_enhancers': n_enhancers,
        'n_te_bearing_strict': n_te_bearing['strict'],
        'n_te_bearing_moderate': n_te_bearing['moderate'],
        'n_te_bearing_relaxed': n_te_bearing['relaxed'],
        'n_te_free': n_te_free,
        'total_hits_written': total_hits_written,
        'total_hits_strict': total_hits_by_tier['strict'],
        'total_hits_moderate': total_hits_by_tier['moderate'],
        'total_hits_relaxed': total_hits_by_tier['relaxed'],
        'total_bp': total_bp,
        'te_bp_strict': te_bp_strict_total,
        'te_bp_moderate': te_bp_moderate_total,
        'te_bp_relaxed': te_bp_relaxed_total,
        'te_fraction_strict': te_bp_strict_total / total_bp if total_bp > 0 else 0,
        'te_fraction_moderate': te_bp_moderate_total / total_bp if total_bp > 0 else 0,
        'te_fraction_relaxed': te_bp_relaxed_total / total_bp if total_bp > 0 else 0,
        'evalue_thresholds': {
            'strict': args.evalue_strict,
            'moderate': args.evalue_moderate,
            'relaxed': args.evalue_relaxed,
        },
        'output_files': {
            'regions': 'regions.tsv',
            'hits': 'hits.tsv',
            'sequences': 'sequences.fasta',
            'masks': 'masks.npz',
        },
    }

    stats_path = args.output_dir / 'dataset_stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)

    # ── Summary ──────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Dataset Summary")
    print("-" * 60)
    print(f"  Regions:          {stats['n_regions']:,}")
    print(f"    Promoters:      {stats['n_promoters']:,}")
    print(f"    Enhancers:      {stats['n_enhancers']:,}")
    print(f"  TE-bearing (strict):   {n_te_bearing['strict']:,}")
    print(f"  TE-bearing (moderate): {n_te_bearing['moderate']:,}")
    print(f"  TE-bearing (relaxed):  {n_te_bearing['relaxed']:,}")
    print(f"  TE-free controls:      {n_te_free:,}")
    print(f"  Total hits:       {total_hits_written:,}")
    print(f"  Total bp:         {total_bp:,}")
    print(f"  TE bp (strict):   {te_bp_strict_total:,} ({stats['te_fraction_strict']:.4%})")
    print(f"  TE bp (moderate): {te_bp_moderate_total:,} ({stats['te_fraction_moderate']:.4%})")
    print(f"  TE bp (relaxed):  {te_bp_relaxed_total:,} ({stats['te_fraction_relaxed']:.4%})")
    print(f"\nOutput directory: {args.output_dir}")
    print(f"  regions.tsv:       {regions_path}")
    print(f"  hits.tsv:          {hits_path}")
    print(f"  sequences.fasta:   {fasta_path}")
    print(f"  masks.npz:         {masks_path}")
    print(f"  dataset_stats.json:{stats_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
