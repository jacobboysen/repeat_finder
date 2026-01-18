#!/usr/bin/env python3
"""
Unified TE hit deduplication for exons, 3'UTRs, and 5'UTRs.

When the same genomic region from different transcript isoforms produces
identical BLAST hits (same TE, same TE coordinates), count it only once.

Deduplication key: (fbgn, chrom, genomic_start, genomic_end, te_id, sstart, send)

Usage:
    python scripts/deduplicate_te_hits.py --seq-type exon
    python scripts/deduplicate_te_hits.py --seq-type 3utr
    python scripts/deduplicate_te_hits.py --seq-type 5utr
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir, get_references_dir


def load_exon_metadata(path):
    """Load exon metadata with genomic coordinates."""
    meta = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            meta[row['exon_id']] = {
                'fbgn': row['fbgn'],
                'fbtr': row['fbtr'],
                'chrom': row['chrom'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
            }
    return meta


def load_utr_metadata_from_gff(gff_path, utr_type):
    """
    Load UTR metadata from GFF file.

    Args:
        gff_path: Path to GFF file
        utr_type: 'three_prime_UTR' or 'five_prime_UTR'

    Returns:
        dict mapping FBtr ID to metadata
    """
    # First pass: build FBtr -> FBgn mapping from mRNA records
    fbtr_to_fbgn = {}

    print(f"  Parsing GFF for transcript-gene mapping...")
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'mRNA':
                attrs = parts[8]
                id_match = re.search(r'ID=([^;]+)', attrs)
                parent_match = re.search(r'Parent=([^;]+)', attrs)

                if id_match and parent_match:
                    fbtr = id_match.group(1)
                    fbgn = parent_match.group(1)
                    fbtr_to_fbgn[fbtr] = fbgn

    print(f"    Found {len(fbtr_to_fbgn):,} transcript-gene mappings")

    # Second pass: collect UTR records
    # Note: A transcript can have multiple UTR segments (due to introns in UTR)
    # We'll merge them into a single span for deduplication purposes
    utr_segments = defaultdict(list)

    print(f"  Parsing GFF for {utr_type} records...")
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == utr_type:
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]

                # Handle comma-separated Parent values (shared UTRs)
                parent_match = re.search(r'Parent=([^;]+)', attrs)
                if parent_match:
                    parent_ids = parent_match.group(1).split(',')
                    for fbtr in parent_ids:
                        utr_segments[fbtr].append({
                            'chrom': chrom,
                            'start': start,
                            'end': end,
                            'strand': strand,
                        })

    # Merge UTR segments per transcript (take full span)
    meta = {}
    for fbtr, segments in utr_segments.items():
        if not segments:
            continue

        # All segments should be on same chrom/strand
        chrom = segments[0]['chrom']
        strand = segments[0]['strand']

        # Get full span
        all_starts = [s['start'] for s in segments]
        all_ends = [s['end'] for s in segments]

        meta[fbtr] = {
            'fbgn': fbtr_to_fbgn.get(fbtr, fbtr),
            'fbtr': fbtr,
            'chrom': chrom,
            'start': min(all_starts),
            'end': max(all_ends),
            'strand': strand,
        }

    print(f"    Found {len(meta):,} transcripts with {utr_type}")

    return meta


def query_to_genomic_coords(qstart, qend, seq_meta):
    """Convert query-relative coordinates to genomic coordinates."""
    if seq_meta['strand'] == '+':
        genomic_start = seq_meta['start'] + qstart - 1
        genomic_end = seq_meta['start'] + qend - 1
    else:
        # Minus strand: query coords are relative to reverse complement
        genomic_end = seq_meta['end'] - qstart + 1
        genomic_start = seq_meta['end'] - qend + 1
    return genomic_start, genomic_end


def normalize_qseqid(qseqid):
    """Strip shuffle suffix from qseqid if present (e.g., FBtr0070000_shuf1 -> FBtr0070000)."""
    if '_shuf' in qseqid:
        return qseqid.split('_shuf')[0]
    return qseqid


def deduplicate_hits(blast_path, seq_meta, output_hits_file=None):
    """
    Load and deduplicate BLAST hits.

    Args:
        blast_path: Path to BLAST TSV file
        seq_meta: Dictionary mapping qseqid to metadata
        output_hits_file: If provided, write deduplicated hits to this file

    Returns:
        - unique_hits: count of deduplicated hits
        - stats: deduplication statistics
        - gene_stats: per-gene statistics
    """
    seen_keys = {}  # key -> first hit line (for output)
    duplicate_count = 0
    total_count = 0
    skipped_no_meta = 0

    # Track per-gene stats
    gene_stats = defaultdict(lambda: {'total': 0, 'unique': 0, 'duplicates': 0})

    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            total_count += 1

            qseqid = parts[0]
            # Handle shuffled sequences (strip _shufN suffix for metadata lookup)
            lookup_id = normalize_qseqid(qseqid)

            if lookup_id not in seq_meta:
                skipped_no_meta += 1
                continue

            meta = seq_meta[lookup_id]
            fbgn = meta['fbgn']
            chrom = meta['chrom']

            # Parse hit info
            te_id = parts[1]
            qstart = int(parts[6])
            qend = int(parts[7])
            sstart = int(parts[8])
            send = int(parts[9])

            # Convert to genomic coordinates
            genomic_start, genomic_end = query_to_genomic_coords(qstart, qend, meta)

            # Normalize TE coordinates (always start < end)
            te_start = min(sstart, send)
            te_end = max(sstart, send)

            # Create deduplication key
            dedup_key = (fbgn, chrom, genomic_start, genomic_end, te_id, te_start, te_end)

            gene_stats[fbgn]['total'] += 1

            if dedup_key in seen_keys:
                duplicate_count += 1
                gene_stats[fbgn]['duplicates'] += 1
                continue

            # Store first occurrence (keep original line for output)
            seen_keys[dedup_key] = line.strip()
            gene_stats[fbgn]['unique'] += 1

    unique_count = len(seen_keys)

    # Write deduplicated hits if requested
    if output_hits_file:
        with open(output_hits_file, 'w') as f:
            for hit_line in seen_keys.values():
                f.write(hit_line + '\n')

    stats = {
        'total_raw_hits': total_count,
        'unique_hits': unique_count,
        'duplicates_removed': duplicate_count,
        'duplication_rate': 100 * duplicate_count / total_count if total_count > 0 else 0,
        'skipped_no_meta': skipped_no_meta,
        'genes_with_hits': len(gene_stats),
    }

    return unique_count, stats, dict(gene_stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seq-type', required=True, choices=['exon', '3utr', '5utr'],
                        help='Sequence type to deduplicate')
    parser.add_argument('--blast-file', type=Path,
                        help='Override default BLAST file path')
    parser.add_argument('--output-dir', type=Path,
                        help='Override default output directory')
    parser.add_argument('--output-hits', action='store_true',
                        help='Output full deduplicated hit file (same format as input)')
    parser.add_argument('--label', type=str, default='',
                        help='Label for output files (e.g., "real" or "shuffled_rep1")')
    args = parser.parse_args()

    # Set default paths based on sequence type
    if args.seq_type == 'exon':
        default_blast = get_results_dir() / 'exon_analysis' / 'genome_wide_exons.tsv'
        default_output = get_results_dir() / 'exon_analysis' / 'deduplicated'
        metadata_file = Path('data/queries/genome_wide/exon_metadata.tsv')
    elif args.seq_type == '3utr':
        default_blast = get_results_dir() / 'genome_wide_all_3utrs.tsv'
        default_output = get_results_dir() / '3utr_deduplicated'
        metadata_file = None  # Will use GFF
    else:  # 5utr
        default_blast = get_results_dir() / 'genome_wide_all_5utrs.tsv'
        default_output = get_results_dir() / '5utr_deduplicated'
        metadata_file = None  # Will use GFF

    blast_file = args.blast_file or default_blast
    output_dir = args.output_dir or default_output
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"TE Hit Deduplication: {args.seq_type.upper()}")
    print("=" * 70)
    print(f"\nDeduplication key: (fbgn, chrom, genomic_start, genomic_end, te_id, te_start, te_end)")

    # Load metadata
    print(f"\nLoading metadata...")
    if args.seq_type == 'exon':
        seq_meta = load_exon_metadata(metadata_file)
        print(f"  Loaded {len(seq_meta):,} exon records")
    else:
        gff_path = get_references_dir() / 'dmel-all-r6.66.gff'
        utr_type = 'three_prime_UTR' if args.seq_type == '3utr' else 'five_prime_UTR'
        seq_meta = load_utr_metadata_from_gff(gff_path, utr_type)

    # Determine output file names
    label = args.label or args.seq_type
    hits_output = output_dir / f'{label}_deduplicated_hits.tsv' if args.output_hits else None

    # Deduplicate
    print(f"\nProcessing BLAST hits: {blast_file}")
    unique_count, stats, gene_stats = deduplicate_hits(blast_file, seq_meta, hits_output)

    # Print results
    print(f"\n{'=' * 70}")
    print("DEDUPLICATION RESULTS")
    print("=" * 70)
    print(f"\n  BEFORE (raw hits):     {stats['total_raw_hits']:>12,}")
    print(f"  AFTER (unique hits):   {stats['unique_hits']:>12,}")
    print(f"  ─────────────────────────────────────")
    print(f"  Duplicates removed:    {stats['duplicates_removed']:>12,}")
    print(f"  Duplication rate:      {stats['duplication_rate']:>11.2f}%")
    print(f"  Genes with hits:       {stats['genes_with_hits']:>12,}")
    if stats['skipped_no_meta'] > 0:
        print(f"  Skipped (no metadata): {stats['skipped_no_meta']:>12,}")

    # Top genes by duplication
    print(f"\n{'=' * 70}")
    print("TOP 20 GENES BY DUPLICATION COUNT")
    print("=" * 70)

    sorted_genes = sorted(gene_stats.items(), key=lambda x: x[1]['duplicates'], reverse=True)

    print(f"\n  {'Gene':<15} {'Raw':>10} {'Unique':>10} {'Removed':>10} {'Rate':>8}")
    print(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    for fbgn, counts in sorted_genes[:20]:
        rate = 100 * counts['duplicates'] / counts['total'] if counts['total'] > 0 else 0
        print(f"  {fbgn:<15} {counts['total']:>10,} {counts['unique']:>10,} "
              f"{counts['duplicates']:>10,} {rate:>7.1f}%")

    # Save stats
    print(f"\n{'=' * 70}")
    print("Saving results...")
    print("=" * 70)

    # Save deduplication stats
    stats_file = output_dir / f'{label}_deduplication_stats.json'
    with open(stats_file, 'w') as f:
        json.dump({
            'seq_type': args.seq_type,
            'label': label,
            'blast_file': str(blast_file),
            'overall': stats,
            'per_gene_sample': {fbgn: counts for fbgn, counts in sorted_genes[:100]},
        }, f, indent=2)
    print(f"  Saved statistics: {stats_file}")

    # Save per-gene comparison
    compare_file = output_dir / f'{label}_original_vs_deduplicated.tsv'
    with open(compare_file, 'w') as f:
        f.write('fbgn\traw_hits\tunique_hits\tduplicates_removed\tduplication_rate\n')
        for fbgn, counts in sorted(gene_stats.items(), key=lambda x: -x[1]['total']):
            rate = 100 * counts['duplicates'] / counts['total'] if counts['total'] > 0 else 0
            f.write(f"{fbgn}\t{counts['total']}\t{counts['unique']}\t"
                    f"{counts['duplicates']}\t{rate:.2f}\n")
    print(f"  Saved comparison: {compare_file}")

    # Report on hits file if generated
    if hits_output:
        print(f"  Saved deduplicated hits: {hits_output}")

    print(f"\n{'=' * 70}")
    print("Deduplication complete!")
    print("=" * 70)

    return stats


if __name__ == '__main__':
    main()
