#!/usr/bin/env python3
"""
Deduplicate exon TE BLAST hits by genomic coordinates.

When the same genomic region from different transcript isoforms produces
identical BLAST hits (same TE, same TE coordinates), count it only once.

Deduplication key: (fbgn, chrom, genomic_start, genomic_end, te_id, sstart, send)

This ensures:
- Same genomic position hitting same TE region = 1 hit
- Same genomic position hitting different TE regions = multiple hits
- Different genomic positions = separate hits (no merging)
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir
from utils.te_domain_classifier import classify_te_domain


def load_exon_metadata(path):
    """Load exon metadata with genomic coordinates."""
    exons = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            exons[row['exon_id']] = {
                'fbgn': row['fbgn'],
                'gene_symbol': row['gene_symbol'],
                'chrom': row['chrom'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'length': int(row['length']),
                'position': row['position'],
                'utr_overlap': row['utr_overlap'],
            }
    return exons


def query_to_genomic_coords(qstart, qend, exon_meta):
    """Convert query-relative coordinates to genomic coordinates."""
    if exon_meta['strand'] == '+':
        genomic_start = exon_meta['start'] + qstart - 1
        genomic_end = exon_meta['start'] + qend - 1
    else:
        # Minus strand: query coords are relative to reverse complement
        genomic_end = exon_meta['end'] - qstart + 1
        genomic_start = exon_meta['end'] - qend + 1
    return genomic_start, genomic_end


def deduplicate_hits(blast_path, exon_meta):
    """
    Load and deduplicate BLAST hits.

    Returns:
        - unique_hits: list of deduplicated hit dicts
        - stats: deduplication statistics
    """
    seen_keys = {}  # dedup_key -> first hit
    duplicate_count = 0
    total_count = 0
    skipped_no_exon = 0

    # Track duplicates by gene for reporting
    gene_dup_counts = defaultdict(lambda: {'total': 0, 'unique': 0, 'duplicates': 0})

    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 16:
                continue

            total_count += 1

            exon_id = parts[0]
            if exon_id not in exon_meta:
                skipped_no_exon += 1
                continue

            meta = exon_meta[exon_id]
            fbgn = meta['fbgn']
            chrom = meta['chrom']

            # Parse hit info
            te_id = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            mismatch = int(parts[4])
            gapopen = int(parts[5])
            qstart = int(parts[6])
            qend = int(parts[7])
            sstart = int(parts[8])
            send = int(parts[9])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qlen = int(parts[12])
            slen = int(parts[13])
            qseq = parts[14] if len(parts) > 14 else ''
            sseq = parts[15] if len(parts) > 15 else ''

            # Convert to genomic coordinates
            genomic_start, genomic_end = query_to_genomic_coords(qstart, qend, meta)

            # Normalize TE coordinates (always start < end)
            te_start = min(sstart, send)
            te_end = max(sstart, send)

            # Create deduplication key
            dedup_key = (fbgn, chrom, genomic_start, genomic_end, te_id, te_start, te_end)

            gene_dup_counts[fbgn]['total'] += 1

            if dedup_key in seen_keys:
                duplicate_count += 1
                gene_dup_counts[fbgn]['duplicates'] += 1
                continue

            # Classify TE domain
            domain_info = classify_te_domain(te_id, sstart, send, slen)

            # Store unique hit
            hit = {
                'exon_id': exon_id,
                'fbgn': fbgn,
                'gene_symbol': meta['gene_symbol'],
                'chrom': chrom,
                'genomic_start': genomic_start,
                'genomic_end': genomic_end,
                'exon_strand': meta['strand'],
                'position': meta['position'],
                'utr_overlap': meta['utr_overlap'],
                'te_id': te_id,
                'te_start': te_start,
                'te_end': te_end,
                'te_strand': '+' if sstart < send else '-',
                'pident': pident,
                'length': length,
                'evalue': evalue,
                'bitscore': bitscore,
                'te_class': domain_info['te_class'],
                'te_domain': domain_info['domain'],
                'qseq': qseq,
                'sseq': sseq,
            }

            seen_keys[dedup_key] = hit
            gene_dup_counts[fbgn]['unique'] += 1

    unique_hits = list(seen_keys.values())

    stats = {
        'total_raw_hits': total_count,
        'unique_hits': len(unique_hits),
        'duplicates_removed': duplicate_count,
        'duplication_rate': 100 * duplicate_count / total_count if total_count > 0 else 0,
        'skipped_no_exon_meta': skipped_no_exon,
        'genes_with_hits': len(gene_dup_counts),
    }

    return unique_hits, stats, dict(gene_dup_counts)


def aggregate_by_gene(unique_hits):
    """Aggregate deduplicated hits by gene."""
    gene_stats = defaultdict(lambda: {
        'symbol': '',
        'hits': 0,
        'hit_bp': 0,
        'sense_hits': 0,
        'antisense_hits': 0,
        'te_families': set(),
        'domains': defaultdict(int),
    })

    for hit in unique_hits:
        fbgn = hit['fbgn']
        gene_stats[fbgn]['symbol'] = hit['gene_symbol']
        gene_stats[fbgn]['hits'] += 1
        gene_stats[fbgn]['hit_bp'] += hit['length']
        gene_stats[fbgn]['te_families'].add(hit['te_id'])
        gene_stats[fbgn]['domains'][hit['te_domain']] += 1

        if hit['te_strand'] == '+':
            gene_stats[fbgn]['sense_hits'] += 1
        else:
            gene_stats[fbgn]['antisense_hits'] += 1

    return dict(gene_stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--blast-file', type=Path,
                        default=get_results_dir() / 'exon_analysis' / 'genome_wide_exons.tsv')
    parser.add_argument('--exon-metadata', type=Path,
                        default=Path('data/queries/genome_wide/exon_metadata.tsv'))
    parser.add_argument('--output-dir', type=Path,
                        default=get_results_dir() / 'exon_analysis' / 'deduplicated')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Exon TE Hit Deduplication")
    print("=" * 70)
    print("\nDeduplication key: (fbgn, chrom, genomic_start, genomic_end, te_id, te_start, te_end)")
    print("Same genomic region + same TE region = 1 hit (regardless of transcript)")

    # Load exon metadata
    print(f"\nLoading exon metadata: {args.exon_metadata}")
    exon_meta = load_exon_metadata(args.exon_metadata)
    print(f"  Loaded {len(exon_meta):,} exons")

    # Deduplicate hits
    print(f"\nProcessing BLAST hits: {args.blast_file}")
    unique_hits, stats, gene_dup_counts = deduplicate_hits(args.blast_file, exon_meta)

    print(f"\n{'=' * 70}")
    print("DEDUPLICATION RESULTS")
    print("=" * 70)
    print(f"\n  Raw hits:              {stats['total_raw_hits']:>12,}")
    print(f"  Unique hits:           {stats['unique_hits']:>12,}")
    print(f"  Duplicates removed:    {stats['duplicates_removed']:>12,}")
    print(f"  Duplication rate:      {stats['duplication_rate']:>11.2f}%")
    print(f"  Genes with hits:       {stats['genes_with_hits']:>12,}")

    # Aggregate by gene
    print("\nAggregating by gene...")
    gene_stats = aggregate_by_gene(unique_hits)

    # Compare with original counts
    print(f"\n{'=' * 70}")
    print("GENE-LEVEL COMPARISON (Top 20 by duplication)")
    print("=" * 70)

    # Sort genes by duplication count
    sorted_genes = sorted(gene_dup_counts.items(),
                          key=lambda x: x[1]['duplicates'], reverse=True)

    print(f"\n  {'Gene':<15} {'Symbol':<12} {'Raw':>10} {'Unique':>10} {'Removed':>10} {'Rate':>8}")
    print(f"  {'-'*15} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*8}")

    for fbgn, counts in sorted_genes[:20]:
        symbol = gene_stats.get(fbgn, {}).get('symbol', fbgn)[:12]
        rate = 100 * counts['duplicates'] / counts['total'] if counts['total'] > 0 else 0
        print(f"  {fbgn:<15} {symbol:<12} {counts['total']:>10,} {counts['unique']:>10,} "
              f"{counts['duplicates']:>10,} {rate:>7.1f}%")

    # Summary statistics
    dup_rates = [100 * g['duplicates'] / g['total'] for g in gene_dup_counts.values() if g['total'] > 0]
    if dup_rates:
        avg_dup_rate = sum(dup_rates) / len(dup_rates)
        max_dup_rate = max(dup_rates)
        genes_with_dups = sum(1 for r in dup_rates if r > 0)

        print(f"\n  Summary:")
        print(f"    Average gene duplication rate: {avg_dup_rate:.1f}%")
        print(f"    Max gene duplication rate: {max_dup_rate:.1f}%")
        print(f"    Genes with any duplicates: {genes_with_dups:,} / {len(dup_rates):,} ({100*genes_with_dups/len(dup_rates):.1f}%)")

    # Save deduplicated hits
    print(f"\n{'=' * 70}")
    print("Saving results...")
    print("=" * 70)

    # Save full deduplicated hits
    hits_file = args.output_dir / 'exon_te_hits_deduplicated.tsv'
    with open(hits_file, 'w') as f:
        f.write('fbgn\tgene_symbol\tchrom\tgenomic_start\tgenomic_end\texon_id\t'
                'position\tutr_overlap\tte_id\tte_start\tte_end\tte_strand\t'
                'pident\tlength\tevalue\tbitscore\tte_class\tte_domain\n')
        for hit in sorted(unique_hits, key=lambda x: (x['fbgn'], x['genomic_start'])):
            f.write(f"{hit['fbgn']}\t{hit['gene_symbol']}\t{hit['chrom']}\t"
                    f"{hit['genomic_start']}\t{hit['genomic_end']}\t{hit['exon_id']}\t"
                    f"{hit['position']}\t{hit['utr_overlap']}\t{hit['te_id']}\t"
                    f"{hit['te_start']}\t{hit['te_end']}\t{hit['te_strand']}\t"
                    f"{hit['pident']:.1f}\t{hit['length']}\t{hit['evalue']:.2e}\t"
                    f"{hit['bitscore']:.1f}\t{hit['te_class']}\t{hit['te_domain']}\n")
    print(f"  Saved deduplicated hits: {hits_file}")

    # Save gene summary
    gene_file = args.output_dir / 'gene_te_summary_deduplicated.tsv'
    with open(gene_file, 'w') as f:
        f.write('rank\tfbgn\tgene_symbol\tunique_hits\thit_bp\tsense_hits\tantisense_hits\t'
                'sense_pct\tn_te_families\ttop_domain\n')

        sorted_gene_stats = sorted(gene_stats.items(), key=lambda x: -x[1]['hits'])
        for rank, (fbgn, gs) in enumerate(sorted_gene_stats, 1):
            sense_pct = 100 * gs['sense_hits'] / gs['hits'] if gs['hits'] > 0 else 0
            top_domain = max(gs['domains'], key=gs['domains'].get) if gs['domains'] else 'unknown'
            f.write(f"{rank}\t{fbgn}\t{gs['symbol']}\t{gs['hits']}\t{gs['hit_bp']}\t"
                    f"{gs['sense_hits']}\t{gs['antisense_hits']}\t{sense_pct:.1f}\t"
                    f"{len(gs['te_families'])}\t{top_domain}\n")
    print(f"  Saved gene summary: {gene_file}")

    # Save deduplication stats
    stats_file = args.output_dir / 'deduplication_stats.json'
    with open(stats_file, 'w') as f:
        json.dump({
            'overall': stats,
            'per_gene_sample': {fbgn: counts for fbgn, counts in sorted_genes[:100]},
        }, f, indent=2)
    print(f"  Saved statistics: {stats_file}")

    # Save comparison with original
    compare_file = args.output_dir / 'original_vs_deduplicated.tsv'
    with open(compare_file, 'w') as f:
        f.write('fbgn\tgene_symbol\traw_hits\tunique_hits\tduplicates_removed\tduplication_rate\n')
        for fbgn, counts in sorted(gene_dup_counts.items(), key=lambda x: -x[1]['total']):
            symbol = gene_stats.get(fbgn, {}).get('symbol', '')
            rate = 100 * counts['duplicates'] / counts['total'] if counts['total'] > 0 else 0
            f.write(f"{fbgn}\t{symbol}\t{counts['total']}\t{counts['unique']}\t"
                    f"{counts['duplicates']}\t{rate:.2f}\n")
    print(f"  Saved comparison: {compare_file}")

    print(f"\n{'=' * 70}")
    print("Deduplication complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
