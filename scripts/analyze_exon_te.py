#!/usr/bin/env python3
"""
Analyze TE similarity in individual exons.

Complements the 5'UTR and 3'UTR analyses by examining coding exons for TE-like
sequences. Key analyses:
- Compare TE density: internal exons vs UTR-overlapping exons
- TE domain distribution: are exon hits enriched in gag/pol/env regions?
- Strand bias: sense vs antisense TE orientation in exons
- Comparison with UTR results

Outputs:
- exon_te_summary.tsv: Per-exon TE statistics
- gene_exon_te_summary.tsv: Per-gene aggregation
- te_domain_distribution.tsv: Hits by TE domain type
- utr_overlap_comparison.tsv: Terminal vs internal exon comparison
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).parent))

from utils.blast_io import load_blast_results, classify_strand, BLAST_COLUMNS
from utils.data_loaders import load_te_database, parse_te_class
from utils.te_domain_classifier import classify_te_domain, infer_te_class, is_coding_domain


def load_exon_metadata(metadata_path: Path) -> Dict[str, dict]:
    """Load exon metadata from TSV file."""
    metadata = {}
    with open(metadata_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            metadata[row['exon_id']] = {
                'fbtr': row['fbtr'],
                'fbgn': row['fbgn'],
                'gene_symbol': row['gene_symbol'],
                'exon_num': int(row['exon_num']),
                'total_exons': int(row['total_exons']),
                'position': row['position'],
                'utr_overlap': row['utr_overlap'],
                'chrom': row['chrom'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': row['strand'],
                'length': int(row['length'])
            }
    return metadata


def analyze_blast_results(
    blast_path: Path,
    metadata_path: Path,
    te_db_path: Optional[Path] = None,
    verbose: bool = False
) -> Tuple[dict, dict, dict, dict]:
    """
    Analyze BLAST results for exon TE content.

    Returns:
        exon_stats: Per-exon statistics
        gene_stats: Per-gene aggregated statistics
        domain_stats: TE domain distribution
        position_comparison: Terminal vs internal exon comparison
    """
    # Load exon metadata
    if verbose:
        print(f"Loading exon metadata from: {metadata_path}")
    exon_metadata = load_exon_metadata(metadata_path)
    if verbose:
        print(f"  Loaded {len(exon_metadata):,} exons")

    # Load TE database for domain classification
    te_info = {}
    if te_db_path and te_db_path.exists():
        if verbose:
            print(f"Loading TE database from: {te_db_path}")
        te_info = load_te_database(te_db_path)
        if verbose:
            print(f"  Loaded {len(te_info):,} TEs")

    # Initialize counters
    exon_stats = defaultdict(lambda: {
        'hits': 0,
        'hit_bp': 0,
        'sense_hits': 0,
        'antisense_hits': 0,
        'te_families': set(),
        'te_classes': defaultdict(int),
        'te_domains': defaultdict(int),
    })

    domain_counts = defaultdict(int)
    class_counts = defaultdict(int)
    position_hits = defaultdict(lambda: {'total': 0, 'sense': 0, 'antisense': 0, 'exons': set()})

    total_hits = 0
    processed_exons = set()

    # Process BLAST results
    if verbose:
        print(f"\nProcessing BLAST results from: {blast_path}")

    with open(blast_path) as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 500_000 == 0 and verbose:
                print(f"  Processed {line_num:,} hits...")

            parts = line.strip().split('\t')
            if len(parts) < 16:
                continue

            total_hits += 1

            # Parse hit info
            qseqid = parts[0]  # exon ID
            sseqid = parts[1]  # TE ID
            pident = float(parts[2])
            length = int(parts[3])
            sstart = int(parts[8])
            send = int(parts[9])
            slen = int(parts[13])

            # Classify strand
            strand = classify_strand(sstart, send)

            # Get exon metadata
            if qseqid not in exon_metadata:
                continue

            exon_meta = exon_metadata[qseqid]
            processed_exons.add(qseqid)

            # Update exon stats
            stats = exon_stats[qseqid]
            stats['hits'] += 1
            stats['hit_bp'] += length

            if strand == 'plus':
                stats['sense_hits'] += 1
            else:
                stats['antisense_hits'] += 1

            # Get TE info
            te_class = 'Unknown'
            if sseqid in te_info:
                te_class = te_info[sseqid].get('class', 'Unknown')
                stats['te_families'].add(te_info[sseqid].get('name', sseqid))
            else:
                te_class = infer_te_class(sseqid)

            stats['te_classes'][te_class] += 1
            class_counts[te_class] += 1

            # Classify TE domain
            domain_info = classify_te_domain(sseqid, sstart, send, slen, te_class)
            domain = domain_info['domain']
            stats['te_domains'][domain] += 1
            domain_counts[domain] += 1

            # Track by exon position
            position = exon_meta['position']
            position_hits[position]['total'] += 1
            position_hits[position]['exons'].add(qseqid)
            if strand == 'plus':
                position_hits[position]['sense'] += 1
            else:
                position_hits[position]['antisense'] += 1

    if verbose:
        print(f"\n  Total hits: {total_hits:,}")
        print(f"  Exons with hits: {len(processed_exons):,}")

    # Finalize exon stats
    for exon_id, stats in exon_stats.items():
        stats['te_families'] = list(stats['te_families'])
        stats['te_classes'] = dict(stats['te_classes'])
        stats['te_domains'] = dict(stats['te_domains'])

        # Add metadata
        if exon_id in exon_metadata:
            meta = exon_metadata[exon_id]
            stats['length'] = meta['length']
            stats['position'] = meta['position']
            stats['utr_overlap'] = meta['utr_overlap']
            stats['fbgn'] = meta['fbgn']
            stats['fbtr'] = meta['fbtr']
            stats['gene_symbol'] = meta['gene_symbol']

            # Calculate density
            if meta['length'] > 0:
                stats['density'] = (stats['hit_bp'] / meta['length']) * 100
            else:
                stats['density'] = 0

    # Aggregate by gene
    gene_stats = defaultdict(lambda: {
        'exons': 0,
        'exons_with_hits': 0,
        'total_hits': 0,
        'total_hit_bp': 0,
        'total_exon_bp': 0,
        'sense_hits': 0,
        'antisense_hits': 0,
        'te_families': set(),
    })

    for exon_id, stats in exon_stats.items():
        if 'fbgn' not in stats:
            continue

        fbgn = stats['fbgn']
        gene = gene_stats[fbgn]
        gene['exons'] += 1
        gene['gene_symbol'] = stats.get('gene_symbol', fbgn)

        if stats['hits'] > 0:
            gene['exons_with_hits'] += 1

        gene['total_hits'] += stats['hits']
        gene['total_hit_bp'] += stats['hit_bp']
        gene['total_exon_bp'] += stats.get('length', 0)
        gene['sense_hits'] += stats['sense_hits']
        gene['antisense_hits'] += stats['antisense_hits']
        gene['te_families'].update(stats.get('te_families', []))

    # Finalize gene stats
    for fbgn, stats in gene_stats.items():
        stats['te_families'] = list(stats['te_families'])
        if stats['total_exon_bp'] > 0:
            stats['density'] = (stats['total_hit_bp'] / stats['total_exon_bp']) * 100
        else:
            stats['density'] = 0

        total = stats['total_hits']
        if total > 0:
            stats['sense_pct'] = 100 * stats['sense_hits'] / total
            stats['antisense_pct'] = 100 * stats['antisense_hits'] / total
        else:
            stats['sense_pct'] = 50
            stats['antisense_pct'] = 50

    # Build position comparison
    position_comparison = {}
    for position, data in position_hits.items():
        total = data['total']
        position_comparison[position] = {
            'total_hits': total,
            'exons_with_hits': len(data['exons']),
            'sense_hits': data['sense'],
            'antisense_hits': data['antisense'],
            'sense_pct': 100 * data['sense'] / total if total > 0 else 50,
        }

    return dict(exon_stats), dict(gene_stats), dict(domain_counts), position_comparison


def write_results(
    output_dir: Path,
    exon_stats: dict,
    gene_stats: dict,
    domain_stats: dict,
    position_comparison: dict,
    verbose: bool = False
):
    """Write analysis results to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write per-exon summary
    exon_file = output_dir / 'exon_te_summary.tsv'
    with open(exon_file, 'w') as f:
        header = [
            'exon_id', 'fbtr', 'fbgn', 'gene_symbol', 'position', 'utr_overlap',
            'length', 'hits', 'hit_bp', 'density', 'sense_hits', 'antisense_hits',
            'n_te_families', 'top_te_class'
        ]
        f.write('\t'.join(header) + '\n')

        for exon_id, stats in sorted(exon_stats.items(), key=lambda x: -x[1].get('hits', 0)):
            if stats.get('hits', 0) == 0:
                continue

            # Get top TE class
            te_classes = stats.get('te_classes', {})
            top_class = max(te_classes, key=te_classes.get) if te_classes else 'Unknown'

            row = [
                exon_id,
                stats.get('fbtr', ''),
                stats.get('fbgn', ''),
                stats.get('gene_symbol', ''),
                stats.get('position', ''),
                stats.get('utr_overlap', ''),
                str(stats.get('length', 0)),
                str(stats.get('hits', 0)),
                str(stats.get('hit_bp', 0)),
                f"{stats.get('density', 0):.2f}",
                str(stats.get('sense_hits', 0)),
                str(stats.get('antisense_hits', 0)),
                str(len(stats.get('te_families', []))),
                top_class
            ]
            f.write('\t'.join(row) + '\n')

    if verbose:
        print(f"Wrote per-exon summary to: {exon_file}")

    # Write per-gene summary
    gene_file = output_dir / 'gene_exon_te_summary.tsv'
    with open(gene_file, 'w') as f:
        header = [
            'rank', 'fbgn', 'gene_symbol', 'exons', 'exons_with_hits',
            'total_hits', 'total_hit_bp', 'total_exon_bp', 'density',
            'sense_hits', 'antisense_hits', 'sense_pct', 'n_te_families'
        ]
        f.write('\t'.join(header) + '\n')

        sorted_genes = sorted(gene_stats.items(), key=lambda x: -x[1].get('density', 0))

        for rank, (fbgn, stats) in enumerate(sorted_genes, 1):
            row = [
                str(rank),
                fbgn,
                stats.get('gene_symbol', ''),
                str(stats.get('exons', 0)),
                str(stats.get('exons_with_hits', 0)),
                str(stats.get('total_hits', 0)),
                str(stats.get('total_hit_bp', 0)),
                str(stats.get('total_exon_bp', 0)),
                f"{stats.get('density', 0):.2f}",
                str(stats.get('sense_hits', 0)),
                str(stats.get('antisense_hits', 0)),
                f"{stats.get('sense_pct', 50):.1f}",
                str(len(stats.get('te_families', [])))
            ]
            f.write('\t'.join(row) + '\n')

    if verbose:
        print(f"Wrote per-gene summary to: {gene_file}")

    # Write TE domain distribution
    domain_file = output_dir / 'te_domain_distribution.tsv'
    with open(domain_file, 'w') as f:
        f.write('domain\thits\tpct\tis_coding\n')
        total = sum(domain_stats.values())
        for domain, count in sorted(domain_stats.items(), key=lambda x: -x[1]):
            pct = 100 * count / total if total > 0 else 0
            coding = 'yes' if is_coding_domain(domain) else 'no'
            f.write(f"{domain}\t{count}\t{pct:.2f}\t{coding}\n")

    if verbose:
        print(f"Wrote TE domain distribution to: {domain_file}")

    # Write position comparison
    position_file = output_dir / 'utr_overlap_comparison.tsv'
    with open(position_file, 'w') as f:
        f.write('position\ttotal_hits\texons_with_hits\tsense_hits\tantisense_hits\tsense_pct\n')
        for position, data in sorted(position_comparison.items()):
            f.write(f"{position}\t{data['total_hits']}\t{data['exons_with_hits']}\t"
                    f"{data['sense_hits']}\t{data['antisense_hits']}\t{data['sense_pct']:.1f}\n")

    if verbose:
        print(f"Wrote position comparison to: {position_file}")


def print_summary(
    exon_stats: dict,
    gene_stats: dict,
    domain_stats: dict,
    position_comparison: dict
):
    """Print analysis summary to console."""
    print("\n" + "=" * 70)
    print("EXON TE ANALYSIS SUMMARY")
    print("=" * 70)

    # Overall stats
    exons_with_hits = sum(1 for s in exon_stats.values() if s.get('hits', 0) > 0)
    total_hits = sum(s.get('hits', 0) for s in exon_stats.values())
    genes_with_hits = sum(1 for s in gene_stats.values() if s.get('total_hits', 0) > 0)

    print(f"\nOverall Statistics:")
    print(f"  Total exons analyzed: {len(exon_stats):,}")
    print(f"  Exons with TE hits: {exons_with_hits:,} ({100*exons_with_hits/len(exon_stats):.1f}%)")
    print(f"  Total BLAST hits: {total_hits:,}")
    print(f"  Genes with exon TE hits: {genes_with_hits:,}")

    # Position comparison
    print(f"\nTE Hits by Exon Position:")
    print(f"  {'Position':<15} {'Hits':>10} {'Exons':>10} {'Sense%':>10}")
    print(f"  {'-'*15} {'-'*10} {'-'*10} {'-'*10}")
    for position in ['first_exon', 'internal_exon', 'last_exon', 'single_exon']:
        if position in position_comparison:
            data = position_comparison[position]
            print(f"  {position:<15} {data['total_hits']:>10,} {data['exons_with_hits']:>10,} "
                  f"{data['sense_pct']:>9.1f}%")

    # TE domain distribution
    print(f"\nTE Domain Distribution (top 10):")
    print(f"  {'Domain':<20} {'Hits':>10} {'%':>8} {'Coding':>8}")
    print(f"  {'-'*20} {'-'*10} {'-'*8} {'-'*8}")
    total = sum(domain_stats.values())
    for domain, count in sorted(domain_stats.items(), key=lambda x: -x[1])[:10]:
        pct = 100 * count / total if total > 0 else 0
        coding = 'yes' if is_coding_domain(domain) else 'no'
        print(f"  {domain:<20} {count:>10,} {pct:>7.1f}% {coding:>8}")

    # Top genes by density
    print(f"\nTop 15 Genes by Exon TE Density:")
    print(f"  {'Rank':<5} {'Gene':<12} {'Density':>10} {'Hits':>8} {'Exons':>8}")
    print(f"  {'-'*5} {'-'*12} {'-'*10} {'-'*8} {'-'*8}")

    sorted_genes = sorted(gene_stats.items(), key=lambda x: -x[1].get('density', 0))[:15]
    for rank, (fbgn, stats) in enumerate(sorted_genes, 1):
        symbol = stats.get('gene_symbol', fbgn)[:12]
        print(f"  {rank:<5} {symbol:<12} {stats['density']:>10.1f} "
              f"{stats['total_hits']:>8,} {stats['exons']:>8}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--blast',
        type=Path,
        default=Path('results/exon_analysis/genome_wide_exons.tsv'),
        help='BLAST results file (default: results/exon_analysis/genome_wide_exons.tsv)'
    )
    parser.add_argument(
        '--metadata',
        type=Path,
        default=Path('data/queries/genome_wide/exon_metadata.tsv'),
        help='Exon metadata file (default: data/queries/genome_wide/exon_metadata.tsv)'
    )
    parser.add_argument(
        '--te-db',
        type=Path,
        default=Path('data/references/dmel_te_flybase.fasta'),
        help='TE database FASTA for family/class info'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('results/exon_analysis'),
        help='Output directory (default: results/exon_analysis)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.blast.exists():
        print(f"Error: BLAST results not found: {args.blast}", file=sys.stderr)
        return 1

    if not args.metadata.exists():
        print(f"Error: Metadata file not found: {args.metadata}", file=sys.stderr)
        return 1

    print("Analyzing Exon TE Content")
    print("=" * 60)

    # Run analysis
    exon_stats, gene_stats, domain_stats, position_comparison = analyze_blast_results(
        args.blast,
        args.metadata,
        args.te_db,
        args.verbose
    )

    # Write results
    write_results(
        args.output_dir,
        exon_stats,
        gene_stats,
        domain_stats,
        position_comparison,
        args.verbose
    )

    # Print summary
    print_summary(exon_stats, gene_stats, domain_stats, position_comparison)

    return 0


if __name__ == '__main__':
    sys.exit(main())
