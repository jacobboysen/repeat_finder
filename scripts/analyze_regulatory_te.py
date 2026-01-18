#!/usr/bin/env python3
"""
Analyze TE content in regulatory regions with RepeatMasker overlap annotation.

This script processes BLAST results from promoters, enhancers, silencers, and
TFBS against the TE database, annotating which hits overlap known RepeatMasker
annotations.

For each region type:
1. Load BLAST hits and region metadata
2. Convert query-relative to genomic coordinates
3. Check overlap with RepeatMasker annotations
4. Aggregate by region and gene
5. Generate summary statistics

Outputs per region type:
- {region}_te_hits_annotated.tsv: Hits with RepeatMasker overlap flag
- {region}_te_summary.tsv: Per-region summary
- gene_{region}_te_summary.tsv: Per-gene summary
- {region}_analysis_stats.json: Overall statistics
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from bisect import bisect_left

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_references_dir, get_results_dir
from utils.blast_io import BLAST_COLUMNS


# BLAST columns for 16-column output (no strand column in output)
BLAST_COLS_16 = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'slen', 'qseq', 'sseq'
]


def parse_repeatmasker_out(rm_file):
    """Parse RepeatMasker .out file into indexed structure."""
    repeats_by_chrom = defaultdict(list)

    with open(rm_file) as f:
        # Skip header lines
        for _ in range(3):
            next(f)

        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 14:
                continue

            try:
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                repeat_name = parts[9]
                repeat_class = parts[10]

                # Normalize chromosome names
                if chrom.startswith('chr'):
                    chrom = chrom[3:]

                repeats_by_chrom[chrom].append({
                    'start': start,
                    'end': end,
                    'name': repeat_name,
                    'class': repeat_class
                })
            except (ValueError, IndexError):
                continue

    # Sort by start position
    for chrom in repeats_by_chrom:
        repeats_by_chrom[chrom].sort(key=lambda x: x['start'])

    return dict(repeats_by_chrom)


def build_rm_starts_index(rm_index):
    """Build sorted start position lists for binary search."""
    starts_by_chrom = {}
    for chrom, intervals in rm_index.items():
        starts_by_chrom[chrom] = [iv['start'] for iv in intervals]
    return starts_by_chrom


def check_rm_overlap(chrom, start, end, rm_index, rm_starts):
    """
    Check if a region overlaps any RepeatMasker annotation.

    Returns (overlaps: bool, rm_info: dict or None)
    """
    if chrom not in rm_index:
        return False, None

    intervals = rm_index[chrom]
    starts = rm_starts[chrom]

    # Binary search for potential overlaps
    idx = bisect_left(starts, start)

    # Check nearby intervals
    for i in range(max(0, idx - 5), min(len(intervals), idx + 10)):
        iv = intervals[i]
        if iv['end'] < start:
            continue
        if iv['start'] > end:
            break
        # Overlap found
        return True, iv

    return False, None


def load_metadata(metadata_path, region_type):
    """Load region metadata with appropriate column names."""
    df = pd.read_csv(metadata_path, sep='\t')
    return df


def classify_strand(sstart, send):
    """Classify strand based on subject coordinates."""
    return 'minus' if sstart > send else 'plus'


def analyze_region_type(
    blast_path,
    metadata_path,
    rm_index,
    rm_starts,
    output_dir,
    region_type,
    verbose=False
):
    """
    Analyze TE hits for a single region type.

    Returns statistics dictionary.
    """
    if verbose:
        print(f"\nAnalyzing {region_type}...")

    # Load metadata
    metadata_df = load_metadata(metadata_path, region_type)

    # Build metadata lookup
    if region_type == 'promoter':
        id_col = 'promoter_id'
        chrom_col = 'chrom'
        start_col = 'region_start'
        end_col = 'region_end'
        strand_col = 'strand'
        gene_col = 'fbgn'
        symbol_col = 'gene_symbol'
    elif region_type == 'enhancer':
        id_col = 'enhancer_id'
        chrom_col = 'chrom'
        start_col = 'start'
        end_col = 'end'
        strand_col = None  # Enhancers don't have strand
        gene_col = 'nearest_fbgn'
        symbol_col = 'nearest_symbol'
    elif region_type == 'silencer':
        id_col = 'silencer_id'
        chrom_col = 'chrom'
        start_col = 'start'
        end_col = 'end'
        strand_col = None
        gene_col = 'nearest_fbgn'
        symbol_col = 'nearest_symbol'
    elif region_type == 'tfbs':
        id_col = 'tfbs_id'
        chrom_col = 'chrom'
        start_col = 'start'
        end_col = 'end'
        strand_col = None
        gene_col = 'nearest_fbgn'
        symbol_col = 'nearest_symbol'
    else:
        raise ValueError(f"Unknown region type: {region_type}")

    metadata_lookup = {}
    for _, row in metadata_df.iterrows():
        region_id = row[id_col]
        metadata_lookup[region_id] = {
            'chrom': row[chrom_col],
            'start': int(row[start_col]),
            'end': int(row[end_col]),
            'strand': row[strand_col] if strand_col and strand_col in row else '+',
            'fbgn': row[gene_col] if pd.notna(row[gene_col]) else None,
            'symbol': row[symbol_col] if pd.notna(row.get(symbol_col, None)) else None,
            'length': int(row['length'])
        }

    # Load BLAST results
    if verbose:
        print(f"  Loading BLAST results from {blast_path}")

    blast_df = pd.read_csv(blast_path, sep='\t', names=BLAST_COLS_16)

    if verbose:
        print(f"  Total hits: {len(blast_df):,}")

    # Process hits
    annotated_hits = []
    stats = {
        'total_hits': len(blast_df),
        'hits_in_rm': 0,
        'hits_novel': 0,
        'regions_with_hits': set(),
        'genes_with_hits': set(),
        'te_families': defaultdict(int),
        'quality_bins': {
            'high': 0,  # >=80% identity, >=50bp
            'medium': 0,  # >=70% identity, >=30bp
            'low': 0  # rest
        }
    }

    for _, hit in blast_df.iterrows():
        qseqid = hit['qseqid']

        # Skip if region not in metadata (e.g., filtered out)
        if qseqid not in metadata_lookup:
            continue

        region = metadata_lookup[qseqid]
        chrom = region['chrom']
        region_start = region['start']
        region_strand = region['strand']

        # Convert query coordinates to genomic coordinates
        qstart, qend = int(hit['qstart']), int(hit['qend'])

        if region_strand == '+':
            genomic_start = region_start + qstart - 1
            genomic_end = region_start + qend - 1
        else:
            # Minus strand: coordinates are reversed
            genomic_end = region['end'] - qstart + 1
            genomic_start = region['end'] - qend + 1

        # Check RepeatMasker overlap
        overlaps_rm, rm_info = check_rm_overlap(
            chrom, genomic_start, genomic_end, rm_index, rm_starts
        )

        if overlaps_rm:
            stats['hits_in_rm'] += 1
        else:
            stats['hits_novel'] += 1

        # Track regions and genes
        stats['regions_with_hits'].add(qseqid)
        if region['fbgn']:
            stats['genes_with_hits'].add(region['fbgn'])

        # Track TE family
        te_id = hit['sseqid']
        te_family = te_id.split('\\')[0] if '\\' in te_id else te_id.split('_')[0]
        stats['te_families'][te_family] += 1

        # Quality bin
        pident = float(hit['pident'])
        length = int(hit['length'])
        if pident >= 80 and length >= 50:
            stats['quality_bins']['high'] += 1
        elif pident >= 70 and length >= 30:
            stats['quality_bins']['medium'] += 1
        else:
            stats['quality_bins']['low'] += 1

        # Classify strand
        sstart, send = int(hit['sstart']), int(hit['send'])
        te_strand = classify_strand(sstart, send)

        # Build annotated hit record
        annotated_hits.append({
            'region_id': qseqid,
            'chrom': chrom,
            'genomic_start': genomic_start,
            'genomic_end': genomic_end,
            'query_start': qstart,
            'query_end': qend,
            'te_id': hit['sseqid'],
            'te_start': min(sstart, send),
            'te_end': max(sstart, send),
            'te_strand': te_strand,
            'pident': pident,
            'length': length,
            'evalue': hit['evalue'],
            'bitscore': hit['bitscore'],
            'in_repeatmasker': overlaps_rm,
            'rm_name': rm_info['name'] if rm_info else '',
            'rm_class': rm_info['class'] if rm_info else '',
            'fbgn': region['fbgn'] or '',
            'symbol': region['symbol'] or ''
        })

    # Convert sets to counts
    stats['regions_with_hits'] = len(stats['regions_with_hits'])
    stats['genes_with_hits'] = len(stats['genes_with_hits'])
    stats['te_families'] = dict(stats['te_families'])

    if verbose:
        print(f"  Hits in RepeatMasker: {stats['hits_in_rm']:,} ({100*stats['hits_in_rm']/stats['total_hits']:.1f}%)")
        print(f"  Novel hits: {stats['hits_novel']:,} ({100*stats['hits_novel']/stats['total_hits']:.1f}%)")
        print(f"  Regions with hits: {stats['regions_with_hits']:,}")
        print(f"  Genes with hits: {stats['genes_with_hits']:,}")

    # Write annotated hits
    output_dir.mkdir(parents=True, exist_ok=True)

    hits_df = pd.DataFrame(annotated_hits)
    hits_path = output_dir / f'{region_type}_te_hits_annotated.tsv'
    hits_df.to_csv(hits_path, sep='\t', index=False)
    if verbose:
        print(f"  Wrote annotated hits to: {hits_path}")

    # Create per-region summary
    region_summary = hits_df.groupby('region_id').agg({
        'length': ['count', 'sum'],
        'pident': 'mean',
        'in_repeatmasker': 'sum',
        'fbgn': 'first',
        'symbol': 'first'
    }).reset_index()
    region_summary.columns = ['region_id', 'n_hits', 'total_hit_bp', 'mean_pident',
                              'n_in_rm', 'fbgn', 'symbol']

    # Add region length for density calculation
    region_summary['region_length'] = region_summary['region_id'].map(
        lambda x: metadata_lookup.get(x, {}).get('length', 0)
    )
    region_summary['density'] = (region_summary['total_hit_bp'] / region_summary['region_length']) * 1000
    region_summary['pct_in_rm'] = (region_summary['n_in_rm'] / region_summary['n_hits']) * 100

    summary_path = output_dir / f'{region_type}_te_summary.tsv'
    region_summary.to_csv(summary_path, sep='\t', index=False)
    if verbose:
        print(f"  Wrote region summary to: {summary_path}")

    # Create per-gene summary
    gene_summary = hits_df[hits_df['fbgn'] != ''].groupby('fbgn').agg({
        'region_id': 'nunique',
        'length': ['count', 'sum'],
        'pident': 'mean',
        'in_repeatmasker': 'sum',
        'symbol': 'first'
    }).reset_index()
    gene_summary.columns = ['fbgn', 'n_regions', 'n_hits', 'total_hit_bp',
                            'mean_pident', 'n_in_rm', 'symbol']
    gene_summary['pct_in_rm'] = (gene_summary['n_in_rm'] / gene_summary['n_hits']) * 100

    gene_path = output_dir / f'gene_{region_type}_te_summary.tsv'
    gene_summary.to_csv(gene_path, sep='\t', index=False)
    if verbose:
        print(f"  Wrote gene summary to: {gene_path}")

    # Write stats
    stats_path = output_dir / f'{region_type}_analysis_stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)

    return stats


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--results-dir',
        type=Path,
        default=get_results_dir() / 'regulatory_analysis',
        help='Results directory containing BLAST outputs'
    )
    parser.add_argument(
        '--queries-dir',
        type=Path,
        default=Path('data/queries/regulatory'),
        help='Queries directory containing metadata files'
    )
    parser.add_argument(
        '--rm-file',
        type=Path,
        default=get_references_dir() / 'dm6.fa.out',
        help='RepeatMasker output file'
    )
    parser.add_argument(
        '--region-types',
        nargs='+',
        default=['promoter', 'enhancer', 'silencer', 'tfbs'],
        help='Region types to analyze'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    print("Regulatory Region TE Analysis")
    print("=" * 60)

    # Load RepeatMasker annotations
    if args.verbose:
        print(f"\nLoading RepeatMasker annotations from: {args.rm_file}")

    rm_index = parse_repeatmasker_out(args.rm_file)
    rm_starts = build_rm_starts_index(rm_index)

    total_rm = sum(len(v) for v in rm_index.values())
    if args.verbose:
        print(f"  Loaded {total_rm:,} RepeatMasker annotations")

    # Analyze each region type
    all_stats = {}

    for region_type in args.region_types:
        blast_path = args.results_dir / region_type / f'{region_type}_te_hits.tsv'
        # Handle plural directory names (promoters, enhancers, silencers, tfbs)
        dir_name = region_type if region_type == 'tfbs' else f'{region_type}s'
        metadata_path = args.queries_dir / dir_name / f'{region_type}_metadata.tsv'

        if not blast_path.exists():
            print(f"\nSkipping {region_type}: BLAST results not found at {blast_path}")
            continue

        if not metadata_path.exists():
            print(f"\nSkipping {region_type}: Metadata not found at {metadata_path}")
            continue

        output_dir = args.results_dir / region_type

        stats = analyze_region_type(
            blast_path,
            metadata_path,
            rm_index,
            rm_starts,
            output_dir,
            region_type,
            args.verbose
        )

        all_stats[region_type] = stats

    # Write integrated summary
    print("\n" + "=" * 60)
    print("Summary")
    print("-" * 60)

    summary_rows = []
    for region_type, stats in all_stats.items():
        pct_novel = 100 * stats['hits_novel'] / stats['total_hits'] if stats['total_hits'] > 0 else 0
        summary_rows.append({
            'region_type': region_type,
            'total_hits': stats['total_hits'],
            'hits_in_rm': stats['hits_in_rm'],
            'hits_novel': stats['hits_novel'],
            'pct_novel': pct_novel,
            'regions_with_hits': stats['regions_with_hits'],
            'genes_with_hits': stats['genes_with_hits'],
            'high_quality_hits': stats['quality_bins']['high']
        })

        print(f"\n{region_type.upper()}:")
        print(f"  Total hits: {stats['total_hits']:,}")
        print(f"  In RepeatMasker: {stats['hits_in_rm']:,} ({100-pct_novel:.1f}%)")
        print(f"  Novel: {stats['hits_novel']:,} ({pct_novel:.1f}%)")
        print(f"  Regions with hits: {stats['regions_with_hits']:,}")
        print(f"  Genes with hits: {stats['genes_with_hits']:,}")

    # Write summary table
    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.results_dir / 'integrated_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False)
    print(f"\nWrote integrated summary to: {summary_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
