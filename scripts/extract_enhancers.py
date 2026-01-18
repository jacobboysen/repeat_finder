#!/usr/bin/env python3
"""
Extract enhancer/CRM sequences from regulatory_region annotations.

Uses REDfly-curated CRMs (cis-regulatory modules) integrated into FlyBase GFF.
Assigns nearest gene and calculates distance for regions without existing
gene associations.

Sources include:
- REDfly_CRMs: Literature-curated CRMs
- GMR_Brain_exp_1_REDfly_CRMs: Janelia GAL4 lines
- VDRC_VT_REDfly_CRMs: Vienna Tile library
- STARR-seq validated regions

Outputs:
- enhancers_sense.fasta: Enhancer sequences (sense strand)
- enhancer_metadata.tsv: Detailed metadata including nearest gene
- extraction_stats.json: Summary statistics

Header format:
>FBsf0000926668 name=Unspecified_Kc167-AF_1 nearest_gene=l(2)gl:FBgn0002121 dist=500 loc=2L:5468-5967 length=500
"""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from bisect import bisect_left

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_references_dir, get_queries_dir


class SequenceCache:
    """Cache for deduplicating identical sequences."""

    def __init__(self):
        self.seen_hashes = {}
        self.duplicates = defaultdict(list)
        self.duplicates_skipped = 0

    def add(self, record: SeqRecord) -> bool:
        seq_hash = hashlib.md5(str(record.seq).upper().encode()).hexdigest()
        if seq_hash in self.seen_hashes:
            self.duplicates[seq_hash].append(record.id)
            self.duplicates_skipped += 1
            return False
        else:
            self.seen_hashes[seq_hash] = record
            self.duplicates[seq_hash].append(record.id)
            return True

    def get_unique_records(self) -> List[SeqRecord]:
        return list(self.seen_hashes.values())

    def get_stats(self) -> dict:
        return {
            'unique_sequences': len(self.seen_hashes),
            'duplicates_skipped': self.duplicates_skipped,
            'total_regions': sum(len(v) for v in self.duplicates.values())
        }


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF9 attributes column into dictionary."""
    attrs = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs


class GeneIndex:
    """
    Index for fast nearest-gene lookup.

    Stores gene coordinates sorted by midpoint for binary search.
    """

    def __init__(self):
        self.genes_by_chrom = defaultdict(list)  # chrom -> [(midpoint, start, end, fbgn, symbol)]
        self.built = False

    def add_gene(self, chrom: str, start: int, end: int, fbgn: str, symbol: str):
        """Add a gene to the index."""
        midpoint = (start + end) // 2
        self.genes_by_chrom[chrom].append((midpoint, start, end, fbgn, symbol))

    def build(self):
        """Sort genes by midpoint for binary search."""
        for chrom in self.genes_by_chrom:
            self.genes_by_chrom[chrom].sort(key=lambda x: x[0])
        self.built = True

    def find_nearest(self, chrom: str, region_start: int, region_end: int) -> Tuple[Optional[str], Optional[str], int]:
        """
        Find the nearest gene to a region.

        Returns:
            (fbgn, symbol, distance)
            Distance is 0 if overlapping, positive otherwise.
        """
        if not self.built:
            self.build()

        genes = self.genes_by_chrom.get(chrom, [])
        if not genes:
            return None, None, -1

        region_mid = (region_start + region_end) // 2

        # Binary search for closest midpoint
        midpoints = [g[0] for g in genes]
        idx = bisect_left(midpoints, region_mid)

        # Check genes around the insertion point
        candidates = []
        for i in range(max(0, idx - 2), min(len(genes), idx + 3)):
            gene_mid, gene_start, gene_end, fbgn, symbol = genes[i]

            # Calculate distance
            if region_end < gene_start:
                # Region is upstream of gene
                dist = gene_start - region_end
            elif region_start > gene_end:
                # Region is downstream of gene
                dist = region_start - gene_end
            else:
                # Overlapping
                dist = 0

            candidates.append((dist, fbgn, symbol))

        # Return closest
        if candidates:
            candidates.sort(key=lambda x: x[0])
            dist, fbgn, symbol = candidates[0]
            return fbgn, symbol, dist

        return None, None, -1


def load_genes_and_regulatory_regions(gff_path: Path, verbose: bool = False) -> Tuple[GeneIndex, List[dict]]:
    """
    Parse GFF to extract gene coordinates and regulatory_region features.

    Returns:
        gene_index: Index for nearest-gene lookup
        regulatory_regions: List of regulatory region records
    """
    gene_index = GeneIndex()
    regulatory_regions = []

    line_count = 0
    gene_count = 0
    reg_count = 0

    if verbose:
        print(f"Parsing GFF: {gff_path}")

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            line_count += 1
            if verbose and line_count % 1_000_000 == 0:
                print(f"  Processed {line_count:,} lines...")

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            attrs = parse_gff_attributes(attributes)

            if feature_type == 'gene':
                fbgn = attrs.get('ID', '')
                symbol = attrs.get('Name', fbgn)
                if fbgn.startswith('FBgn'):
                    gene_index.add_gene(chrom, start, end, fbgn, symbol)
                    gene_count += 1

            elif feature_type == 'regulatory_region':
                reg_id = attrs.get('ID', '')
                name = attrs.get('Name', '')

                # Parse associated genes if present
                assoc_fbgn = None
                assoc_symbol = None
                if 'associated_genes' in attrs:
                    assoc = attrs['associated_genes'].split(',')[0]
                    if ':' in assoc:
                        assoc_symbol, assoc_fbgn = assoc.split(':', 1)

                # Check if STARR-seq
                is_starr = 'STARR' in name or 'STARR' in source

                regulatory_regions.append({
                    'reg_id': reg_id,
                    'name': name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand if strand in ['+', '-'] else '+',
                    'source': source,
                    'assoc_fbgn': assoc_fbgn,
                    'assoc_symbol': assoc_symbol,
                    'is_starr': is_starr
                })
                reg_count += 1

    gene_index.build()

    if verbose:
        print(f"  Total lines: {line_count:,}")
        print(f"  Genes: {gene_count:,}")
        print(f"  Regulatory regions: {reg_count:,}")

    return gene_index, regulatory_regions


def load_genome(genome_path: Path, verbose: bool = False) -> Dict[str, str]:
    """Load genome sequences into dictionary."""
    if verbose:
        print(f"Loading genome: {genome_path}")

    genome = {}
    for record in SeqIO.parse(genome_path, 'fasta'):
        genome[record.id] = str(record.seq)

    if verbose:
        print(f"  Loaded {len(genome)} chromosomes/scaffolds")

    return genome


def extract_enhancers(
    gff_path: Path,
    genome_path: Path,
    output_dir: Path,
    verbose: bool = False,
    deduplicate: bool = True,
    major_chroms_only: bool = True,
    exclude_repressed: bool = True
) -> dict:
    """
    Extract enhancer/CRM sequences with nearest gene annotation.

    Args:
        gff_path: Path to GFF annotation file
        genome_path: Path to genome FASTA
        output_dir: Output directory
        verbose: Print progress
        deduplicate: Remove duplicate sequences
        major_chroms_only: Only use major chromosomes
        exclude_repressed: Exclude regions with "Repressed" in name (silencers)

    Returns:
        Statistics dictionary
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    major_chroms = {'2L', '2R', '3L', '3R', 'X', '4'}

    # Load annotations
    gene_index, regulatory_regions = load_genes_and_regulatory_regions(gff_path, verbose)

    # Load genome
    genome = load_genome(genome_path, verbose)

    if verbose:
        print(f"\nExtracting enhancer sequences...")

    all_records = []
    metadata_rows = []
    cache = SequenceCache()

    stats = {
        'total_regions': len(regulatory_regions),
        'extracted': 0,
        'with_assoc_gene': 0,
        'with_nearest_gene': 0,
        'starr_seq': 0,
        'skipped_chrom': 0,
        'skipped_repressed': 0,
        'sources': defaultdict(int),
        'distance_distribution': {
            'overlapping': 0,
            'within_1kb': 0,
            'within_10kb': 0,
            'within_100kb': 0,
            'beyond_100kb': 0
        }
    }

    for reg in regulatory_regions:
        chrom = reg['chrom']

        # Filter repressed regions (silencers) if requested
        if exclude_repressed and 'Repressed' in reg['name']:
            stats['skipped_repressed'] += 1
            continue

        # Filter chromosomes
        if major_chroms_only and chrom not in major_chroms:
            stats['skipped_chrom'] += 1
            continue

        if chrom not in genome:
            stats['skipped_chrom'] += 1
            continue

        chrom_seq = genome[chrom]
        start, end = reg['start'], reg['end']

        # Check bounds
        if start < 1 or end > len(chrom_seq):
            continue

        # Extract sequence
        seq_start = start - 1
        enhancer_seq = chrom_seq[seq_start:end]

        stats['extracted'] += 1
        stats['sources'][reg['source']] += 1

        if reg['is_starr']:
            stats['starr_seq'] += 1

        # Determine gene association
        if reg['assoc_fbgn']:
            # Use existing association
            fbgn = reg['assoc_fbgn']
            symbol = reg['assoc_symbol']
            distance = 0  # Assumed associated
            gene_source = 'annotated'
            stats['with_assoc_gene'] += 1
        else:
            # Find nearest gene
            fbgn, symbol, distance = gene_index.find_nearest(chrom, start, end)
            gene_source = 'nearest'
            if fbgn:
                stats['with_nearest_gene'] += 1

        # Track distance distribution
        if fbgn:
            if distance == 0:
                stats['distance_distribution']['overlapping'] += 1
            elif distance <= 1000:
                stats['distance_distribution']['within_1kb'] += 1
            elif distance <= 10000:
                stats['distance_distribution']['within_10kb'] += 1
            elif distance <= 100000:
                stats['distance_distribution']['within_100kb'] += 1
            else:
                stats['distance_distribution']['beyond_100kb'] += 1

        # Build record
        gene_str = f"{symbol}:{fbgn}" if fbgn else "none"
        dist_str = str(distance) if distance >= 0 else "NA"

        description = (
            f"name={reg['name']} nearest_gene={gene_str} dist={dist_str} "
            f"gene_source={gene_source} loc={chrom}:{start}-{end} "
            f"is_starr={reg['is_starr']} length={len(enhancer_seq)}"
        )

        record = SeqRecord(
            Seq(enhancer_seq),
            id=reg['reg_id'],
            description=description
        )

        if deduplicate:
            if cache.add(record):
                all_records.append(record)
        else:
            all_records.append(record)

        # Store metadata
        metadata_rows.append({
            'enhancer_id': reg['reg_id'],
            'name': reg['name'],
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': len(enhancer_seq),
            'source': reg['source'],
            'is_starr_seq': reg['is_starr'],
            'assoc_fbgn': reg['assoc_fbgn'] or '',
            'assoc_symbol': reg['assoc_symbol'] or '',
            'nearest_fbgn': fbgn or '',
            'nearest_symbol': symbol or '',
            'distance_to_gene': distance if distance >= 0 else '',
            'gene_source': gene_source if fbgn else ''
        })

    if verbose:
        print(f"\nExtraction summary:")
        print(f"  Total regions: {stats['total_regions']:,}")
        print(f"  Extracted: {stats['extracted']:,}")
        print(f"  STARR-seq validated: {stats['starr_seq']:,}")
        print(f"  With annotated gene: {stats['with_assoc_gene']:,}")
        print(f"  With nearest gene assigned: {stats['with_nearest_gene']:,}")
        print(f"  Skipped (chromosome): {stats['skipped_chrom']:,}")
        print(f"  Skipped (repressed/silencer): {stats['skipped_repressed']:,}")
        print(f"\n  Distance to gene:")
        for dist_cat, count in stats['distance_distribution'].items():
            print(f"    {dist_cat}: {count:,}")
        print(f"\n  By source (top 10):")
        for src, count in sorted(stats['sources'].items(), key=lambda x: -x[1])[:10]:
            print(f"    {src}: {count:,}")

    # Get final records
    if deduplicate:
        dedup_stats = cache.get_stats()
        final_records = cache.get_unique_records()
        if verbose:
            print(f"\n  Deduplication:")
            print(f"    Total regions: {dedup_stats['total_regions']:,}")
            print(f"    Unique sequences: {dedup_stats['unique_sequences']:,}")
            print(f"    Duplicates removed: {dedup_stats['duplicates_skipped']:,}")
    else:
        final_records = all_records

    # Write FASTA
    fasta_path = output_dir / 'enhancers_sense.fasta'
    SeqIO.write(final_records, fasta_path, 'fasta')
    if verbose:
        print(f"\nWrote {len(final_records):,} sequences to: {fasta_path}")

    # Write metadata TSV
    metadata_path = output_dir / 'enhancer_metadata.tsv'
    with open(metadata_path, 'w') as f:
        header = ['enhancer_id', 'name', 'chrom', 'start', 'end', 'length',
                  'source', 'is_starr_seq', 'assoc_fbgn', 'assoc_symbol',
                  'nearest_fbgn', 'nearest_symbol', 'distance_to_gene', 'gene_source']
        f.write('\t'.join(header) + '\n')
        for row in metadata_rows:
            f.write('\t'.join(str(row[col]) for col in header) + '\n')
    if verbose:
        print(f"Wrote metadata to: {metadata_path}")

    # Write stats
    final_stats = {
        'total_regions': stats['total_regions'],
        'extracted': stats['extracted'],
        'unique_sequences': len(final_records),
        'starr_seq_validated': stats['starr_seq'],
        'with_annotated_gene': stats['with_assoc_gene'],
        'with_nearest_gene': stats['with_nearest_gene'],
        'skipped_chromosome': stats['skipped_chrom'],
        'skipped_repressed': stats['skipped_repressed'],
        'distance_distribution': stats['distance_distribution'],
        'sources': dict(stats['sources']),
        'deduplicated': deduplicate
    }

    if deduplicate:
        final_stats['duplicates_removed'] = dedup_stats['duplicates_skipped']

    stats_path = output_dir / 'extraction_stats.json'
    with open(stats_path, 'w') as f:
        json.dump(final_stats, f, indent=2)
    if verbose:
        print(f"Wrote stats to: {stats_path}")

    return final_stats


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--gff',
        type=Path,
        default=get_references_dir() / 'dmel-all-r6.66.gff',
        help='GFF annotation file'
    )
    parser.add_argument(
        '--genome',
        type=Path,
        default=get_references_dir() / 'dmel_genome.fasta',
        help='Genome FASTA file'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=get_queries_dir() / 'regulatory' / 'enhancers',
        help='Output directory'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--no-deduplicate',
        action='store_true',
        help='Disable deduplication'
    )
    parser.add_argument(
        '--all-chroms',
        action='store_true',
        help='Include all chromosomes/scaffolds'
    )
    parser.add_argument(
        '--include-repressed',
        action='store_true',
        help='Include repressed regions (normally excluded as silencers)'
    )

    args = parser.parse_args()

    if not args.gff.exists():
        print(f"Error: GFF file not found: {args.gff}", file=sys.stderr)
        return 1

    if not args.genome.exists():
        print(f"Error: Genome file not found: {args.genome}", file=sys.stderr)
        return 1

    print("Extracting Enhancer/CRM Regions")
    print("=" * 60)
    print()

    stats = extract_enhancers(
        args.gff,
        args.genome,
        args.output_dir,
        verbose=args.verbose,
        deduplicate=not args.no_deduplicate,
        major_chroms_only=not args.all_chroms,
        exclude_repressed=not args.include_repressed
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  Regions extracted: {stats['extracted']:,}")
    print(f"  Unique sequences: {stats['unique_sequences']:,}")
    print(f"  STARR-seq validated: {stats['starr_seq_validated']:,}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
