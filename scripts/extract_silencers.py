#!/usr/bin/env python3
"""
Extract silencer/repressive element sequences.

Uses regulatory_region features with "Repressed" in name, which come from
STARR-seq experiments identifying regions that decrease reporter activity.

IMPORTANT CAVEAT:
This dataset only includes 331 repressed regions from STARR-seq negative
controls in Kc167 cells. For comprehensive genome-wide silencer/repressive
chromatin annotation, H3K27me3 ChIP-seq data would be more appropriate.
The current dataset should be considered a limited but high-confidence set.

Sources:
- Kc167-ST (STARR-seq) Repressed regions

Outputs:
- silencers_sense.fasta: Silencer sequences
- silencer_metadata.tsv: Detailed metadata
- extraction_stats.json: Summary statistics with caveats noted

Header format:
>FBsf0000XXXXXX name=Unspecified_Kc167-ST_46_Repressed001 nearest_gene=XXX:FBgnYYY dist=500 loc=2L:1000-1500 length=500
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
    """Index for fast nearest-gene lookup."""

    def __init__(self):
        self.genes_by_chrom = defaultdict(list)
        self.built = False

    def add_gene(self, chrom: str, start: int, end: int, fbgn: str, symbol: str):
        midpoint = (start + end) // 2
        self.genes_by_chrom[chrom].append((midpoint, start, end, fbgn, symbol))

    def build(self):
        for chrom in self.genes_by_chrom:
            self.genes_by_chrom[chrom].sort(key=lambda x: x[0])
        self.built = True

    def find_nearest(self, chrom: str, region_start: int, region_end: int) -> Tuple[Optional[str], Optional[str], int]:
        if not self.built:
            self.build()

        genes = self.genes_by_chrom.get(chrom, [])
        if not genes:
            return None, None, -1

        region_mid = (region_start + region_end) // 2
        midpoints = [g[0] for g in genes]
        idx = bisect_left(midpoints, region_mid)

        candidates = []
        for i in range(max(0, idx - 2), min(len(genes), idx + 3)):
            gene_mid, gene_start, gene_end, fbgn, symbol = genes[i]

            if region_end < gene_start:
                dist = gene_start - region_end
            elif region_start > gene_end:
                dist = region_start - gene_end
            else:
                dist = 0

            candidates.append((dist, fbgn, symbol))

        if candidates:
            candidates.sort(key=lambda x: x[0])
            dist, fbgn, symbol = candidates[0]
            return fbgn, symbol, dist

        return None, None, -1


def load_genes_and_silencers(gff_path: Path, verbose: bool = False) -> Tuple[GeneIndex, List[dict]]:
    """
    Parse GFF to extract gene coordinates and silencer (Repressed) regions.
    """
    gene_index = GeneIndex()
    silencer_regions = []

    line_count = 0
    gene_count = 0
    silencer_count = 0

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
                name = attrs.get('Name', '')

                # Only include "Repressed" regions (silencers)
                if 'Repressed' not in name:
                    continue

                reg_id = attrs.get('ID', '')

                silencer_regions.append({
                    'reg_id': reg_id,
                    'name': name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'source': source
                })
                silencer_count += 1

    gene_index.build()

    if verbose:
        print(f"  Total lines: {line_count:,}")
        print(f"  Genes: {gene_count:,}")
        print(f"  Silencer (Repressed) regions: {silencer_count:,}")

    return gene_index, silencer_regions


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


def extract_silencers(
    gff_path: Path,
    genome_path: Path,
    output_dir: Path,
    verbose: bool = False,
    deduplicate: bool = True,
    major_chroms_only: bool = True
) -> dict:
    """
    Extract silencer sequences with nearest gene annotation.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    major_chroms = {'2L', '2R', '3L', '3R', 'X', '4'}

    # Load annotations
    gene_index, silencer_regions = load_genes_and_silencers(gff_path, verbose)

    # Load genome
    genome = load_genome(genome_path, verbose)

    if verbose:
        print(f"\nExtracting silencer sequences...")
        print(f"  WARNING: Only {len(silencer_regions)} repressed regions available.")
        print(f"  For comprehensive silencer annotation, consider H3K27me3 ChIP-seq data.")

    all_records = []
    metadata_rows = []
    cache = SequenceCache()

    stats = {
        'total_regions': len(silencer_regions),
        'extracted': 0,
        'with_nearest_gene': 0,
        'skipped_chrom': 0,
        'sources': defaultdict(int)
    }

    for reg in silencer_regions:
        chrom = reg['chrom']

        if major_chroms_only and chrom not in major_chroms:
            stats['skipped_chrom'] += 1
            continue

        if chrom not in genome:
            stats['skipped_chrom'] += 1
            continue

        chrom_seq = genome[chrom]
        start, end = reg['start'], reg['end']

        if start < 1 or end > len(chrom_seq):
            continue

        # Extract sequence
        seq_start = start - 1
        silencer_seq = chrom_seq[seq_start:end]

        stats['extracted'] += 1
        stats['sources'][reg['source']] += 1

        # Find nearest gene
        fbgn, symbol, distance = gene_index.find_nearest(chrom, start, end)
        if fbgn:
            stats['with_nearest_gene'] += 1

        # Build record
        gene_str = f"{symbol}:{fbgn}" if fbgn else "none"
        dist_str = str(distance) if distance >= 0 else "NA"

        description = (
            f"name={reg['name']} nearest_gene={gene_str} dist={dist_str} "
            f"loc={chrom}:{start}-{end} length={len(silencer_seq)}"
        )

        record = SeqRecord(
            Seq(silencer_seq),
            id=reg['reg_id'],
            description=description
        )

        if deduplicate:
            if cache.add(record):
                all_records.append(record)
        else:
            all_records.append(record)

        metadata_rows.append({
            'silencer_id': reg['reg_id'],
            'name': reg['name'],
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': len(silencer_seq),
            'source': reg['source'],
            'nearest_fbgn': fbgn or '',
            'nearest_symbol': symbol or '',
            'distance_to_gene': distance if distance >= 0 else ''
        })

    if verbose:
        print(f"\nExtraction summary:")
        print(f"  Total regions: {stats['total_regions']:,}")
        print(f"  Extracted: {stats['extracted']:,}")
        print(f"  With nearest gene: {stats['with_nearest_gene']:,}")
        print(f"  Skipped (chromosome): {stats['skipped_chrom']:,}")

    # Get final records
    if deduplicate:
        dedup_stats = cache.get_stats()
        final_records = cache.get_unique_records()
        if verbose:
            print(f"\n  Deduplication:")
            print(f"    Unique sequences: {dedup_stats['unique_sequences']:,}")
            print(f"    Duplicates removed: {dedup_stats['duplicates_skipped']:,}")
    else:
        final_records = all_records

    # Write FASTA
    fasta_path = output_dir / 'silencers_sense.fasta'
    SeqIO.write(final_records, fasta_path, 'fasta')
    if verbose:
        print(f"\nWrote {len(final_records):,} sequences to: {fasta_path}")

    # Write metadata TSV
    metadata_path = output_dir / 'silencer_metadata.tsv'
    with open(metadata_path, 'w') as f:
        header = ['silencer_id', 'name', 'chrom', 'start', 'end', 'length',
                  'source', 'nearest_fbgn', 'nearest_symbol', 'distance_to_gene']
        f.write('\t'.join(header) + '\n')
        for row in metadata_rows:
            f.write('\t'.join(str(row[col]) for col in header) + '\n')
    if verbose:
        print(f"Wrote metadata to: {metadata_path}")

    # Write stats with caveat
    final_stats = {
        'total_regions': stats['total_regions'],
        'extracted': stats['extracted'],
        'unique_sequences': len(final_records),
        'with_nearest_gene': stats['with_nearest_gene'],
        'skipped_chromosome': stats['skipped_chrom'],
        'sources': dict(stats['sources']),
        'deduplicated': deduplicate,
        'caveats': [
            'Limited to STARR-seq repressed regions only (n=331)',
            'For comprehensive silencer/repressive chromatin annotation, H3K27me3 ChIP-seq data is recommended',
            'These regions represent high-confidence silencers but likely miss many genome-wide'
        ]
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
        default=get_queries_dir() / 'regulatory' / 'silencers',
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

    args = parser.parse_args()

    if not args.gff.exists():
        print(f"Error: GFF file not found: {args.gff}", file=sys.stderr)
        return 1

    if not args.genome.exists():
        print(f"Error: Genome file not found: {args.genome}", file=sys.stderr)
        return 1

    print("Extracting Silencer (Repressed) Regions")
    print("=" * 60)
    print()
    print("NOTE: This dataset is limited to STARR-seq repressed regions.")
    print("For comprehensive silencer annotation, H3K27me3 ChIP-seq data")
    print("would provide better genome-wide coverage.")
    print()

    stats = extract_silencers(
        args.gff,
        args.genome,
        args.output_dir,
        verbose=args.verbose,
        deduplicate=not args.no_deduplicate,
        major_chroms_only=not args.all_chroms
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  Regions extracted: {stats['extracted']:,}")
    print(f"  Unique sequences: {stats['unique_sequences']:,}")
    print()
    print("Caveats:")
    for caveat in stats.get('caveats', []):
        print(f"  - {caveat}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
