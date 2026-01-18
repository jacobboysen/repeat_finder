#!/usr/bin/env python3
"""
Extract transcription factor binding site sequences.

Uses TF_binding_site features from FlyBase GFF, which includes ChIP-seq
and ChIP-chip data from multiple sources:

- BDTNP: Berkeley Drosophila Transcription Network Project (early embryo)
- modENCODE (mE1): Various TFs across developmental stages
- Timed embryo ChIP: E2-4h, E4-6h, E6-8h, E8-10h, E10-12h

Key embryo TFs:
- dl (dorsal), twi (twist), sna (snail), hb (hunchback) - early patterning
- Mef2 - muscle development
- pnr (pannier) - heart/dorsal development

Outputs:
- tfbs_sense.fasta: TFBS sequences (deduplicated by coordinates)
- tfbs_metadata.tsv: Detailed metadata per binding site
- extraction_stats.json: Summary statistics

Header format:
>FBsf0000295671 tf=dl tf_fbgn=FBgn0260632 library=BDTNP1_TFBS_dl nearest_gene=gene:FBgn loc=2L:5203-6186 length=984
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


class CoordinateDeduplicator:
    """
    Deduplicate by genomic coordinates.

    Multiple TFs can bind the same region; we deduplicate by exact coordinates
    but keep track of which TFs bind each region.
    """

    def __init__(self):
        self.seen_coords = {}  # (chrom, start, end) -> record
        self.coord_to_tfs = defaultdict(set)  # (chrom, start, end) -> set of TF names
        self.duplicates_skipped = 0

    def add(self, record: SeqRecord, chrom: str, start: int, end: int, tf: str) -> bool:
        """
        Add a region, returning True if new coordinates, False if duplicate.
        """
        key = (chrom, start, end)
        self.coord_to_tfs[key].add(tf)

        if key in self.seen_coords:
            self.duplicates_skipped += 1
            return False
        else:
            self.seen_coords[key] = record
            return True

    def get_unique_records(self) -> List[SeqRecord]:
        return list(self.seen_coords.values())

    def get_tf_count(self, chrom: str, start: int, end: int) -> int:
        """Get number of TFs binding this region."""
        return len(self.coord_to_tfs.get((chrom, start, end), set()))

    def get_stats(self) -> dict:
        return {
            'unique_regions': len(self.seen_coords),
            'duplicates_skipped': self.duplicates_skipped,
            'multi_tf_regions': sum(1 for tfs in self.coord_to_tfs.values() if len(tfs) > 1)
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


def load_genes_and_tfbs(gff_path: Path, verbose: bool = False) -> Tuple[GeneIndex, List[dict]]:
    """
    Parse GFF to extract gene coordinates and TF_binding_site features.
    """
    gene_index = GeneIndex()
    tfbs_regions = []

    line_count = 0
    gene_count = 0
    tfbs_count = 0

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

            elif feature_type == 'TF_binding_site':
                tfbs_id = attrs.get('ID', '')
                name = attrs.get('Name', '')

                # Parse bound_moiety (the TF)
                tf_symbol = None
                tf_fbgn = None
                if 'bound_moiety' in attrs:
                    # Format: "tf_symbol:FBgnXXXXXXX"
                    bound = attrs['bound_moiety'].split(',')[0]  # Take first if multiple
                    if ':' in bound:
                        tf_symbol, tf_fbgn = bound.split(':', 1)

                # Parse library to get source/timepoint
                library = attrs.get('library', '')

                # Extract embryo timepoint if present
                timepoint = None
                if '_E' in source and 'h_' in source:
                    # e.g., ChIP_chip_twi_E2_4h_organism -> E2-4h
                    parts_src = source.split('_')
                    for i, p in enumerate(parts_src):
                        if p.startswith('E') and i + 1 < len(parts_src) and 'h' in parts_src[i + 1]:
                            timepoint = f"{p}-{parts_src[i + 1].replace('_', '')}"
                            break

                tfbs_regions.append({
                    'tfbs_id': tfbs_id,
                    'name': name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'source': source,
                    'library': library,
                    'tf_symbol': tf_symbol,
                    'tf_fbgn': tf_fbgn,
                    'timepoint': timepoint
                })
                tfbs_count += 1

    gene_index.build()

    if verbose:
        print(f"  Total lines: {line_count:,}")
        print(f"  Genes: {gene_count:,}")
        print(f"  TFBS regions: {tfbs_count:,}")

    return gene_index, tfbs_regions


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


def extract_tfbs(
    gff_path: Path,
    genome_path: Path,
    output_dir: Path,
    verbose: bool = False,
    deduplicate: bool = True,
    major_chroms_only: bool = True
) -> dict:
    """
    Extract TFBS sequences with coordinate-based deduplication.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    major_chroms = {'2L', '2R', '3L', '3R', 'X', '4'}

    # Load annotations
    gene_index, tfbs_regions = load_genes_and_tfbs(gff_path, verbose)

    # Load genome
    genome = load_genome(genome_path, verbose)

    if verbose:
        print(f"\nExtracting TFBS sequences...")

    all_records = []
    metadata_rows = []
    deduplicator = CoordinateDeduplicator()

    stats = {
        'total_regions': len(tfbs_regions),
        'extracted': 0,
        'unique_coords': 0,
        'with_nearest_gene': 0,
        'skipped_chrom': 0,
        'tfs': defaultdict(int),
        'sources': defaultdict(int),
        'timepoints': defaultdict(int)
    }

    for tfbs in tfbs_regions:
        chrom = tfbs['chrom']

        if major_chroms_only and chrom not in major_chroms:
            stats['skipped_chrom'] += 1
            continue

        if chrom not in genome:
            stats['skipped_chrom'] += 1
            continue

        chrom_seq = genome[chrom]
        start, end = tfbs['start'], tfbs['end']

        if start < 1 or end > len(chrom_seq):
            continue

        # Extract sequence
        seq_start = start - 1
        tfbs_seq = chrom_seq[seq_start:end]

        stats['extracted'] += 1
        stats['tfs'][tfbs['tf_symbol'] or 'unknown'] += 1
        stats['sources'][tfbs['source']] += 1
        if tfbs['timepoint']:
            stats['timepoints'][tfbs['timepoint']] += 1

        # Find nearest gene
        fbgn, symbol, distance = gene_index.find_nearest(chrom, start, end)
        if fbgn:
            stats['with_nearest_gene'] += 1

        # Build record
        tf_str = tfbs['tf_symbol'] or 'unknown'
        tf_fbgn_str = tfbs['tf_fbgn'] or 'none'
        gene_str = f"{symbol}:{fbgn}" if fbgn else "none"
        dist_str = str(distance) if distance >= 0 else "NA"

        description = (
            f"tf={tf_str} tf_fbgn={tf_fbgn_str} library={tfbs['source']} "
            f"nearest_gene={gene_str} dist={dist_str} "
            f"loc={chrom}:{start}-{end} length={len(tfbs_seq)}"
        )

        record = SeqRecord(
            Seq(tfbs_seq),
            id=tfbs['tfbs_id'],
            description=description
        )

        if deduplicate:
            if deduplicator.add(record, chrom, start, end, tf_str):
                all_records.append(record)
                stats['unique_coords'] += 1
        else:
            all_records.append(record)

        # Store metadata
        metadata_rows.append({
            'tfbs_id': tfbs['tfbs_id'],
            'name': tfbs['name'],
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': len(tfbs_seq),
            'tf_symbol': tfbs['tf_symbol'] or '',
            'tf_fbgn': tfbs['tf_fbgn'] or '',
            'source': tfbs['source'],
            'library': tfbs['library'],
            'timepoint': tfbs['timepoint'] or '',
            'nearest_fbgn': fbgn or '',
            'nearest_symbol': symbol or '',
            'distance_to_gene': distance if distance >= 0 else ''
        })

    if verbose:
        print(f"\nExtraction summary:")
        print(f"  Total TFBS: {stats['total_regions']:,}")
        print(f"  Extracted: {stats['extracted']:,}")
        print(f"  With nearest gene: {stats['with_nearest_gene']:,}")
        print(f"  Skipped (chromosome): {stats['skipped_chrom']:,}")
        print(f"\n  By TF (top 15):")
        for tf, count in sorted(stats['tfs'].items(), key=lambda x: -x[1])[:15]:
            print(f"    {tf}: {count:,}")
        if stats['timepoints']:
            print(f"\n  By timepoint:")
            for tp, count in sorted(stats['timepoints'].items()):
                print(f"    {tp}: {count:,}")

    # Get final records
    if deduplicate:
        dedup_stats = deduplicator.get_stats()
        final_records = deduplicator.get_unique_records()
        if verbose:
            print(f"\n  Coordinate deduplication:")
            print(f"    Unique regions: {dedup_stats['unique_regions']:,}")
            print(f"    Duplicate coordinates: {dedup_stats['duplicates_skipped']:,}")
            print(f"    Multi-TF regions: {dedup_stats['multi_tf_regions']:,}")
    else:
        final_records = all_records

    # Write FASTA
    fasta_path = output_dir / 'tfbs_sense.fasta'
    SeqIO.write(final_records, fasta_path, 'fasta')
    if verbose:
        print(f"\nWrote {len(final_records):,} sequences to: {fasta_path}")

    # Write metadata TSV
    metadata_path = output_dir / 'tfbs_metadata.tsv'
    with open(metadata_path, 'w') as f:
        header = ['tfbs_id', 'name', 'chrom', 'start', 'end', 'length',
                  'tf_symbol', 'tf_fbgn', 'source', 'library', 'timepoint',
                  'nearest_fbgn', 'nearest_symbol', 'distance_to_gene']
        f.write('\t'.join(header) + '\n')
        for row in metadata_rows:
            f.write('\t'.join(str(row[col]) for col in header) + '\n')
    if verbose:
        print(f"Wrote metadata to: {metadata_path}")

    # Write stats
    final_stats = {
        'total_regions': stats['total_regions'],
        'extracted': stats['extracted'],
        'unique_coordinate_regions': len(final_records) if deduplicate else stats['extracted'],
        'with_nearest_gene': stats['with_nearest_gene'],
        'skipped_chromosome': stats['skipped_chrom'],
        'tfs': dict(stats['tfs']),
        'sources': dict(stats['sources']),
        'timepoints': dict(stats['timepoints']),
        'deduplicated': deduplicate
    }

    if deduplicate:
        final_stats['duplicates_removed'] = dedup_stats['duplicates_skipped']
        final_stats['multi_tf_regions'] = dedup_stats['multi_tf_regions']

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
        default=get_queries_dir() / 'regulatory' / 'tfbs',
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
        help='Disable coordinate deduplication'
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

    print("Extracting TF Binding Site Regions")
    print("=" * 60)
    print()

    stats = extract_tfbs(
        args.gff,
        args.genome,
        args.output_dir,
        verbose=args.verbose,
        deduplicate=not args.no_deduplicate,
        major_chroms_only=not args.all_chroms
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  TFBS extracted: {stats['extracted']:,}")
    print(f"  Unique coordinate regions: {stats['unique_coordinate_regions']:,}")
    print(f"  TFs represented: {len(stats['tfs'])}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
