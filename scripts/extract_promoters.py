#!/usr/bin/env python3
"""
Extract extended promoter regions from experimentally determined TSS data.

Uses RAMPAGE and modENCODE TSS annotations from FlyBase GFF to define
promoter regions. Extended promoters are defined as TSS-2000 to TSS+500
(2.5kb total), with proper strand orientation.

Outputs:
- promoters_sense.fasta: Promoter sequences (sense strand)
- promoter_metadata.tsv: Detailed metadata for each promoter
- extraction_stats.json: Summary statistics

Header format:
>TSS_RAMPAGE_004769 gene=CR11023 FBgn=FBgn0031208 loc=2L:7435-9934:+ tss=7435-7519 length=2500
"""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add parent directory to path for utils
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_references_dir, get_queries_dir


class SequenceCache:
    """Cache for deduplicating identical sequences."""

    def __init__(self):
        self.seen_hashes = {}
        self.duplicates = defaultdict(list)
        self.duplicates_skipped = 0

    def add(self, record: SeqRecord) -> bool:
        """Add a sequence, returning True if unique, False if duplicate."""
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
        """Return list of unique sequences."""
        return list(self.seen_hashes.values())

    def get_stats(self) -> dict:
        """Return deduplication statistics."""
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


def load_tss_annotations(gff_path: Path, verbose: bool = False) -> List[dict]:
    """
    Parse GFF file to extract TSS features.

    Returns list of TSS records with:
        - tss_id: Feature ID (e.g., FBsf0000549311)
        - name: TSS name (e.g., TSS_RAMPAGE_004769)
        - chrom: Chromosome
        - tss_start, tss_end: TSS region coordinates
        - strand: + or -
        - source: Library source (TSS_RAMPAGE or mE_Transcription_Start_Sites)
        - fbgn: Associated gene FBgn (if available)
        - gene_symbol: Associated gene symbol (if available)
    """
    tss_records = []

    line_count = 0
    tss_count = 0

    if verbose:
        print(f"Parsing GFF for TSS features: {gff_path}")

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

            if feature_type != 'TSS':
                continue

            start, end = int(start), int(end)
            attrs = parse_gff_attributes(attributes)

            tss_id = attrs.get('ID', '')
            name = attrs.get('Name', '')

            # Parse associated genes if present
            fbgn = None
            gene_symbol = None
            if 'associated_genes' in attrs:
                # Format: "gene_symbol:FBgnXXXXXXX" or "gene_symbol:FBgnXXXXXXX,gene2:FBgnYYYYYYY"
                assoc = attrs['associated_genes'].split(',')[0]  # Take first gene
                if ':' in assoc:
                    gene_symbol, fbgn = assoc.split(':', 1)

            tss_records.append({
                'tss_id': tss_id,
                'name': name,
                'chrom': chrom,
                'tss_start': start,
                'tss_end': end,
                'strand': strand,
                'source': source,
                'fbgn': fbgn,
                'gene_symbol': gene_symbol
            })
            tss_count += 1

    if verbose:
        print(f"  Total lines: {line_count:,}")
        print(f"  TSS features found: {tss_count:,}")

        # Count by source
        sources = defaultdict(int)
        for rec in tss_records:
            sources[rec['source']] += 1
        print(f"  By source:")
        for src, count in sorted(sources.items(), key=lambda x: -x[1]):
            print(f"    {src}: {count:,}")

    return tss_records


def load_genome(genome_path: Path, verbose: bool = False) -> Dict[str, str]:
    """Load genome sequences into dictionary."""
    if verbose:
        print(f"Loading genome: {genome_path}")

    genome = {}
    for record in SeqIO.parse(genome_path, 'fasta'):
        genome[record.id] = str(record.seq)

    if verbose:
        print(f"  Loaded {len(genome)} chromosomes/scaffolds")
        total_bp = sum(len(s) for s in genome.values())
        print(f"  Total: {total_bp:,} bp")

    return genome


def extract_promoters(
    gff_path: Path,
    genome_path: Path,
    output_dir: Path,
    upstream: int = 2000,
    downstream: int = 500,
    verbose: bool = False,
    deduplicate: bool = True,
    major_chroms_only: bool = True
) -> dict:
    """
    Extract extended promoter sequences around TSS.

    Args:
        gff_path: Path to GFF annotation file
        genome_path: Path to genome FASTA
        output_dir: Output directory
        upstream: bp upstream of TSS (default: 2000)
        downstream: bp downstream of TSS (default: 500)
        verbose: Print progress
        deduplicate: Remove duplicate sequences
        major_chroms_only: Only use major chromosomes (2L, 2R, 3L, 3R, X, 4)

    Returns:
        Statistics dictionary
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    major_chroms = {'2L', '2R', '3L', '3R', 'X', '4'}

    # Load TSS annotations
    tss_records = load_tss_annotations(gff_path, verbose)

    # Load genome
    genome = load_genome(genome_path, verbose)

    if verbose:
        print(f"\nExtracting promoter sequences...")
        print(f"  Region: TSS-{upstream} to TSS+{downstream}")

    # Collect records
    all_records = []
    metadata_rows = []
    cache = SequenceCache()

    stats = {
        'total_tss': len(tss_records),
        'extracted': 0,
        'with_gene': 0,
        'skipped_chrom': 0,
        'skipped_bounds': 0,
        'sources': defaultdict(int)
    }

    for tss in tss_records:
        chrom = tss['chrom']
        strand = tss['strand']
        tss_start = tss['tss_start']
        tss_end = tss['tss_end']

        # Filter to major chromosomes if requested
        if major_chroms_only and chrom not in major_chroms:
            stats['skipped_chrom'] += 1
            continue

        if chrom not in genome:
            stats['skipped_chrom'] += 1
            continue

        chrom_seq = genome[chrom]
        chrom_len = len(chrom_seq)

        # Use TSS start as anchor point for promoter region
        # For + strand: promoter is upstream (left) of TSS
        # For - strand: promoter is upstream (right) of TSS
        tss_pos = tss_start  # Use start of TSS region as reference

        if strand == '+':
            # Upstream is to the left, downstream is to the right
            region_start = tss_pos - upstream
            region_end = tss_pos + downstream
        else:
            # For minus strand, upstream is to the right
            region_start = tss_pos - downstream
            region_end = tss_pos + upstream

        # Check bounds
        if region_start < 1 or region_end > chrom_len:
            stats['skipped_bounds'] += 1
            continue

        # Extract sequence (GFF is 1-based)
        seq_start = region_start - 1  # Convert to 0-based
        seq_end = region_end
        promoter_seq = chrom_seq[seq_start:seq_end]

        # Reverse complement for minus strand
        if strand == '-':
            promoter_seq = str(Seq(promoter_seq).reverse_complement())

        stats['extracted'] += 1
        stats['sources'][tss['source']] += 1

        if tss['fbgn']:
            stats['with_gene'] += 1

        # Build record
        gene_str = tss['gene_symbol'] if tss['gene_symbol'] else 'unknown'
        fbgn_str = tss['fbgn'] if tss['fbgn'] else 'none'

        description = (
            f"gene={gene_str} FBgn={fbgn_str} "
            f"loc={chrom}:{region_start}-{region_end}:{strand} "
            f"tss={tss_start}-{tss_end} length={len(promoter_seq)}"
        )

        record = SeqRecord(
            Seq(promoter_seq),
            id=tss['name'],
            description=description
        )

        if deduplicate:
            if cache.add(record):
                all_records.append(record)
        else:
            all_records.append(record)

        # Store metadata
        metadata_rows.append({
            'promoter_id': tss['name'],
            'tss_id': tss['tss_id'],
            'fbgn': tss['fbgn'] or '',
            'gene_symbol': tss['gene_symbol'] or '',
            'chrom': chrom,
            'region_start': region_start,
            'region_end': region_end,
            'tss_start': tss_start,
            'tss_end': tss_end,
            'strand': strand,
            'source': tss['source'],
            'length': len(promoter_seq),
            'upstream_bp': upstream,
            'downstream_bp': downstream
        })

    if verbose:
        print(f"\nExtraction summary:")
        print(f"  Total TSS: {stats['total_tss']:,}")
        print(f"  Extracted: {stats['extracted']:,}")
        print(f"  With gene association: {stats['with_gene']:,}")
        print(f"  Skipped (chromosome): {stats['skipped_chrom']:,}")
        print(f"  Skipped (out of bounds): {stats['skipped_bounds']:,}")
        print(f"\n  By source:")
        for src, count in sorted(stats['sources'].items(), key=lambda x: -x[1]):
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
    fasta_path = output_dir / 'promoters_sense.fasta'
    SeqIO.write(final_records, fasta_path, 'fasta')
    if verbose:
        print(f"\nWrote {len(final_records):,} sequences to: {fasta_path}")

    # Write metadata TSV
    metadata_path = output_dir / 'promoter_metadata.tsv'
    with open(metadata_path, 'w') as f:
        header = ['promoter_id', 'tss_id', 'fbgn', 'gene_symbol', 'chrom',
                  'region_start', 'region_end', 'tss_start', 'tss_end',
                  'strand', 'source', 'length', 'upstream_bp', 'downstream_bp']
        f.write('\t'.join(header) + '\n')
        for row in metadata_rows:
            f.write('\t'.join(str(row[col]) for col in header) + '\n')
    if verbose:
        print(f"Wrote metadata to: {metadata_path}")

    # Compile final stats
    final_stats = {
        'total_tss': stats['total_tss'],
        'extracted': stats['extracted'],
        'with_gene_association': stats['with_gene'],
        'unique_sequences': len(final_records),
        'skipped_chromosome': stats['skipped_chrom'],
        'skipped_bounds': stats['skipped_bounds'],
        'sources': dict(stats['sources']),
        'region_definition': f'TSS-{upstream} to TSS+{downstream}',
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
        default=get_queries_dir() / 'regulatory' / 'promoters',
        help='Output directory'
    )
    parser.add_argument(
        '--upstream',
        type=int,
        default=2000,
        help='bp upstream of TSS (default: 2000)'
    )
    parser.add_argument(
        '--downstream',
        type=int,
        default=500,
        help='bp downstream of TSS (default: 500)'
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
        help='Include all chromosomes/scaffolds (default: major chroms only)'
    )

    args = parser.parse_args()

    if not args.gff.exists():
        print(f"Error: GFF file not found: {args.gff}", file=sys.stderr)
        return 1

    if not args.genome.exists():
        print(f"Error: Genome file not found: {args.genome}", file=sys.stderr)
        return 1

    print("Extracting Extended Promoter Regions")
    print("=" * 60)
    print(f"Region: TSS-{args.upstream} to TSS+{args.downstream}")
    print()

    stats = extract_promoters(
        args.gff,
        args.genome,
        args.output_dir,
        upstream=args.upstream,
        downstream=args.downstream,
        verbose=args.verbose,
        deduplicate=not args.no_deduplicate,
        major_chroms_only=not args.all_chroms
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  TSS processed: {stats['extracted']:,}")
    print(f"  Unique sequences: {stats['unique_sequences']:,}")
    print(f"  With gene associations: {stats['with_gene_association']:,}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
