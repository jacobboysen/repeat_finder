#!/usr/bin/env python3
"""
Extract individual exons from GFF annotations and genome FASTA.

Unlike the UTR extraction scripts which use pre-extracted FlyBase FASTAs,
this script parses the GFF directly to get individual exon coordinates,
enabling position tracking (first/last/internal) and UTR overlap annotation.

Outputs:
- exons_sense.fasta: Individual exon sequences with position annotations
- exon_metadata.tsv: Detailed metadata for each exon
- extraction_stats.json: Summary statistics

Header format:
>FBtr0346766_exon1 gene=DIP-lambda FBgn=FBgn0267428 position=first_exon utr_overlap=5utr loc=2R:326768..326939:+ length=172
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


class SequenceCache:
    """Cache for deduplicating identical sequences from different exons."""

    def __init__(self):
        self.seen_hashes = {}  # hash -> first SeqRecord
        self.duplicates = defaultdict(list)  # hash -> list of duplicate IDs
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
            'total_exons': sum(len(v) for v in self.duplicates.values())
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


def load_gff_annotations(gff_path: Path, verbose: bool = False) -> Tuple[dict, dict, dict, dict, dict]:
    """
    Parse GFF file to extract exon, mRNA, and UTR information.

    Returns:
        exons_by_transcript: {FBtr: [(chrom, start, end, strand), ...]}
        transcript_to_gene: {FBtr: FBgn}
        gene_symbols: {FBgn: symbol}
        utr5_by_transcript: {FBtr: [(start, end), ...]}
        utr3_by_transcript: {FBtr: [(start, end), ...]}
    """
    exons_by_transcript = defaultdict(list)
    transcript_to_gene = {}
    gene_symbols = {}
    utr5_by_transcript = defaultdict(list)
    utr3_by_transcript = defaultdict(list)

    line_count = 0
    exon_count = 0
    mrna_count = 0
    gene_count = 0

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

            if feature_type == 'exon':
                parent = attrs.get('Parent', '')
                if parent.startswith('FBtr'):
                    exons_by_transcript[parent].append((chrom, start, end, strand))
                    exon_count += 1

            elif feature_type == 'mRNA':
                fbtr = attrs.get('ID', '')
                parent = attrs.get('Parent', '')
                name = attrs.get('Name', '')
                if fbtr.startswith('FBtr') and parent.startswith('FBgn'):
                    transcript_to_gene[fbtr] = parent
                    mrna_count += 1

            elif feature_type == 'gene':
                fbgn = attrs.get('ID', '')
                name = attrs.get('Name', '')
                if fbgn.startswith('FBgn'):
                    gene_symbols[fbgn] = name
                    gene_count += 1

            elif feature_type == 'five_prime_UTR':
                parent = attrs.get('Parent', '')
                if parent.startswith('FBtr'):
                    utr5_by_transcript[parent].append((start, end))

            elif feature_type == 'three_prime_UTR':
                parent = attrs.get('Parent', '')
                if parent.startswith('FBtr'):
                    utr3_by_transcript[parent].append((start, end))

    if verbose:
        print(f"  Total lines: {line_count:,}")
        print(f"  Genes: {gene_count:,}")
        print(f"  Transcripts (mRNA): {mrna_count:,}")
        print(f"  Exons: {exon_count:,}")
        print(f"  5'UTRs: {sum(len(v) for v in utr5_by_transcript.values()):,}")
        print(f"  3'UTRs: {sum(len(v) for v in utr3_by_transcript.values()):,}")

    return exons_by_transcript, transcript_to_gene, gene_symbols, utr5_by_transcript, utr3_by_transcript


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


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Check if two intervals overlap."""
    return a_start <= b_end and b_start <= a_end


def classify_exon_position(exon_idx: int, total_exons: int, strand: str) -> str:
    """
    Classify exon as first, last, internal, or single.

    For + strand: first exon is lowest coordinate
    For - strand: first exon is highest coordinate (5' end of gene)
    """
    if total_exons == 1:
        return 'single_exon'

    if strand == '+':
        if exon_idx == 0:
            return 'first_exon'
        elif exon_idx == total_exons - 1:
            return 'last_exon'
        else:
            return 'internal_exon'
    else:  # strand == '-'
        if exon_idx == total_exons - 1:
            return 'first_exon'
        elif exon_idx == 0:
            return 'last_exon'
        else:
            return 'internal_exon'


def get_utr_overlap(exon_start: int, exon_end: int,
                    utr5_regions: List[Tuple[int, int]],
                    utr3_regions: List[Tuple[int, int]]) -> str:
    """Determine if exon overlaps with UTR regions."""
    overlaps_5utr = any(intervals_overlap(exon_start, exon_end, u[0], u[1]) for u in utr5_regions)
    overlaps_3utr = any(intervals_overlap(exon_start, exon_end, u[0], u[1]) for u in utr3_regions)

    if overlaps_5utr and overlaps_3utr:
        return 'both'
    elif overlaps_5utr:
        return '5utr'
    elif overlaps_3utr:
        return '3utr'
    else:
        return 'none'


def extract_exons(
    gff_path: Path,
    genome_path: Path,
    output_dir: Path,
    verbose: bool = False,
    deduplicate: bool = True
) -> dict:
    """
    Extract individual exon sequences with annotations.

    Returns summary statistics.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GFF annotations
    exons_by_transcript, transcript_to_gene, gene_symbols, utr5_by_transcript, utr3_by_transcript = \
        load_gff_annotations(gff_path, verbose)

    # Load genome
    genome = load_genome(genome_path, verbose)

    if verbose:
        print("\nExtracting exon sequences...")

    # Collect all exon records
    all_records = []
    metadata_rows = []
    cache = SequenceCache()

    stats = {
        'total_transcripts': 0,
        'total_exons': 0,
        'position_counts': defaultdict(int),
        'utr_overlap_counts': defaultdict(int),
        'extraction_errors': 0
    }

    for fbtr, exons in exons_by_transcript.items():
        if fbtr not in transcript_to_gene:
            continue

        fbgn = transcript_to_gene[fbtr]
        gene_symbol = gene_symbols.get(fbgn, fbgn)

        # Sort exons by start position
        sorted_exons = sorted(exons, key=lambda x: x[1])
        total_exons = len(sorted_exons)

        # Get strand from first exon
        strand = sorted_exons[0][3] if sorted_exons else '+'

        # Get UTR regions for this transcript
        utr5_regions = utr5_by_transcript.get(fbtr, [])
        utr3_regions = utr3_by_transcript.get(fbtr, [])

        stats['total_transcripts'] += 1

        for exon_idx, (chrom, start, end, exon_strand) in enumerate(sorted_exons):
            stats['total_exons'] += 1

            # Get chromosome sequence
            if chrom not in genome:
                stats['extraction_errors'] += 1
                continue

            chrom_seq = genome[chrom]

            # Extract sequence (GFF is 1-based, inclusive)
            # Convert to 0-based for Python slicing
            seq_start = start - 1
            seq_end = end

            if seq_start < 0 or seq_end > len(chrom_seq):
                stats['extraction_errors'] += 1
                continue

            exon_seq = chrom_seq[seq_start:seq_end]

            # Reverse complement if on minus strand
            if exon_strand == '-':
                exon_seq = str(Seq(exon_seq).reverse_complement())

            # Classify exon position
            position = classify_exon_position(exon_idx, total_exons, exon_strand)
            stats['position_counts'][position] += 1

            # Check UTR overlap
            utr_overlap = get_utr_overlap(start, end, utr5_regions, utr3_regions)
            stats['utr_overlap_counts'][utr_overlap] += 1

            # Create exon ID (1-indexed)
            exon_num = exon_idx + 1
            exon_id = f"{fbtr}_exon{exon_num}"

            # Build header
            description = (
                f"gene={gene_symbol} FBgn={fbgn} position={position} "
                f"utr_overlap={utr_overlap} loc={chrom}:{start}..{end}:{exon_strand} "
                f"length={len(exon_seq)}"
            )

            record = SeqRecord(
                Seq(exon_seq),
                id=exon_id,
                description=description
            )

            # Add to cache (handles deduplication)
            if deduplicate:
                if cache.add(record):
                    all_records.append(record)
            else:
                all_records.append(record)

            # Store metadata
            metadata_rows.append({
                'exon_id': exon_id,
                'fbtr': fbtr,
                'fbgn': fbgn,
                'gene_symbol': gene_symbol,
                'exon_num': exon_num,
                'total_exons': total_exons,
                'position': position,
                'utr_overlap': utr_overlap,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': exon_strand,
                'length': len(exon_seq)
            })

    if verbose:
        print(f"\nExtraction summary:")
        print(f"  Transcripts processed: {stats['total_transcripts']:,}")
        print(f"  Total exons: {stats['total_exons']:,}")
        print(f"  Extraction errors: {stats['extraction_errors']}")
        print(f"\n  Position distribution:")
        for pos, count in sorted(stats['position_counts'].items()):
            print(f"    {pos}: {count:,}")
        print(f"\n  UTR overlap distribution:")
        for overlap, count in sorted(stats['utr_overlap_counts'].items()):
            print(f"    {overlap}: {count:,}")

    # Get records (deduplicated or all)
    if deduplicate:
        dedup_stats = cache.get_stats()
        final_records = cache.get_unique_records()
        if verbose:
            print(f"\n  Deduplication:")
            print(f"    Total exons: {dedup_stats['total_exons']:,}")
            print(f"    Unique sequences: {dedup_stats['unique_sequences']:,}")
            print(f"    Duplicates removed: {dedup_stats['duplicates_skipped']:,}")
    else:
        final_records = all_records

    # Write FASTA
    fasta_path = output_dir / 'exons_sense.fasta'
    SeqIO.write(final_records, fasta_path, 'fasta')
    if verbose:
        print(f"\nWrote {len(final_records):,} sequences to: {fasta_path}")

    # Write metadata TSV
    metadata_path = output_dir / 'exon_metadata.tsv'
    with open(metadata_path, 'w') as f:
        header = ['exon_id', 'fbtr', 'fbgn', 'gene_symbol', 'exon_num', 'total_exons',
                  'position', 'utr_overlap', 'chrom', 'start', 'end', 'strand', 'length']
        f.write('\t'.join(header) + '\n')
        for row in metadata_rows:
            f.write('\t'.join(str(row[col]) for col in header) + '\n')
    if verbose:
        print(f"Wrote metadata to: {metadata_path}")

    # Write stats JSON
    final_stats = {
        'total_transcripts': stats['total_transcripts'],
        'total_exons': stats['total_exons'],
        'unique_sequences': len(final_records),
        'extraction_errors': stats['extraction_errors'],
        'position_counts': dict(stats['position_counts']),
        'utr_overlap_counts': dict(stats['utr_overlap_counts']),
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
        default=Path('data/references/dmel-all-r6.66.gff'),
        help='GFF annotation file (default: data/references/dmel-all-r6.66.gff)'
    )
    parser.add_argument(
        '--genome',
        type=Path,
        default=Path('data/references/dmel_genome.fasta'),
        help='Genome FASTA file (default: data/references/dmel_genome.fasta)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/queries/genome_wide'),
        help='Output directory (default: data/queries/genome_wide)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--no-deduplicate',
        action='store_true',
        help='Disable deduplication (keep all exons even if identical)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.gff.exists():
        print(f"Error: GFF file not found: {args.gff}", file=sys.stderr)
        return 1

    if not args.genome.exists():
        print(f"Error: Genome file not found: {args.genome}", file=sys.stderr)
        return 1

    print("Extracting Individual Exons")
    print("=" * 60)

    stats = extract_exons(
        args.gff,
        args.genome,
        args.output_dir,
        args.verbose,
        deduplicate=not args.no_deduplicate
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  Transcripts: {stats['total_transcripts']:,}")
    print(f"  Total exons: {stats['total_exons']:,}")
    print(f"  Unique sequences: {stats['unique_sequences']:,}")
    if stats.get('duplicates_removed'):
        print(f"  Duplicates removed: {stats['duplicates_removed']:,}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
