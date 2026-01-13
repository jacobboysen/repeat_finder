#!/usr/bin/env python3
"""
Extract 3'UTRs for germ plasm genes from FlyBase reference.

Reads the consolidated gene list and extracts corresponding 3'UTR sequences
from the FlyBase reference FASTA. Handles multiple transcript isoforms per gene.

Features:
- Deduplication: Identical sequences from different isoforms are cached
  and only included once (with all isoform IDs recorded)
- Supports both germ plasm and housekeeping gene extraction

Outputs:
- 3UTR_sense.fasta: Original strand sequences (deduplicated)
- 3UTR_antisense.fasta: Reverse complement (negative control)
- Tier1 subsets of both
"""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def hash_sequence(seq):
    """Generate hash for sequence deduplication."""
    return hashlib.md5(str(seq).upper().encode()).hexdigest()


class SequenceCache:
    """Cache for deduplicating identical sequences from different isoforms."""

    def __init__(self):
        self.seen_hashes = {}  # hash -> first SeqRecord
        self.isoform_map = defaultdict(list)  # hash -> list of isoform IDs
        self.duplicates_skipped = 0

    def add(self, record, gene_symbol):
        """Add a sequence, returning True if unique, False if duplicate."""
        seq_hash = hash_sequence(record.seq)

        if seq_hash in self.seen_hashes:
            # Duplicate sequence - record the isoform ID but don't add
            self.isoform_map[seq_hash].append(record.id)
            self.duplicates_skipped += 1
            return False
        else:
            # New unique sequence
            self.seen_hashes[seq_hash] = record
            self.isoform_map[seq_hash].append(record.id)
            return True

    def get_unique_records(self):
        """Return list of unique sequences."""
        return list(self.seen_hashes.values())

    def get_stats(self):
        """Return deduplication statistics."""
        return {
            'unique_sequences': len(self.seen_hashes),
            'duplicates_skipped': self.duplicates_skipped,
            'total_isoforms': sum(len(v) for v in self.isoform_map.values())
        }

    def get_isoform_map(self):
        """Return map of sequence hash to isoform IDs."""
        # Convert to use record ID as key for JSON output
        result = {}
        for seq_hash, record in self.seen_hashes.items():
            result[record.id] = {
                'hash': seq_hash,
                'all_isoforms': self.isoform_map[seq_hash],
                'num_isoforms': len(self.isoform_map[seq_hash])
            }
        return result


def load_gene_list(gene_list_file):
    """Load gene list from TSV file."""
    genes = {}
    tier1_genes = set()

    with open(gene_list_file, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                symbol = fields[0]
                fbgn_id = fields[1]
                tier = int(fields[3])

                genes[fbgn_id] = {
                    'symbol': symbol,
                    'tier': tier
                }

                if tier == 1:
                    tier1_genes.add(fbgn_id)

    return genes, tier1_genes


def extract_gene_id(description):
    """Extract FlyBase gene ID (FBgn) from FASTA description."""
    # Look for parent=FBgn pattern in description
    if 'parent=' in description:
        parts = description.split(';')
        for part in parts:
            if 'parent=' in part:
                fbgn = part.split('=')[1].strip()
                if fbgn.startswith('FBgn'):
                    return fbgn

    # Direct FBgn in ID
    for word in description.split():
        if word.startswith('FBgn'):
            return word.rstrip(';,')

    return None


def extract_transcript_id(record_id):
    """Extract transcript ID from record ID."""
    if record_id.startswith('FBtr'):
        return record_id
    return None


def reverse_complement(seq_record, suffix='_antisense'):
    """Create reverse complement of a SeqRecord."""
    rc_seq = seq_record.seq.reverse_complement()
    new_id = seq_record.id + suffix
    new_desc = seq_record.description.replace(seq_record.id, new_id)

    return SeqRecord(
        rc_seq,
        id=new_id,
        description=new_desc
    )


def extract_3utrs(reference_fasta, gene_list, output_dir, verbose=False, deduplicate=True):
    """Extract 3'UTRs for genes in the list with optional deduplication."""
    output_dir.mkdir(parents=True, exist_ok=True)

    genes, tier1_genes = load_gene_list(gene_list)

    if verbose:
        print(f"Loaded {len(genes)} genes from gene list")
        print(f"  Tier 1 genes: {len(tier1_genes)}")

    # Track results by gene
    extracted = defaultdict(list)
    total_seqs = 0

    # Deduplication caches
    all_cache = SequenceCache()
    tier1_cache = SequenceCache()

    # Parse 3'UTR reference
    print(f"\nParsing 3'UTRs from: {reference_fasta}")

    for record in SeqIO.parse(reference_fasta, 'fasta'):
        gene_id = extract_gene_id(record.description)

        if gene_id and gene_id in genes:
            gene_info = genes[gene_id]
            symbol = gene_info['symbol']

            # Rename record to include gene symbol
            transcript_id = extract_transcript_id(record.id)
            if transcript_id:
                new_id = f"{symbol}_{transcript_id}"
            else:
                new_id = f"{symbol}_{record.id}"

            new_record = SeqRecord(
                record.seq,
                id=new_id,
                description=f"gene={symbol} {record.description}"
            )

            extracted[gene_id].append(new_record)
            total_seqs += 1

            # Add to caches (handles deduplication)
            all_cache.add(new_record, symbol)
            if gene_id in tier1_genes:
                tier1_cache.add(new_record, symbol)

    print(f"  Found {total_seqs} 3'UTR sequences for {len(extracted)} genes")

    # Check for missing genes
    missing = set(genes.keys()) - set(extracted.keys())
    if missing:
        print(f"\n  Warning: {len(missing)} genes not found in reference:")
        for fbgn in missing:
            print(f"    - {genes[fbgn]['symbol']} ({fbgn})")

    # Get unique records (deduplicated) or all records
    if deduplicate:
        all_records = all_cache.get_unique_records()
        tier1_records = tier1_cache.get_unique_records()

        dedup_stats = all_cache.get_stats()
        tier1_dedup_stats = tier1_cache.get_stats()

        print(f"\nDeduplication:")
        print(f"  All: {dedup_stats['total_isoforms']} isoforms -> {dedup_stats['unique_sequences']} unique sequences")
        print(f"       ({dedup_stats['duplicates_skipped']} duplicates removed)")
        print(f"  Tier1: {tier1_dedup_stats['total_isoforms']} isoforms -> {tier1_dedup_stats['unique_sequences']} unique sequences")
        print(f"       ({tier1_dedup_stats['duplicates_skipped']} duplicates removed)")

        # Save isoform mapping for reference
        isoform_map_file = output_dir / 'isoform_map.json'
        with open(isoform_map_file, 'w') as f:
            json.dump({
                'all': all_cache.get_isoform_map(),
                'tier1': tier1_cache.get_isoform_map()
            }, f, indent=2)
        print(f"  Saved isoform mapping to: {isoform_map_file.name}")
    else:
        # No deduplication - flatten all records
        all_records = []
        tier1_records = []
        for gene_id, records in extracted.items():
            all_records.extend(records)
            if gene_id in tier1_genes:
                tier1_records.extend(records)

    # Write sense FASTA files
    sense_file = output_dir / '3UTR_sense.fasta'
    SeqIO.write(all_records, sense_file, 'fasta')
    print(f"\nWrote {len(all_records)} sequences to: {sense_file.name}")

    tier1_sense_file = output_dir / '3UTR_sense_tier1.fasta'
    SeqIO.write(tier1_records, tier1_sense_file, 'fasta')
    print(f"Wrote {len(tier1_records)} tier1 sequences to: {tier1_sense_file.name}")

    # Create antisense (reverse complement) versions
    antisense_records = [reverse_complement(r) for r in all_records]
    tier1_antisense_records = [reverse_complement(r) for r in tier1_records]

    antisense_file = output_dir / '3UTR_antisense.fasta'
    SeqIO.write(antisense_records, antisense_file, 'fasta')
    print(f"Wrote {len(antisense_records)} sequences to: {antisense_file.name}")

    tier1_antisense_file = output_dir / '3UTR_antisense_tier1.fasta'
    SeqIO.write(tier1_antisense_records, tier1_antisense_file, 'fasta')
    print(f"Wrote {len(tier1_antisense_records)} tier1 sequences to: {tier1_antisense_file.name}")

    # Summary statistics
    print("\n" + "=" * 60)
    print("Extraction Summary by Gene:")
    print("-" * 60)
    print(f"{'Gene':<12} {'FBgn ID':<15} {'Tier':<6} {'Isoforms':<10} {'Total bp'}")
    print("-" * 60)

    for gene_id in sorted(extracted.keys(), key=lambda x: genes[x]['symbol']):
        records = extracted[gene_id]
        symbol = genes[gene_id]['symbol']
        tier = genes[gene_id]['tier']
        total_bp = sum(len(r.seq) for r in records)
        print(f"{symbol:<12} {gene_id:<15} {tier:<6} {len(records):<10} {total_bp:,}")

    return {
        'total_genes': len(extracted),
        'total_sequences': total_seqs,
        'unique_sequences': len(all_records),
        'tier1_unique': len(tier1_records),
        'missing_genes': list(missing),
        'duplicates_removed': total_seqs - len(all_records) if deduplicate else 0
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--reference',
        type=Path,
        default=Path('data/references/dmel_3utr.fasta'),
        help='FlyBase 3\'UTR reference FASTA (default: data/references/dmel_3utr.fasta)'
    )
    parser.add_argument(
        '--gene-list',
        type=Path,
        default=Path('data/gene_lists/germ_plasm_genes_consolidated.tsv'),
        help='Consolidated gene list TSV (default: data/gene_lists/germ_plasm_genes_consolidated.tsv)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/queries/germ_plasm'),
        help='Output directory (default: data/queries/germ_plasm)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--no-deduplicate',
        action='store_true',
        help='Disable deduplication (keep all isoforms even if identical)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.reference.exists():
        print(f"Error: Reference file not found: {args.reference}", file=sys.stderr)
        return 1

    if not args.gene_list.exists():
        print(f"Error: Gene list not found: {args.gene_list}", file=sys.stderr)
        print("Run build_germ_plasm_genelist.py first", file=sys.stderr)
        return 1

    print("Extracting Germ Plasm 3'UTRs")
    print("=" * 60)

    stats = extract_3utrs(
        args.reference,
        args.gene_list,
        args.output_dir,
        args.verbose,
        deduplicate=not args.no_deduplicate
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  Genes with 3'UTRs: {stats['total_genes']}")
    print(f"  Total isoforms found: {stats['total_sequences']}")
    print(f"  Unique sequences (after dedup): {stats['unique_sequences']}")
    print(f"  Tier1 unique sequences: {stats['tier1_unique']}")
    if stats['duplicates_removed'] > 0:
        print(f"  Duplicates removed: {stats['duplicates_removed']}")
    if stats['missing_genes']:
        print(f"  Missing genes: {len(stats['missing_genes'])}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
