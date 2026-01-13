#!/usr/bin/env python3
"""
Extract 3'UTRs for germ plasm genes from FlyBase reference.

Reads the consolidated gene list and extracts corresponding 3'UTR sequences
from the FlyBase reference FASTA. Handles multiple transcript isoforms per gene.

Outputs:
- 3UTR_sense.fasta: Original strand sequences
- 3UTR_antisense.fasta: Reverse complement (negative control)
- Tier1 subsets of both
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def extract_3utrs(reference_fasta, gene_list, output_dir, verbose=False):
    """Extract 3'UTRs for genes in the list."""
    output_dir.mkdir(parents=True, exist_ok=True)

    genes, tier1_genes = load_gene_list(gene_list)

    if verbose:
        print(f"Loaded {len(genes)} genes from gene list")
        print(f"  Tier 1 genes: {len(tier1_genes)}")

    # Track results by gene
    extracted = defaultdict(list)
    total_seqs = 0

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

    print(f"  Found {total_seqs} 3'UTR sequences for {len(extracted)} genes")

    # Check for missing genes
    missing = set(genes.keys()) - set(extracted.keys())
    if missing:
        print(f"\n  Warning: {len(missing)} genes not found in reference:")
        for fbgn in missing:
            print(f"    - {genes[fbgn]['symbol']} ({fbgn})")

    # Flatten all records
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
        'tier1_sequences': len(tier1_records),
        'missing_genes': list(missing)
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
        args.verbose
    )

    print("\n" + "=" * 60)
    print("Extraction complete!")
    print(f"  Genes with 3'UTRs: {stats['total_genes']}")
    print(f"  Total sequences: {stats['total_sequences']}")
    print(f"  Tier1 sequences: {stats['tier1_sequences']}")
    if stats['missing_genes']:
        print(f"  Missing genes: {len(stats['missing_genes'])}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
