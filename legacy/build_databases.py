#!/usr/bin/env python3
"""
Build BLAST databases from reference FASTA files.

Creates nucleotide BLAST databases for:
- D. melanogaster genome
- FlyBase transposons
- Dfam transposons
- Combined (deduplicated) transposon database
"""

import argparse
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO


def run_makeblastdb(input_fasta, output_db, dbtype='nucl', title=None):
    """
    Run makeblastdb to create a BLAST database.

    Args:
        input_fasta: Path to input FASTA file
        output_db: Path to output database (without extension)
        dbtype: Database type ('nucl' or 'prot')
        title: Database title (optional)

    Returns:
        True if successful, False otherwise
    """
    input_fasta = Path(input_fasta)
    output_db = Path(output_db)

    if not input_fasta.exists():
        print(f"✗ Error: Input file not found: {input_fasta}", file=sys.stderr)
        return False

    output_db.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        'makeblastdb',
        '-in', str(input_fasta),
        '-dbtype', dbtype,
        '-out', str(output_db),
        '-parse_seqids'
    ]

    if title:
        cmd.extend(['-title', title])

    print(f"Building database: {output_db.name}")
    print(f"  Input: {input_fasta}")
    print(f"  Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        if result.stdout:
            print(result.stdout)

        # Get database statistics
        stats = get_db_stats(output_db)
        if stats:
            print(f"  Sequences: {stats['num_seqs']}")
            print(f"  Total length: {stats['total_length']:,} bp")

        print("✓ Database created successfully\n")
        return True

    except subprocess.CalledProcessError as e:
        print(f"✗ Error running makeblastdb:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return False
    except FileNotFoundError:
        print("✗ Error: makeblastdb not found. Is BLAST+ installed?", file=sys.stderr)
        print("  Install with: conda install -c bioconda blast", file=sys.stderr)
        return False


def get_db_stats(db_path):
    """
    Get statistics about a BLAST database.

    Args:
        db_path: Path to database (without extension)

    Returns:
        Dictionary with num_seqs and total_length, or None if error
    """
    try:
        result = subprocess.run(
            ['blastdbcmd', '-db', str(db_path), '-info'],
            capture_output=True,
            text=True,
            check=True
        )

        # Parse output for statistics
        stats = {}
        for line in result.stdout.split('\n'):
            if 'sequences;' in line:
                parts = line.split()
                stats['num_seqs'] = int(parts[0].replace(',', ''))
                stats['total_length'] = int(parts[2].replace(',', ''))
                break

        return stats if stats else None

    except Exception:
        return None


def merge_fasta_files(input_files, output_file, deduplicate=True):
    """
    Merge multiple FASTA files, optionally deduplicating by sequence ID.

    Args:
        input_files: List of input FASTA file paths
        output_file: Output FASTA file path
        deduplicate: If True, keep only first occurrence of each sequence ID

    Returns:
        Number of sequences written
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    seen_ids = set()
    seq_count = 0

    print(f"Merging FASTA files to {output_file.name}")

    with open(output_file, 'w') as out_f:
        for fasta_file in input_files:
            fasta_file = Path(fasta_file)

            if not fasta_file.exists():
                print(f"  Warning: {fasta_file} not found, skipping")
                continue

            print(f"  Processing: {fasta_file.name}")
            file_count = 0

            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id

                if deduplicate and seq_id in seen_ids:
                    continue

                SeqIO.write(record, out_f, 'fasta')
                seen_ids.add(seq_id)
                seq_count += 1
                file_count += 1

            print(f"    Added: {file_count} sequences")

    print(f"  Total: {seq_count} sequences written")
    if deduplicate:
        print(f"  Duplicates removed: {len(seen_ids) - seq_count}")

    return seq_count


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--reference-dir',
        type=Path,
        default=Path('data/references'),
        help='Directory containing reference FASTA files (default: data/references)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/blastdb'),
        help='Output directory for BLAST databases (default: data/blastdb)'
    )
    parser.add_argument(
        '--databases',
        nargs='+',
        choices=['genome', 'te_flybase', 'te_dfam', 'te_combined', 'all'],
        default=['all'],
        help='Databases to build (default: all)'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Rebuild databases even if they exist'
    )

    args = parser.parse_args()

    # Determine which databases to build
    if 'all' in args.databases:
        databases_to_build = ['genome', 'te_flybase', 'te_dfam', 'te_combined']
    else:
        databases_to_build = args.databases

    # Database definitions
    db_config = {
        'genome': {
            'input': args.reference_dir / 'dmel_genome.fasta',
            'output': args.output_dir / 'dmel_genome',
            'title': 'D. melanogaster genome'
        },
        'te_flybase': {
            'input': args.reference_dir / 'dmel_te_flybase.fasta',
            'output': args.output_dir / 'dmel_te_flybase',
            'title': 'D. melanogaster TEs (FlyBase)'
        },
        'te_dfam': {
            'input': args.reference_dir / 'dmel_te_dfam.fasta',
            'output': args.output_dir / 'dmel_te_dfam',
            'title': 'D. melanogaster TEs (Dfam)'
        }
    }

    success_count = 0
    fail_count = 0

    # Build individual databases
    for db_name in databases_to_build:
        if db_name == 'te_combined':
            continue  # Handle separately below

        config = db_config[db_name]

        # Check if database exists
        if not args.force and (config['output'].with_suffix('.nhr').exists()):
            print(f"Database {db_name} already exists (use --force to rebuild)\n")
            success_count += 1
            continue

        # Build database
        if run_makeblastdb(
            config['input'],
            config['output'],
            title=config['title']
        ):
            success_count += 1
        else:
            fail_count += 1

    # Build combined TE database if requested
    if 'te_combined' in databases_to_build:
        combined_output = args.output_dir / 'dmel_te_combined'

        if not args.force and combined_output.with_suffix('.nhr').exists():
            print("Database te_combined already exists (use --force to rebuild)\n")
            success_count += 1
        else:
            # Merge FlyBase and Dfam TEs
            temp_fasta = args.output_dir / 'dmel_te_combined.fasta'

            te_files = [
                args.reference_dir / 'dmel_te_flybase.fasta',
                args.reference_dir / 'dmel_te_dfam.fasta'
            ]

            if merge_fasta_files(te_files, temp_fasta, deduplicate=True) > 0:
                if run_makeblastdb(
                    temp_fasta,
                    combined_output,
                    title='D. melanogaster TEs (Combined)'
                ):
                    success_count += 1
                else:
                    fail_count += 1
            else:
                print("✗ Failed to merge TE files\n", file=sys.stderr)
                fail_count += 1

    # Summary
    print("=" * 60)
    print(f"Database build complete: {success_count} succeeded, {fail_count} failed")

    if fail_count > 0:
        print("\n⚠ Some databases failed to build. Check errors above.")
        return 1
    else:
        print("\n✓ All databases built successfully!")
        return 0


if __name__ == '__main__':
    sys.exit(main())
