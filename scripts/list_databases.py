#!/usr/bin/env python3
"""
List available BLAST databases with statistics.

Shows database name, type, number of sequences, and total length.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def get_blast_databases(db_dir):
    """
    Find all BLAST databases in directory.

    Args:
        db_dir: Directory to search for databases

    Returns:
        List of database paths (without extensions)
    """
    db_dir = Path(db_dir)

    if not db_dir.exists():
        return []

    # Find .nhr files (nucleotide) or .phr files (protein)
    databases = []

    for ext in ['.nhr', '.phr']:
        for db_file in db_dir.glob(f"*{ext}"):
            db_path = db_file.with_suffix('')  # Remove extension
            if db_path not in databases:
                databases.append(db_path)

    return sorted(databases)


def get_database_info(db_path):
    """
    Get information about a BLAST database.

    Args:
        db_path: Path to database (without extension)

    Returns:
        Dictionary with database information, or None if error
    """
    try:
        result = subprocess.run(
            ['blastdbcmd', '-db', str(db_path), '-info'],
            capture_output=True,
            text=True,
            check=True,
            timeout=10
        )

        info = {
            'name': db_path.name,
            'path': str(db_path),
            'type': None,
            'title': None,
            'num_seqs': None,
            'total_length': None,
            'date': None
        }

        # Parse output
        for line in result.stdout.split('\n'):
            line = line.strip()

            if line.startswith('Database:'):
                info['title'] = line.replace('Database:', '').strip()

            elif 'protein' in line.lower():
                info['type'] = 'protein'

            elif 'nucleotide' in line.lower():
                info['type'] = 'nucleotide'

            elif 'sequences;' in line:
                # Format: "X sequences; Y total bases/residues"
                parts = line.split()
                if len(parts) >= 4:
                    info['num_seqs'] = int(parts[0].replace(',', ''))
                    info['total_length'] = int(parts[2].replace(',', ''))

            elif line.startswith('Date:'):
                info['date'] = line.replace('Date:', '').strip()

        # Infer type from file extensions if not found
        if info['type'] is None:
            if db_path.with_suffix('.nhr').exists():
                info['type'] = 'nucleotide'
            elif db_path.with_suffix('.phr').exists():
                info['type'] = 'protein'

        return info

    except subprocess.CalledProcessError:
        return None
    except FileNotFoundError:
        print("Error: blastdbcmd not found. Is BLAST+ installed?", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error getting info for {db_path}: {e}", file=sys.stderr)
        return None


def format_size(num_bases):
    """Format number of bases with appropriate units."""
    if num_bases is None:
        return "N/A"

    if num_bases < 1000:
        return f"{num_bases} bp"
    elif num_bases < 1_000_000:
        return f"{num_bases / 1000:.1f} Kb"
    elif num_bases < 1_000_000_000:
        return f"{num_bases / 1_000_000:.1f} Mb"
    else:
        return f"{num_bases / 1_000_000_000:.1f} Gb"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--db-dir',
        type=Path,
        default=Path('data/blastdb'),
        help='Directory containing BLAST databases (default: data/blastdb)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show additional information'
    )

    args = parser.parse_args()

    if not args.db_dir.exists():
        print(f"Error: Database directory not found: {args.db_dir}", file=sys.stderr)
        print("Run build_databases.py to create databases", file=sys.stderr)
        return 1

    # Find databases
    databases = get_blast_databases(args.db_dir)

    if not databases:
        print(f"No BLAST databases found in {args.db_dir}")
        print("Run build_databases.py to create databases")
        return 0

    print(f"BLAST Databases in {args.db_dir}")
    print("=" * 80)
    print()

    # Get info for each database
    for db_path in databases:
        info = get_database_info(db_path)

        if not info:
            print(f"✗ {db_path.name}: Unable to read database info")
            continue

        # Print basic info
        print(f"Database: {info['name']}")

        if args.verbose and info['title']:
            print(f"  Title: {info['title']}")

        print(f"  Type: {info['type'] or 'unknown'}")
        print(f"  Sequences: {info['num_seqs']:,}" if info['num_seqs'] else "  Sequences: N/A")
        print(f"  Size: {format_size(info['total_length'])}")

        if args.verbose:
            print(f"  Path: {info['path']}")
            if info['date']:
                print(f"  Created: {info['date']}")

        print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
