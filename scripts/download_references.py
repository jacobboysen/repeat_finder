#!/usr/bin/env python3
"""
Download reference files from FlyBase and Dfam for Drosophila melanogaster.

This script downloads genome sequences, annotations, genic segments, gene groups,
and transposable element sequences to data/references/.
"""

import argparse
import gzip
import os
import shutil
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests
from tqdm import tqdm


# FlyBase URLs (AWS S3 - using current release)
# FlyBase moved to AWS S3: https://s3ftp.flybase.org/
FLYBASE_RELEASE = "r6.66"
FLYBASE_FASTA_BASE = f"https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/fasta"
FLYBASE_GFF_BASE = f"https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/current/gff"
FLYBASE_PRECOMPUTED_BASE = f"https://s3ftp.flybase.org/releases/current/precomputed_files/genes"

FLYBASE_URLS = {
    "genome": f"{FLYBASE_FASTA_BASE}/dmel-all-chromosome-{FLYBASE_RELEASE}.fasta.gz",
    "gff": f"{FLYBASE_GFF_BASE}/dmel-all-{FLYBASE_RELEASE}.gff.gz",
    "cds": f"{FLYBASE_FASTA_BASE}/dmel-all-CDS-{FLYBASE_RELEASE}.fasta.gz",
    "5utr": f"{FLYBASE_FASTA_BASE}/dmel-all-five_prime_UTR-{FLYBASE_RELEASE}.fasta.gz",
    "3utr": f"{FLYBASE_FASTA_BASE}/dmel-all-three_prime_UTR-{FLYBASE_RELEASE}.fasta.gz",
    "intron": f"{FLYBASE_FASTA_BASE}/dmel-all-intron-{FLYBASE_RELEASE}.fasta.gz",
    "transposons": f"{FLYBASE_FASTA_BASE}/dmel-all-transposon-{FLYBASE_RELEASE}.fasta.gz",
    "gene_groups": f"{FLYBASE_PRECOMPUTED_BASE}/fbgn_fbtr_fbpp_fb_2025_05.tsv.gz",
    "gene_group_data": f"{FLYBASE_PRECOMPUTED_BASE}/gene_group_data_fb_2025_05.tsv.gz",
}

# Dfam URLs (Updated for Dfam 3.8)
DFAM_URLS = {
    "te_dfam": "https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_drosophila.h5.gz",
}


def download_file(url, output_path, decompress=True, chunk_size=8192):
    """
    Download a file from URL with progress bar.

    Args:
        url: URL to download from
        output_path: Path to save the file
        decompress: If True, decompress .gz files automatically
        chunk_size: Size of chunks to download
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine if we need to decompress
    is_gzipped = url.endswith('.gz')

    # Download to temporary file first
    if is_gzipped and decompress:
        temp_path = output_path.parent / f"{output_path.name}.gz"
    else:
        temp_path = output_path

    print(f"Downloading {url}")
    print(f"  -> {temp_path}")

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with open(temp_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True) as pbar:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

        # Decompress if needed
        if is_gzipped and decompress:
            print(f"Decompressing to {output_path}")
            with gzip.open(temp_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            temp_path.unlink()  # Remove compressed file

        print(f"✓ Downloaded successfully\n")
        return True

    except requests.exceptions.RequestException as e:
        print(f"✗ Error downloading {url}: {e}\n", file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        return False
    except Exception as e:
        print(f"✗ Error processing {url}: {e}\n", file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        return False


def download_dfam_sequences(output_dir):
    """
    Download and extract D. melanogaster sequences from Dfam.

    Note: Dfam 3.8+ provides data in HDF5 format (FamDB), not simple FASTA.
    For basic TE analysis, the FlyBase transposon sequences are usually sufficient.
    If you need Dfam data, consider using the FamDB tools or RepeatMasker.
    """
    print("Note: Dfam TE sequences are in HDF5 format and require special tools.")
    print("FlyBase transposon sequences (already downloaded) are sufficient for most analyses.")
    print("If you need Dfam data, visit: https://www.dfam.org/releases/Dfam_3.8/")
    print("Skipping Dfam download.\n")

    # Return success since we're intentionally skipping
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data/references'),
        help='Output directory for downloaded files (default: data/references)'
    )
    parser.add_argument(
        '--files',
        nargs='+',
        choices=list(FLYBASE_URLS.keys()) + ['dfam', 'all'],
        default=['all'],
        help='Specific files to download (default: all)'
    )
    parser.add_argument(
        '--no-decompress',
        action='store_true',
        help='Keep files compressed (do not decompress .gz files)'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List available files and exit'
    )

    args = parser.parse_args()

    if args.list:
        print("Available files from FlyBase:")
        for key, url in FLYBASE_URLS.items():
            print(f"  {key:20s} {url}")
        print("\nAvailable files from Dfam:")
        print("  dfam                 D. melanogaster TE consensus sequences")
        return 0

    # Determine which files to download
    if 'all' in args.files:
        files_to_download = list(FLYBASE_URLS.keys())
        download_dfam = True
    else:
        files_to_download = [f for f in args.files if f != 'dfam']
        download_dfam = 'dfam' in args.files

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading to: {args.output_dir}")
    print(f"Decompress: {not args.no_decompress}")
    print()

    # Download FlyBase files
    success_count = 0
    fail_count = 0

    for file_key in files_to_download:
        url = FLYBASE_URLS[file_key]

        # Determine output filename
        filename = Path(urlparse(url).path).name
        if not args.no_decompress and filename.endswith('.gz'):
            filename = filename[:-3]  # Remove .gz extension

        # Map to more intuitive names
        name_map = {
            'gene_groups': 'gene_info.tsv',
            'gene_group_data': 'gene_groups.tsv',
            'genome': 'dmel_genome.fasta',
            'gff': 'dmel_annotation.gff',
            'cds': 'dmel_cds.fasta',
            '5utr': 'dmel_5utr.fasta',
            '3utr': 'dmel_3utr.fasta',
            'intron': 'dmel_intron.fasta',
            'transposons': 'dmel_te_flybase.fasta',
        }
        filename = name_map.get(file_key, filename)

        output_path = args.output_dir / filename

        if download_file(url, output_path, decompress=not args.no_decompress):
            success_count += 1
        else:
            fail_count += 1

    # Download Dfam files
    if download_dfam:
        if download_dfam_sequences(args.output_dir):
            success_count += 1
        else:
            fail_count += 1

    # Summary
    print("=" * 60)
    print(f"Download complete: {success_count} succeeded, {fail_count} failed")

    if fail_count > 0:
        print("\n⚠ Some downloads failed. Check errors above.")
        return 1
    else:
        print("\n✓ All downloads completed successfully!")
        return 0


if __name__ == '__main__':
    sys.exit(main())
