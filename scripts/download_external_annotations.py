#!/usr/bin/env python3
"""
Download external annotation files from FlyBase and FlyFISH.

Downloads:
- FlyBase RNA-Seq expression data (tissue-specific RPKM)
- FlyBase GO annotations (Gene Ontology associations)
- FlyBase gene group data (pathway memberships)
- FlyFISH RNA localization patterns

Output: data/annotations/raw/
"""

import argparse
import gzip
import shutil
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests
from tqdm import tqdm

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root

# FlyBase precomputed file URLs (AWS S3)
# Using fb_2025_05 release (adjust if needed)
FLYBASE_BASE = "https://s3ftp.flybase.org/releases/current/precomputed_files"

ANNOTATION_URLS = {
    # RNA-Seq expression data
    "rnaseq_rpkm": {
        "url": f"{FLYBASE_BASE}/genes/gene_rpkm_report_fb_2025_05.tsv.gz",
        "output": "gene_rpkm_report.tsv",
        "description": "RNA-Seq RPKM expression values by tissue/stage",
    },
    "rnaseq_matrix": {
        "url": f"{FLYBASE_BASE}/genes/gene_rpkm_matrix_fb_2025_05.tsv.gz",
        "output": "gene_rpkm_matrix.tsv",
        "description": "RNA-Seq RPKM matrix (genes x samples)",
    },
    # Gene Ontology annotations
    "go_annotations": {
        "url": f"{FLYBASE_BASE}/go/gene_association.fb.gz",
        "output": "gene_association.gaf",
        "description": "GO annotations in GAF 2.2 format",
    },
    # Gene groups and pathways
    "gene_groups": {
        "url": f"{FLYBASE_BASE}/genes/gene_group_data_fb_2025_05.tsv.gz",
        "output": "gene_group_data.tsv",
        "description": "Gene group memberships (pathways, families)",
    },
    "pathway_groups": {
        "url": f"{FLYBASE_BASE}/genes/pathway_group_data_fb_2025_05.tsv.gz",
        "output": "pathway_group_data.tsv",
        "description": "Signaling pathway group data",
    },
    # FlyFISH RNA localization
    "flyfish": {
        "url": "https://fly-fish.ccbr.utoronto.ca/probes.csv",
        "output": "flyfish_localization.csv",
        "description": "FlyFISH RNA subcellular localization patterns",
    },
}


def download_file(url, output_path, decompress=True, chunk_size=8192, timeout=60):
    """
    Download a file from URL with progress bar.

    Args:
        url: URL to download from
        output_path: Path to save the file
        decompress: If True, decompress .gz files automatically
        chunk_size: Size of chunks to download
        timeout: Request timeout in seconds

    Returns:
        True if successful, False otherwise
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    is_gzipped = url.endswith('.gz')

    if is_gzipped and decompress:
        temp_path = output_path.parent / f"{output_path.name}.gz"
    else:
        temp_path = output_path

    print(f"Downloading: {url}")
    print(f"  -> {output_path}")

    try:
        response = requests.get(url, stream=True, timeout=timeout)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with open(temp_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, disable=total_size == 0) as pbar:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

        # Decompress if needed
        if is_gzipped and decompress:
            print(f"  Decompressing...")
            with gzip.open(temp_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            temp_path.unlink()

        # Verify file has content
        if output_path.stat().st_size == 0:
            print(f"  Warning: Downloaded file is empty")
            return False

        print(f"  Done ({output_path.stat().st_size:,} bytes)\n")
        return True

    except requests.exceptions.RequestException as e:
        print(f"  Error: {e}\n", file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        return False
    except Exception as e:
        print(f"  Error: {e}\n", file=sys.stderr)
        if temp_path.exists():
            temp_path.unlink()
        return False


def download_flyfish(output_path):
    """
    Download FlyFISH data from probes.csv endpoint.
    """
    url = ANNOTATION_URLS["flyfish"]["url"]
    if download_file(url, output_path, decompress=False, timeout=120):
        return True

    print("  FlyFISH download failed. You may need to download manually from:")
    print("  https://fly-fish.ccbr.utoronto.ca/")
    return False


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Output directory (default: data/annotations/raw)'
    )
    parser.add_argument(
        '--files',
        nargs='+',
        choices=list(ANNOTATION_URLS.keys()) + ['all', 'flybase', 'expression', 'go'],
        default=['all'],
        help='Specific files to download (default: all)'
    )
    parser.add_argument(
        '--list',
        action='store_true',
        help='List available files and exit'
    )
    parser.add_argument(
        '--skip-flyfish',
        action='store_true',
        help='Skip FlyFISH download (may be slow/unavailable)'
    )

    args = parser.parse_args()

    if args.list:
        print("Available annotation files:\n")
        for key, info in ANNOTATION_URLS.items():
            print(f"  {key:20s} {info['description']}")
            print(f"  {' ':20s} -> {info['output']}")
        print("\nGroups:")
        print("  all         Download all files")
        print("  flybase     FlyBase files only (excludes FlyFISH)")
        print("  expression  RNA-Seq expression data only")
        print("  go          GO annotations only")
        return 0

    # Set output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = get_project_root() / "data" / "annotations" / "raw"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine which files to download
    files_to_download = set()

    for file_spec in args.files:
        if file_spec == 'all':
            files_to_download.update(ANNOTATION_URLS.keys())
        elif file_spec == 'flybase':
            files_to_download.update(k for k in ANNOTATION_URLS.keys() if k != 'flyfish')
        elif file_spec == 'expression':
            files_to_download.update(['rnaseq_rpkm', 'rnaseq_matrix'])
        elif file_spec == 'go':
            files_to_download.add('go_annotations')
        else:
            files_to_download.add(file_spec)

    if args.skip_flyfish:
        files_to_download.discard('flyfish')

    print(f"Output directory: {output_dir}")
    print(f"Files to download: {', '.join(sorted(files_to_download))}\n")

    # Download files
    success_count = 0
    fail_count = 0

    for file_key in sorted(files_to_download):
        info = ANNOTATION_URLS[file_key]
        output_path = output_dir / info['output']

        print(f"[{file_key}] {info['description']}")

        if file_key == 'flyfish':
            success = download_flyfish(output_path)
        else:
            success = download_file(info['url'], output_path)

        if success:
            success_count += 1
        else:
            fail_count += 1

    # Summary
    print("=" * 60)
    print(f"Download complete: {success_count} succeeded, {fail_count} failed")

    if fail_count > 0:
        print("\nSome downloads failed. You may need to:")
        print("1. Check your internet connection")
        print("2. Try again later (FlyBase/FlyFISH servers may be busy)")
        print("3. Download manually from the URLs listed above")
        return 1
    else:
        print(f"\nAll files saved to: {output_dir}")
        return 0


if __name__ == '__main__':
    sys.exit(main())
