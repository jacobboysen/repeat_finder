#!/usr/bin/env python3
"""
Summarize BLAST results with hit statistics and distributions.

Analyzes BLAST output files to provide summary statistics including:
- Total hits and unique queries/subjects
- E-value distribution
- Percent identity distribution
- Alignment length distribution
"""

import argparse
import sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_blast_results(blast_file):
    """
    Parse BLAST results file in tabular format.

    Expected columns (outfmt 6):
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

    Returns:
        pandas DataFrame
    """
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]

    try:
        df = pd.read_csv(blast_file, sep='\t', names=columns, comment='#')
        return df
    except Exception as e:
        print(f"Error reading BLAST file: {e}", file=sys.stderr)
        return None


def calculate_statistics(df):
    """
    Calculate summary statistics from BLAST results.

    Returns:
        Dictionary of statistics
    """
    stats = {}

    # Basic counts
    stats['total_hits'] = len(df)
    stats['unique_queries'] = df['qseqid'].nunique()
    stats['unique_subjects'] = df['sseqid'].nunique()

    # Queries with/without hits
    queries_with_hits = df['qseqid'].nunique()
    stats['queries_with_hits'] = queries_with_hits

    # E-value statistics
    stats['evalue_min'] = df['evalue'].min()
    stats['evalue_max'] = df['evalue'].max()
    stats['evalue_median'] = df['evalue'].median()

    # Percent identity statistics
    stats['pident_min'] = df['pident'].min()
    stats['pident_max'] = df['pident'].max()
    stats['pident_mean'] = df['pident'].mean()
    stats['pident_median'] = df['pident'].median()

    # Alignment length statistics
    stats['length_min'] = df['length'].min()
    stats['length_max'] = df['length'].max()
    stats['length_mean'] = df['length'].mean()
    stats['length_median'] = df['length'].median()

    # Bitscore statistics
    stats['bitscore_min'] = df['bitscore'].min()
    stats['bitscore_max'] = df['bitscore'].max()
    stats['bitscore_mean'] = df['bitscore'].mean()
    stats['bitscore_median'] = df['bitscore'].median()

    # Coverage statistics
    df['query_coverage'] = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100
    df['subject_coverage'] = (abs(df['send'] - df['sstart']) + 1) / df['slen'] * 100

    stats['query_cov_mean'] = df['query_coverage'].mean()
    stats['subject_cov_mean'] = df['subject_coverage'].mean()

    return stats


def print_statistics(stats):
    """Print statistics in a formatted table."""
    print("\nBLAST Results Summary")
    print("=" * 80)

    print("\nHit Counts:")
    print(f"  Total hits:         {stats['total_hits']:,}")
    print(f"  Unique queries:     {stats['unique_queries']:,}")
    print(f"  Unique subjects:    {stats['unique_subjects']:,}")
    print(f"  Queries with hits:  {stats['queries_with_hits']:,}")

    print("\nE-value Distribution:")
    print(f"  Minimum:   {stats['evalue_min']:.2e}")
    print(f"  Median:    {stats['evalue_median']:.2e}")
    print(f"  Maximum:   {stats['evalue_max']:.2e}")

    print("\nPercent Identity:")
    print(f"  Minimum:   {stats['pident_min']:.2f}%")
    print(f"  Mean:      {stats['pident_mean']:.2f}%")
    print(f"  Median:    {stats['pident_median']:.2f}%")
    print(f"  Maximum:   {stats['pident_max']:.2f}%")

    print("\nAlignment Length:")
    print(f"  Minimum:   {stats['length_min']:,} bp")
    print(f"  Mean:      {stats['length_mean']:.1f} bp")
    print(f"  Median:    {stats['length_median']:.1f} bp")
    print(f"  Maximum:   {stats['length_max']:,} bp")

    print("\nBit Score:")
    print(f"  Minimum:   {stats['bitscore_min']:.1f}")
    print(f"  Mean:      {stats['bitscore_mean']:.1f}")
    print(f"  Median:    {stats['bitscore_median']:.1f}")
    print(f"  Maximum:   {stats['bitscore_max']:.1f}")

    print("\nCoverage:")
    print(f"  Query coverage (mean):    {stats['query_cov_mean']:.1f}%")
    print(f"  Subject coverage (mean):  {stats['subject_cov_mean']:.1f}%")


def plot_distributions(df, output_dir):
    """
    Create distribution plots for key statistics.

    Args:
        df: BLAST results DataFrame
        output_dir: Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('BLAST Results Distributions', fontsize=16)

    # E-value distribution (log scale)
    ax = axes[0, 0]
    log_evalue = np.log10(df['evalue'].replace(0, 1e-300))
    ax.hist(log_evalue, bins=50, edgecolor='black')
    ax.set_xlabel('Log10(E-value)')
    ax.set_ylabel('Count')
    ax.set_title('E-value Distribution')

    # Percent identity distribution
    ax = axes[0, 1]
    ax.hist(df['pident'], bins=50, edgecolor='black')
    ax.set_xlabel('Percent Identity (%)')
    ax.set_ylabel('Count')
    ax.set_title('Percent Identity Distribution')

    # Alignment length distribution
    ax = axes[1, 0]
    ax.hist(df['length'], bins=50, edgecolor='black')
    ax.set_xlabel('Alignment Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Alignment Length Distribution')

    # Query coverage distribution
    ax = axes[1, 1]
    query_coverage = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100
    ax.hist(query_coverage, bins=50, edgecolor='black')
    ax.set_xlabel('Query Coverage (%)')
    ax.set_ylabel('Count')
    ax.set_title('Query Coverage Distribution')

    plt.tight_layout()

    # Save plot
    plot_file = output_dir / 'distributions.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Distribution plots saved: {plot_file}")

    plt.close()


def get_top_hits(df, n=10):
    """
    Get top queries and subjects by hit count.

    Returns:
        Tuple of (top_queries, top_subjects)
    """
    query_counts = df['qseqid'].value_counts().head(n)
    subject_counts = df['sseqid'].value_counts().head(n)

    return query_counts, subject_counts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'blast_file',
        type=Path,
        help='BLAST results file (tabular format)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Directory for plots and summary files (default: same as input)'
    )
    parser.add_argument(
        '--plot',
        action='store_true',
        help='Generate distribution plots'
    )
    parser.add_argument(
        '--top',
        type=int,
        default=10,
        help='Number of top queries/subjects to show (default: 10)'
    )

    args = parser.parse_args()

    # Check input file
    if not args.blast_file.exists():
        print(f"Error: BLAST file not found: {args.blast_file}", file=sys.stderr)
        return 1

    # Set output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = args.blast_file.parent

    print(f"Analyzing BLAST results: {args.blast_file}")

    # Parse BLAST results
    df = parse_blast_results(args.blast_file)

    if df is None or df.empty:
        print("Error: No BLAST results found or unable to parse file", file=sys.stderr)
        return 1

    # Calculate statistics
    stats = calculate_statistics(df)

    # Print statistics
    print_statistics(stats)

    # Get top hits
    print(f"\nTop {args.top} Queries by Hit Count:")
    print("-" * 60)
    top_queries, top_subjects = get_top_hits(df, args.top)

    for query, count in top_queries.items():
        print(f"  {count:5d} hits  {query}")

    print(f"\nTop {args.top} Subjects by Hit Count:")
    print("-" * 60)
    for subject, count in top_subjects.items():
        # Truncate long subject names
        display_name = subject
        if len(display_name) > 60:
            display_name = display_name[:57] + "..."
        print(f"  {count:5d} hits  {display_name}")

    # Generate plots if requested
    if args.plot:
        try:
            plot_distributions(df, output_dir)
        except Exception as e:
            print(f"\nWarning: Could not generate plots: {e}", file=sys.stderr)

    print("\n" + "=" * 80)

    return 0


if __name__ == '__main__':
    sys.exit(main())
