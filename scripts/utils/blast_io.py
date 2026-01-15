"""
BLAST I/O utilities for the TE fossil mining pipeline.

Provides consistent parsing of BLAST tabular output across all analysis scripts.
All TSV files use the 17-column format documented in FILE_MAP.md.
"""

from pathlib import Path
from typing import Union

import pandas as pd


# Standard BLAST output columns (17-column format)
BLAST_COLUMNS = [
    'qseqid',    # Query sequence ID
    'sseqid',    # Subject (TE) sequence ID
    'pident',    # Percent identity
    'length',    # Alignment length
    'mismatch',  # Number of mismatches
    'gapopen',   # Number of gap openings
    'qstart',    # Query start position
    'qend',      # Query end position
    'sstart',    # Subject start position
    'send',      # Subject end position
    'evalue',    # E-value
    'bitscore',  # Bit score
    'qlen',      # Query length
    'slen',      # Subject length
    'qseq',      # Query sequence (aligned portion)
    'sseq',      # Subject sequence (aligned portion)
    'strand',    # Strand orientation (plus/minus)
]

# 16-column format (without pre-computed strand)
BLAST_COLUMNS_NO_STRAND = BLAST_COLUMNS[:-1]


def classify_strand(sstart: int, send: int) -> str:
    """
    Classify hit strand based on subject coordinates.

    In BLAST output, when sstart > send, the hit is on the minus strand.

    Args:
        sstart: Subject start position
        send: Subject end position

    Returns:
        'plus' or 'minus'
    """
    return 'plus' if sstart < send else 'minus'


def load_blast_results(
    results_file: Union[str, Path],
    add_strand: bool = True,
    min_length: int = None,
    min_pident: float = None,
    max_evalue: float = None,
) -> pd.DataFrame:
    """
    Load BLAST results from TSV file.

    Automatically detects whether the file has 16 or 17 columns and adds
    strand classification if needed.

    Args:
        results_file: Path to BLAST TSV file
        add_strand: If True, add strand column if not present
        min_length: Minimum alignment length filter
        min_pident: Minimum percent identity filter
        max_evalue: Maximum E-value filter

    Returns:
        DataFrame with BLAST results
    """
    results_file = Path(results_file)

    if not results_file.exists():
        return pd.DataFrame(columns=BLAST_COLUMNS)

    if results_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    # Detect column count from first line
    with open(results_file) as f:
        first_line = f.readline()
    num_cols = len(first_line.strip().split('\t'))

    # Load with appropriate columns
    if num_cols == 17:
        df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS)
    elif num_cols == 16:
        df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS_NO_STRAND)
        if add_strand:
            df['strand'] = df.apply(
                lambda row: classify_strand(row['sstart'], row['send']),
                axis=1
            )
    else:
        # Try to load with basic columns
        basic_cols = BLAST_COLUMNS[:min(num_cols, len(BLAST_COLUMNS))]
        df = pd.read_csv(results_file, sep='\t', names=basic_cols)
        if add_strand and 'sstart' in df.columns and 'send' in df.columns:
            df['strand'] = df.apply(
                lambda row: classify_strand(row['sstart'], row['send']),
                axis=1
            )

    # Apply filters
    if min_length is not None and 'length' in df.columns:
        df = df[df['length'] >= min_length]

    if min_pident is not None and 'pident' in df.columns:
        df = df[df['pident'] >= min_pident]

    if max_evalue is not None and 'evalue' in df.columns:
        df = df[df['evalue'] <= max_evalue]

    return df


def parse_blast_line(line: str) -> dict:
    """
    Parse a single BLAST TSV line into a dictionary.

    Useful for streaming large files without loading into memory.

    Args:
        line: Tab-separated BLAST output line

    Returns:
        Dictionary with parsed values
    """
    parts = line.strip().split('\t')

    result = {}
    for i, col in enumerate(BLAST_COLUMNS):
        if i >= len(parts):
            result[col] = None
            continue

        value = parts[i]

        # Type conversion based on column
        if col in ('pident', 'evalue', 'bitscore'):
            result[col] = float(value) if value else 0.0
        elif col in ('length', 'mismatch', 'gapopen', 'qstart', 'qend',
                     'sstart', 'send', 'qlen', 'slen'):
            result[col] = int(value) if value else 0
        else:
            result[col] = value

    # Add strand if not present
    if 'strand' not in result or result['strand'] is None:
        if result.get('sstart') and result.get('send'):
            result['strand'] = classify_strand(result['sstart'], result['send'])

    return result


def iter_blast_results(results_file: Union[str, Path]):
    """
    Iterate over BLAST results line by line.

    Useful for processing large files without loading into memory.

    Args:
        results_file: Path to BLAST TSV file

    Yields:
        Dictionary for each BLAST hit
    """
    results_file = Path(results_file)

    if not results_file.exists():
        return

    with open(results_file) as f:
        for line in f:
            if line.strip():
                yield parse_blast_line(line)


def summarize_blast_results(df: pd.DataFrame) -> dict:
    """
    Generate summary statistics for BLAST results.

    Args:
        df: DataFrame with BLAST results

    Returns:
        Dictionary with summary statistics
    """
    if df.empty:
        return {
            'total_hits': 0,
            'unique_queries': 0,
            'unique_subjects': 0,
            'mean_pident': 0,
            'mean_length': 0,
            'mean_evalue': 0,
            'strand_plus': 0,
            'strand_minus': 0,
        }

    summary = {
        'total_hits': len(df),
        'unique_queries': df['qseqid'].nunique() if 'qseqid' in df.columns else 0,
        'unique_subjects': df['sseqid'].nunique() if 'sseqid' in df.columns else 0,
        'mean_pident': df['pident'].mean() if 'pident' in df.columns else 0,
        'mean_length': df['length'].mean() if 'length' in df.columns else 0,
        'mean_evalue': df['evalue'].mean() if 'evalue' in df.columns else 0,
    }

    if 'strand' in df.columns:
        strand_counts = df['strand'].value_counts()
        summary['strand_plus'] = strand_counts.get('plus', 0)
        summary['strand_minus'] = strand_counts.get('minus', 0)
    else:
        summary['strand_plus'] = 0
        summary['strand_minus'] = 0

    return summary
