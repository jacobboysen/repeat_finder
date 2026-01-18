"""
Disorder prediction loading and utilities.

Provides functions for loading pre-computed disorder predictions
and integrating them into TE analyses.

Disorder predictions are expected as TSV files with columns:
- fbgn: FlyBase gene ID
- disorder_fraction: Fraction of residues with disorder score > 0.5
- max_disordered_region: Length of longest contiguous disordered stretch
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


def load_disorder_scores(path: Union[str, Path]) -> Dict[str, dict]:
    """
    Load disorder predictions from TSV file.

    Args:
        path: Path to disorder scores TSV

    Returns:
        Dictionary mapping FBgn to disorder info dict with keys:
            - disorder_fraction: float (0-1)
            - max_disordered_region: int (amino acids)
            - total_length: int (protein length)
    """
    path = Path(path)
    scores = {}

    with open(path) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            row = dict(zip(header, parts))

            fbgn = row.get('fbgn', '')
            if not fbgn.startswith('FBgn'):
                continue

            scores[fbgn] = {
                'disorder_fraction': float(row.get('disorder_fraction', 0)),
                'max_disordered_region': int(row.get('max_disordered_region', 0)),
                'total_length': int(row.get('total_length', 0)),
            }

    return scores


def get_disorder_category(fraction: float) -> str:
    """
    Categorize disorder level.

    Args:
        fraction: Fraction of protein that is disordered (0-1)

    Returns:
        Category string: 'low', 'medium', or 'high'
    """
    if fraction < 0.2:
        return 'low'
    elif fraction < 0.5:
        return 'medium'
    else:
        return 'high'


def get_disorder_percentile(fraction: float, all_fractions: List[float]) -> float:
    """
    Calculate percentile rank of a disorder fraction.

    Args:
        fraction: Disorder fraction to rank
        all_fractions: List of all disorder fractions for comparison

    Returns:
        Percentile (0-100)
    """
    if not all_fractions:
        return 50.0

    n_below = sum(1 for f in all_fractions if f < fraction)
    return 100 * n_below / len(all_fractions)


def categorize_genes_by_disorder(
    scores: Dict[str, dict],
    low_threshold: float = 0.2,
    high_threshold: float = 0.5
) -> Dict[str, List[str]]:
    """
    Group genes by disorder level.

    Args:
        scores: Disorder scores dictionary
        low_threshold: Max fraction for 'low' category
        high_threshold: Min fraction for 'high' category

    Returns:
        Dictionary with keys 'low', 'medium', 'high' mapping to gene lists
    """
    categories = {
        'low': [],
        'medium': [],
        'high': []
    }

    for fbgn, info in scores.items():
        fraction = info.get('disorder_fraction', 0)

        if fraction < low_threshold:
            categories['low'].append(fbgn)
        elif fraction >= high_threshold:
            categories['high'].append(fbgn)
        else:
            categories['medium'].append(fbgn)

    return categories


def calculate_disorder_stats(scores: Dict[str, dict]) -> dict:
    """
    Calculate summary statistics for disorder predictions.

    Args:
        scores: Disorder scores dictionary

    Returns:
        Dictionary with summary stats
    """
    if not scores:
        return {
            'n_genes': 0,
            'mean_disorder': 0,
            'median_disorder': 0,
            'n_low': 0,
            'n_medium': 0,
            'n_high': 0
        }

    fractions = [s['disorder_fraction'] for s in scores.values()]
    fractions.sort()

    n = len(fractions)
    median = fractions[n // 2] if n % 2 == 1 else (fractions[n//2 - 1] + fractions[n//2]) / 2

    categories = categorize_genes_by_disorder(scores)

    return {
        'n_genes': n,
        'mean_disorder': sum(fractions) / n,
        'median_disorder': median,
        'min_disorder': min(fractions),
        'max_disorder': max(fractions),
        'n_low': len(categories['low']),
        'n_medium': len(categories['medium']),
        'n_high': len(categories['high']),
    }


def merge_disorder_with_te_data(
    te_data: Dict[str, dict],
    disorder_scores: Dict[str, dict]
) -> Dict[str, dict]:
    """
    Merge disorder predictions with TE analysis data.

    Args:
        te_data: Dictionary of TE statistics per gene (fbgn -> stats)
        disorder_scores: Dictionary of disorder scores per gene

    Returns:
        Merged dictionary with disorder info added to each gene
    """
    merged = {}

    for fbgn, te_stats in te_data.items():
        merged[fbgn] = te_stats.copy()

        if fbgn in disorder_scores:
            disorder_info = disorder_scores[fbgn]
            merged[fbgn]['disorder_fraction'] = disorder_info['disorder_fraction']
            merged[fbgn]['disorder_category'] = get_disorder_category(disorder_info['disorder_fraction'])
            merged[fbgn]['max_disordered_region'] = disorder_info.get('max_disordered_region', 0)
        else:
            merged[fbgn]['disorder_fraction'] = None
            merged[fbgn]['disorder_category'] = None
            merged[fbgn]['max_disordered_region'] = None

    return merged


def correlate_disorder_with_te_density(
    merged_data: Dict[str, dict],
    min_hits: int = 10
) -> Tuple[List[float], List[float], int]:
    """
    Extract disorder fraction and TE density for correlation analysis.

    Args:
        merged_data: Merged TE + disorder data per gene
        min_hits: Minimum TE hits to include gene

    Returns:
        Tuple of (disorder_fractions, te_densities, n_genes)
    """
    disorder_fractions = []
    te_densities = []

    for fbgn, data in merged_data.items():
        # Skip genes without disorder data
        if data.get('disorder_fraction') is None:
            continue

        # Skip genes with too few hits
        if data.get('total_hits', 0) < min_hits:
            continue

        disorder_fractions.append(data['disorder_fraction'])
        te_densities.append(data.get('density', 0))

    return disorder_fractions, te_densities, len(disorder_fractions)
