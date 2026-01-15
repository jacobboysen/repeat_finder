#!/usr/bin/env python3
"""
Analyze TE enrichment in functional gene sets vs genome-wide baseline.

Performs:
- Fisher's exact test for TE presence enrichment (binary)
- Mann-Whitney U test for TE density distribution comparison
- Benjamini-Hochberg multiple testing correction

Output: results/functional_te_enrichment.tsv
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.data_loaders import load_gene_list, load_strand_bias_data


def load_gene_sets(gene_sets_dir: Path) -> dict:
    """
    Load all gene sets from directory.

    Args:
        gene_sets_dir: Directory containing *_fbgn_ids.txt files

    Returns:
        Dictionary mapping set_name -> list of FBgn IDs
    """
    gene_sets = {}

    for path in sorted(gene_sets_dir.glob("*_fbgn_ids.txt")):
        set_name = path.stem.replace('_fbgn_ids', '')
        genes = load_gene_list(path)
        if genes:
            gene_sets[set_name] = genes

    return gene_sets


def compute_te_metrics(genes: set, te_data: dict) -> dict:
    """
    Compute TE metrics for a set of genes.

    Args:
        genes: Set of FBgn IDs
        te_data: Dictionary from load_strand_bias_data()

    Returns:
        Dictionary with TE statistics
    """
    total_hits = 0
    sense_hits = 0
    anti_hits = 0
    total_utr_len = 0
    genes_with_hits = 0
    hit_counts = []
    densities = []

    for fbgn in genes:
        data = te_data.get(fbgn)
        if data:
            hits = data['total_hits']
            utr_len = data['total_utr_len']

            total_hits += hits
            sense_hits += data['sense_hits']
            anti_hits += data['anti_hits']
            total_utr_len += utr_len

            if hits > 0:
                genes_with_hits += 1
                hit_counts.append(hits)

                # Density: hits per kb
                if utr_len > 0:
                    density = (hits / utr_len) * 1000
                    densities.append(density)

    n_genes = len(genes)

    return {
        'n_genes': n_genes,
        'n_with_hits': genes_with_hits,
        'pct_with_hits': 100 * genes_with_hits / n_genes if n_genes > 0 else 0,
        'total_hits': total_hits,
        'total_utr_len': total_utr_len,
        'mean_density': np.mean(densities) if densities else 0,
        'median_density': np.median(densities) if densities else 0,
        'sense_pct': 100 * sense_hits / total_hits if total_hits > 0 else 50,
        'anti_pct': 100 * anti_hits / total_hits if total_hits > 0 else 50,
        'densities': densities,  # For statistical tests
    }


def fisher_exact_te_enrichment(
    gene_set: set,
    all_genes: set,
    genes_with_te: set
) -> dict:
    """
    Fisher's exact test for TE enrichment.

    Contingency table:
                    | TE+ | TE- |
    In gene set     |  a  |  b  |
    Not in gene set |  c  |  d  |

    Args:
        gene_set: Set of FBgn IDs in the functional set
        all_genes: Set of all FBgn IDs in genome
        genes_with_te: Set of FBgn IDs with TE hits

    Returns:
        Dictionary with odds ratio, p-value, and contingency table
    """
    in_set_with_te = len(gene_set & genes_with_te)
    in_set_without_te = len(gene_set - genes_with_te)
    not_in_set_with_te = len((all_genes - gene_set) & genes_with_te)
    not_in_set_without_te = len((all_genes - gene_set) - genes_with_te)

    # Contingency table
    table = [
        [in_set_with_te, in_set_without_te],
        [not_in_set_with_te, not_in_set_without_te]
    ]

    odds_ratio, p_value = stats.fisher_exact(table, alternative='two-sided')

    return {
        'odds_ratio': odds_ratio,
        'p_value': p_value,
        'a': in_set_with_te,
        'b': in_set_without_te,
        'c': not_in_set_with_te,
        'd': not_in_set_without_te,
    }


def mann_whitney_te_density(
    set_densities: list,
    background_densities: list
) -> dict:
    """
    Mann-Whitney U test comparing TE density distributions.

    Args:
        set_densities: TE densities for genes in the set
        background_densities: TE densities for genes not in the set

    Returns:
        Dictionary with U statistic and p-value
    """
    if len(set_densities) < 3 or len(background_densities) < 3:
        return {
            'u_statistic': np.nan,
            'p_value': np.nan,
            'median_in': np.median(set_densities) if set_densities else 0,
            'median_out': np.median(background_densities) if background_densities else 0,
        }

    try:
        statistic, p_value = stats.mannwhitneyu(
            set_densities,
            background_densities,
            alternative='two-sided'
        )
    except ValueError:
        return {
            'u_statistic': np.nan,
            'p_value': np.nan,
            'median_in': np.median(set_densities),
            'median_out': np.median(background_densities),
        }

    return {
        'u_statistic': statistic,
        'p_value': p_value,
        'median_in': np.median(set_densities),
        'median_out': np.median(background_densities),
    }


def apply_fdr_correction(p_values: list, method: str = 'fdr_bh') -> list:
    """
    Apply multiple testing correction.

    Args:
        p_values: List of p-values
        method: Correction method ('fdr_bh' for Benjamini-Hochberg)

    Returns:
        List of adjusted p-values (q-values)
    """
    try:
        from statsmodels.stats.multitest import multipletests
        _, q_values, _, _ = multipletests(p_values, method=method)
        return list(q_values)
    except ImportError:
        # Fallback: simple Bonferroni
        n = len(p_values)
        return [min(p * n, 1.0) for p in p_values]


def analyze_enrichment(
    gene_sets: dict,
    te_data: dict,
    all_genes: set
) -> pd.DataFrame:
    """
    Perform enrichment analysis for all gene sets.

    Args:
        gene_sets: Dictionary of set_name -> gene list
        te_data: TE data from load_strand_bias_data()
        all_genes: Set of all genome genes

    Returns:
        DataFrame with enrichment results
    """
    # Identify genes with TE hits
    genes_with_te = {fbgn for fbgn, data in te_data.items() if data['total_hits'] > 0}

    # Compute background densities (all genes not in any specific set)
    all_densities = []
    for fbgn in all_genes:
        data = te_data.get(fbgn)
        if data and data['total_hits'] > 0 and data['total_utr_len'] > 0:
            density = (data['total_hits'] / data['total_utr_len']) * 1000
            all_densities.append(density)

    results = []

    for set_name, genes in sorted(gene_sets.items()):
        gene_set = set(genes) & all_genes  # Intersect with genes we have data for

        if len(gene_set) < 5:
            continue

        # Basic metrics
        metrics = compute_te_metrics(gene_set, te_data)

        # Fisher's exact test
        fisher = fisher_exact_te_enrichment(gene_set, all_genes, genes_with_te)

        # Mann-Whitney U test
        background_genes = all_genes - gene_set
        background_densities = []
        for fbgn in background_genes:
            data = te_data.get(fbgn)
            if data and data['total_hits'] > 0 and data['total_utr_len'] > 0:
                density = (data['total_hits'] / data['total_utr_len']) * 1000
                background_densities.append(density)

        mw = mann_whitney_te_density(metrics['densities'], background_densities)

        # Categorize
        if set_name.startswith('expr_'):
            category = 'expression'
        elif set_name.startswith('go_'):
            category = 'go'
        elif set_name.startswith('flyfish_'):
            category = 'flyfish'
        elif set_name.startswith('group_'):
            category = 'gene_group'
        else:
            category = 'other'

        results.append({
            'gene_set': set_name,
            'category': category,
            'n_genes': metrics['n_genes'],
            'n_with_hits': metrics['n_with_hits'],
            'pct_with_hits': metrics['pct_with_hits'],
            'mean_density': metrics['mean_density'],
            'median_density': metrics['median_density'],
            'sense_pct': metrics['sense_pct'],
            'fisher_or': fisher['odds_ratio'],
            'fisher_p': fisher['p_value'],
            'mw_median_in': mw['median_in'],
            'mw_median_out': mw['median_out'],
            'mw_p': mw['p_value'],
        })

    df = pd.DataFrame(results)

    # Apply FDR correction
    if len(df) > 0:
        df['fisher_q'] = apply_fdr_correction(df['fisher_p'].tolist())
        mw_p_valid = df['mw_p'].fillna(1.0).tolist()
        df['mw_q'] = apply_fdr_correction(mw_p_valid)

    return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--gene-sets',
        type=Path,
        default=None,
        help='Directory containing gene sets (default: data/gene_lists/functional)'
    )
    parser.add_argument(
        '--te-data',
        type=Path,
        default=None,
        help='TE data file (default: results/strand_bias_by_utr.tsv)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output file (default: results/functional_te_enrichment.tsv)'
    )
    parser.add_argument(
        '--sig-threshold',
        type=float,
        default=0.05,
        help='Significance threshold for q-values (default: 0.05)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()

    # Set paths
    if args.gene_sets:
        gene_sets_dir = args.gene_sets
    else:
        gene_sets_dir = project_root / "data" / "gene_lists" / "functional"

    if args.te_data:
        te_data_path = args.te_data
    else:
        te_data_path = results_dir / "strand_bias_by_utr.tsv"

    if args.output:
        output_path = args.output
    else:
        output_path = results_dir / "functional_te_enrichment.tsv"

    # Check inputs
    if not gene_sets_dir.exists():
        print(f"Error: Gene sets directory not found: {gene_sets_dir}")
        print("Run build_functional_gene_sets.py first.")
        return 1

    if not te_data_path.exists():
        print(f"Error: TE data file not found: {te_data_path}")
        return 1

    # Load data
    print(f"Loading gene sets from: {gene_sets_dir}")
    gene_sets = load_gene_sets(gene_sets_dir)
    print(f"  Loaded {len(gene_sets)} gene sets")

    print(f"Loading TE data from: {te_data_path}")
    te_data = load_strand_bias_data(te_data_path)
    print(f"  Loaded data for {len(te_data):,} genes")

    # Get all genes
    all_genes = set(te_data.keys())
    for genes in gene_sets.values():
        all_genes.update(genes)
    print(f"  Total genes: {len(all_genes):,}")

    # Compute genome-wide baseline
    print()
    print("Genome-wide baseline:")
    baseline = compute_te_metrics(all_genes, te_data)
    print(f"  Genes with hits: {baseline['n_with_hits']:,} / {baseline['n_genes']:,} ({baseline['pct_with_hits']:.1f}%)")
    print(f"  Mean density: {baseline['mean_density']:.1f} hits/kb")
    print(f"  Median density: {baseline['median_density']:.1f} hits/kb")
    print(f"  Strand bias: {baseline['sense_pct']:.1f}% sense / {baseline['anti_pct']:.1f}% antisense")

    # Run enrichment analysis
    print()
    print("Running enrichment analysis...")
    df = analyze_enrichment(gene_sets, te_data, all_genes)

    # Save results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.4g')
    print(f"Results saved to: {output_path}")
    print(f"  {len(df)} gene sets analyzed")

    # Summary
    print()
    print("Significant results (q < {:.2f}):".format(args.sig_threshold))

    sig_fisher = df[df['fisher_q'] < args.sig_threshold].sort_values('fisher_q')
    sig_mw = df[df['mw_q'] < args.sig_threshold].sort_values('mw_q')

    print(f"\nFisher's exact test (TE presence):")
    if len(sig_fisher) > 0:
        for _, row in sig_fisher.head(10).iterrows():
            direction = "enriched" if row['fisher_or'] > 1 else "depleted"
            print(f"  {row['gene_set']}: OR={row['fisher_or']:.2f} ({direction}), q={row['fisher_q']:.2e}")
    else:
        print("  No significant results")

    print(f"\nMann-Whitney U test (TE density):")
    if len(sig_mw) > 0:
        for _, row in sig_mw.head(10).iterrows():
            direction = "higher" if row['mw_median_in'] > row['mw_median_out'] else "lower"
            print(f"  {row['gene_set']}: median={row['mw_median_in']:.1f} vs {row['mw_median_out']:.1f} ({direction}), q={row['mw_q']:.2e}")
    else:
        print("  No significant results")

    # Category summary
    print()
    print("Results by category:")
    for cat in df['category'].unique():
        cat_df = df[df['category'] == cat]
        sig_count = ((cat_df['fisher_q'] < args.sig_threshold) | (cat_df['mw_q'] < args.sig_threshold)).sum()
        print(f"  {cat}: {len(cat_df)} sets, {sig_count} significant")

    return 0


if __name__ == '__main__':
    sys.exit(main())
