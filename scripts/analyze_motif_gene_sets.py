#!/usr/bin/env python3
"""
Test if specific TE motifs are enriched in functional gene sets.

For each top motif, identifies genes with TE hits containing that motif,
then tests enrichment against each gene set using Fisher's exact test.

Output: results/motif_analysis/motif_gene_set_enrichment.tsv
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, List

import numpy as np
import pandas as pd
from scipy import stats

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir, get_references_dir
from utils.blast_io import iter_blast_results
from utils.data_loaders import load_gene_list, parse_fasta_by_parent


def find_motif_in_sequence(sequence: str, motif: str) -> bool:
    """Check if motif exists in sequence (after removing gaps)."""
    clean_seq = sequence.replace('-', '').upper()
    return motif.upper() in clean_seq


def load_gene_sets(gene_sets_dir: Path) -> Dict[str, Set[str]]:
    """
    Load all gene sets from directory.

    Returns:
        Dictionary mapping set_name -> set of FBgn IDs
    """
    gene_sets = {}

    for path in sorted(gene_sets_dir.glob("*_fbgn_ids.txt")):
        set_name = path.stem.replace('_fbgn_ids', '')
        genes = load_gene_list(path)
        if genes:
            gene_sets[set_name] = set(genes)

    return gene_sets


def build_transcript_to_gene_map(fasta_path: Path) -> Dict[str, str]:
    """
    Build transcript ID to gene ID mapping from FASTA.

    Returns:
        Dictionary mapping FBtr -> FBgn
    """
    transcripts = parse_fasta_by_parent(fasta_path)
    return {fbtr: info['parent'] for fbtr, info in transcripts.items() if info['parent']}


def find_genes_with_motif(
    blast_file: Path,
    motifs: List[str],
    transcript_to_gene: Dict[str, str]
) -> Dict[str, Set[str]]:
    """
    Find genes with TE hits containing each motif.

    Args:
        blast_file: BLAST results file
        motifs: List of motifs to search for
        transcript_to_gene: Mapping from transcript to gene IDs

    Returns:
        Dictionary mapping motif -> set of FBgn IDs
    """
    motif_genes = {m: set() for m in motifs}

    for hit in iter_blast_results(blast_file):
        qseqid = hit.get('qseqid', '')
        qseq = hit.get('qseq', '')

        if not qseq:
            continue

        # Extract transcript ID (handle shuffled IDs like FBtr0070000_shuf1)
        if '_shuf' in qseqid:
            fbtr = qseqid.split('_shuf')[0]
        else:
            fbtr = qseqid

        # Get gene ID
        fbgn = transcript_to_gene.get(fbtr)
        if not fbgn:
            continue

        # Check each motif
        for motif in motifs:
            if find_motif_in_sequence(qseq, motif):
                motif_genes[motif].add(fbgn)

    return motif_genes


def fisher_exact_test(
    motif_genes: Set[str],
    gene_set: Set[str],
    all_genes: Set[str]
) -> Dict:
    """
    Fisher's exact test for motif enrichment in gene set.

    Contingency table:
                        | In gene set | Not in set |
    Has motif TE hit   |     a       |     b      |
    No motif TE hit    |     c       |     d      |

    Returns:
        Dictionary with odds ratio, p-value, and overlap count
    """
    a = len(motif_genes & gene_set)
    b = len(motif_genes - gene_set)
    c = len(gene_set - motif_genes)
    d = len((all_genes - motif_genes) - gene_set)

    # Avoid divide by zero
    if (a + b) == 0 or (c + d) == 0:
        return {
            'overlap': a,
            'odds_ratio': np.nan,
            'p_value': 1.0,
        }

    table = [[a, b], [c, d]]
    odds_ratio, p_value = stats.fisher_exact(table, alternative='two-sided')

    return {
        'overlap': a,
        'odds_ratio': odds_ratio,
        'p_value': p_value,
    }


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--motif-results',
        type=Path,
        default=None,
        help='Motif enrichment results (default: results/motif_analysis/motif_enrichment_6mer.tsv)'
    )
    parser.add_argument(
        '--blast',
        type=Path,
        default=None,
        help='BLAST results file (default: results/genome_wide_all_3utrs.tsv)'
    )
    parser.add_argument(
        '--gene-sets',
        type=Path,
        default=None,
        help='Gene sets directory (default: data/gene_lists/functional)'
    )
    parser.add_argument(
        '--utr-fasta',
        type=Path,
        default=None,
        help='3\'UTR FASTA file (default: data/references/dmel_3utr.fasta)'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=30,
        help='Number of top motifs to analyze (default: 30)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output file (default: results/motif_analysis/motif_gene_set_enrichment.tsv)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    refs_dir = get_references_dir()
    motif_dir = results_dir / "motif_analysis"

    # Set paths
    motif_results_path = args.motif_results or motif_dir / "motif_enrichment_6mer.tsv"
    blast_path = args.blast or results_dir / "genome_wide_all_3utrs.tsv"
    gene_sets_dir = args.gene_sets or project_root / "data" / "gene_lists" / "functional"
    utr_fasta_path = args.utr_fasta or refs_dir / "dmel_3utr.fasta"
    output_path = args.output or motif_dir / "motif_gene_set_enrichment.tsv"

    # Create output directory
    motif_dir.mkdir(parents=True, exist_ok=True)

    # Check inputs
    if not motif_results_path.exists():
        print(f"Error: Motif results not found: {motif_results_path}")
        print("Run analyze_te_motifs.py first.")
        return 1

    if not blast_path.exists():
        print(f"Error: BLAST file not found: {blast_path}")
        return 1

    if not gene_sets_dir.exists():
        print(f"Error: Gene sets directory not found: {gene_sets_dir}")
        return 1

    if not utr_fasta_path.exists():
        print(f"Error: UTR FASTA not found: {utr_fasta_path}")
        return 1

    print("=" * 70)
    print("MOTIF-GENE SET ENRICHMENT ANALYSIS")
    print("=" * 70)

    # Load motif results
    print(f"\nLoading motif results from: {motif_results_path}")
    motif_df = pd.read_csv(motif_results_path, sep='\t')
    print(f"  Total motifs: {len(motif_df)}")

    # Select top motifs (by significance)
    sig_motifs = motif_df[motif_df['q_value'] < 0.05].nsmallest(args.top_n * 2, 'p_value')

    # Take top enriched and depleted
    top_enriched = sig_motifs.nlargest(args.top_n, 'enrichment')['motif'].tolist()
    top_depleted = sig_motifs.nsmallest(args.top_n, 'enrichment')['motif'].tolist()
    motifs_to_analyze = list(set(top_enriched + top_depleted))

    print(f"  Analyzing {len(motifs_to_analyze)} motifs")

    # Load transcript-to-gene mapping
    print(f"\nLoading transcript-gene mapping from: {utr_fasta_path}")
    transcript_to_gene = build_transcript_to_gene_map(utr_fasta_path)
    print(f"  {len(transcript_to_gene)} transcripts mapped")

    # Get all genes with any UTR data
    all_genes = set(transcript_to_gene.values())
    print(f"  {len(all_genes)} unique genes")

    # Load gene sets
    print(f"\nLoading gene sets from: {gene_sets_dir}")
    gene_sets = load_gene_sets(gene_sets_dir)
    print(f"  {len(gene_sets)} gene sets loaded")

    # Find genes with each motif
    print(f"\nFinding genes with motif-containing TE hits...")
    motif_genes = find_genes_with_motif(blast_path, motifs_to_analyze, transcript_to_gene)

    for motif in motifs_to_analyze[:5]:
        print(f"  {motif}: {len(motif_genes[motif])} genes")
    if len(motifs_to_analyze) > 5:
        print(f"  ... and {len(motifs_to_analyze) - 5} more motifs")

    # Run enrichment tests
    print(f"\nRunning Fisher's exact tests ({len(motifs_to_analyze)} motifs x {len(gene_sets)} gene sets)...")
    results = []
    n_tests = 0

    for motif in motifs_to_analyze:
        motif_info = motif_df[motif_df['motif'] == motif].iloc[0]
        genes_with_motif = motif_genes[motif]

        if len(genes_with_motif) < 5:
            continue

        for set_name, gene_set in gene_sets.items():
            # Intersect gene set with genes we have data for
            gene_set_filtered = gene_set & all_genes

            if len(gene_set_filtered) < 10:
                continue

            # Run test
            test_result = fisher_exact_test(genes_with_motif, gene_set_filtered, all_genes)
            n_tests += 1

            results.append({
                'motif': motif,
                'motif_enrichment': motif_info['enrichment'],
                'motif_q_value': motif_info['q_value'],
                'gene_set': set_name,
                'n_genes_set': len(gene_set_filtered),
                'n_genes_motif': len(genes_with_motif),
                'overlap_count': test_result['overlap'],
                'fisher_or': test_result['odds_ratio'],
                'fisher_p': test_result['p_value'],
            })

    print(f"  {n_tests} tests performed")

    # Create DataFrame and apply FDR correction
    df = pd.DataFrame(results)

    if len(df) > 0:
        try:
            from statsmodels.stats.multitest import multipletests
            _, q_values, _, _ = multipletests(df['fisher_p'].values, method='fdr_bh')
            df['fisher_q'] = q_values
        except ImportError:
            df['fisher_q'] = np.minimum(df['fisher_p'] * len(df), 1.0)

    # Sort by significance
    df = df.sort_values('fisher_p')

    # Save results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    sig_results = df[df['fisher_q'] < 0.05]
    sig_enriched = sig_results[sig_results['fisher_or'] > 1]
    sig_depleted = sig_results[sig_results['fisher_or'] < 1]

    print(f"\nSignificant associations (q < 0.05):")
    print(f"  Total: {len(sig_results)}")
    print(f"  Enriched (OR > 1): {len(sig_enriched)}")
    print(f"  Depleted (OR < 1): {len(sig_depleted)}")

    # Top enriched associations
    print("\nTop 20 motif-gene set associations (enriched):")
    print("-" * 90)
    top_assoc = sig_enriched.nlargest(20, 'fisher_or')
    if len(top_assoc) > 0:
        print(f"{'Motif':<10} {'Gene Set':<35} {'Overlap':>8} {'OR':>8} {'Q-value':>12}")
        print("-" * 90)
        for _, row in top_assoc.iterrows():
            gene_set_short = row['gene_set'][:33] + '..' if len(row['gene_set']) > 35 else row['gene_set']
            print(f"{row['motif']:<10} {gene_set_short:<35} {row['overlap_count']:>8} "
                  f"{row['fisher_or']:>8.2f} {row['fisher_q']:>12.2e}")
    else:
        print("  No significant enriched associations found")

    # Summarize by motif
    print("\n" + "=" * 70)
    print("MOTIFS WITH MOST GENE SET ASSOCIATIONS")
    print("=" * 70)

    motif_summary = sig_results.groupby('motif').agg({
        'gene_set': 'count',
        'fisher_or': ['mean', 'max'],
        'fisher_q': 'min',
    }).round(3)
    motif_summary.columns = ['n_associations', 'mean_or', 'max_or', 'min_q']
    motif_summary = motif_summary.sort_values('n_associations', ascending=False)

    if len(motif_summary) > 0:
        print(f"\n{'Motif':<10} {'N assoc':>10} {'Mean OR':>10} {'Max OR':>10} {'Min Q':>12}")
        print("-" * 55)
        for motif, row in motif_summary.head(15).iterrows():
            print(f"{motif:<10} {row['n_associations']:>10} {row['mean_or']:>10.2f} "
                  f"{row['max_or']:>10.2f} {row['min_q']:>12.2e}")

    # Summarize by gene set
    print("\n" + "=" * 70)
    print("GENE SETS WITH MOST MOTIF ASSOCIATIONS")
    print("=" * 70)

    geneset_summary = sig_results.groupby('gene_set').agg({
        'motif': 'count',
        'fisher_or': ['mean', 'max'],
        'fisher_q': 'min',
    }).round(3)
    geneset_summary.columns = ['n_motifs', 'mean_or', 'max_or', 'min_q']
    geneset_summary = geneset_summary.sort_values('n_motifs', ascending=False)

    if len(geneset_summary) > 0:
        print(f"\n{'Gene Set':<40} {'N motifs':>10} {'Mean OR':>10} {'Max OR':>10}")
        print("-" * 75)
        for gene_set, row in geneset_summary.head(15).iterrows():
            gene_set_short = gene_set[:38] + '..' if len(gene_set) > 40 else gene_set
            print(f"{gene_set_short:<40} {row['n_motifs']:>10} {row['mean_or']:>10.2f} "
                  f"{row['max_or']:>10.2f}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
