#!/usr/bin/env python3
"""
Analyze isoform-specific TE motifs.

Identifies TE hits with motifs that differ between isoforms of the same gene,
potentially indicating alternative 3'UTR regulation.

Output: results/motif_analysis/isoform_specific_te_motifs.tsv
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, List, Tuple

import pandas as pd

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir, get_references_dir
from utils.blast_io import iter_blast_results
from utils.data_loaders import parse_fasta_by_parent


def find_motifs_in_sequence(sequence: str, motifs: List[str]) -> Set[str]:
    """Find which motifs are present in a sequence."""
    clean_seq = sequence.replace('-', '').upper()
    found = set()
    for motif in motifs:
        if motif.upper() in clean_seq:
            found.add(motif)
    return found


def build_gene_isoform_map(fasta_path: Path) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Build transcript-gene and gene-isoform mappings.

    Returns:
        Tuple of:
        - transcript_to_gene: FBtr -> FBgn
        - gene_to_isoforms: FBgn -> list of FBtr
    """
    transcripts = parse_fasta_by_parent(fasta_path)

    transcript_to_gene = {}
    gene_to_isoforms = defaultdict(list)

    for fbtr, info in transcripts.items():
        fbgn = info.get('parent')
        if fbgn:
            transcript_to_gene[fbtr] = fbgn
            gene_to_isoforms[fbgn].append(fbtr)

    return transcript_to_gene, dict(gene_to_isoforms)


def analyze_isoform_motifs(
    blast_file: Path,
    motifs: List[str],
    transcript_to_gene: Dict[str, str],
    gene_to_isoforms: Dict[str, List[str]]
) -> Dict:
    """
    Analyze which isoforms have TE hits containing each motif.

    Returns:
        Dictionary: gene -> motif -> set of isoforms with that motif
    """
    # gene -> motif -> set of isoforms
    gene_motif_isoforms = defaultdict(lambda: defaultdict(set))

    # Also track all isoforms seen per gene
    gene_isoforms_seen = defaultdict(set)

    for hit in iter_blast_results(blast_file):
        qseqid = hit.get('qseqid', '')
        qseq = hit.get('qseq', '')

        if not qseq:
            continue

        # Handle shuffled IDs
        if '_shuf' in qseqid:
            continue  # Skip shuffled data

        fbtr = qseqid
        fbgn = transcript_to_gene.get(fbtr)
        if not fbgn:
            continue

        # Track that we've seen this isoform
        gene_isoforms_seen[fbgn].add(fbtr)

        # Find motifs in this hit
        found_motifs = find_motifs_in_sequence(qseq, motifs)

        for motif in found_motifs:
            gene_motif_isoforms[fbgn][motif].add(fbtr)

    return gene_motif_isoforms, gene_isoforms_seen


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
        help='Output file (default: results/motif_analysis/isoform_specific_te_motifs.tsv)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    refs_dir = get_references_dir()
    motif_dir = results_dir / "motif_analysis"

    # Set paths
    motif_results_path = args.motif_results or motif_dir / "motif_enrichment_6mer.tsv"
    blast_path = args.blast or results_dir / "genome_wide_all_3utrs.tsv"
    utr_fasta_path = args.utr_fasta or refs_dir / "dmel_3utr.fasta"
    output_path = args.output or motif_dir / "isoform_specific_te_motifs.tsv"

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

    if not utr_fasta_path.exists():
        print(f"Error: UTR FASTA not found: {utr_fasta_path}")
        return 1

    print("=" * 70)
    print("ISOFORM-SPECIFIC TE MOTIF ANALYSIS")
    print("=" * 70)

    # Load motif results
    print(f"\nLoading motif results from: {motif_results_path}")
    motif_df = pd.read_csv(motif_results_path, sep='\t')
    print(f"  Total motifs: {len(motif_df)}")

    # Select top motifs
    sig_motifs = motif_df[motif_df['q_value'] < 0.05]
    top_enriched = sig_motifs.nlargest(args.top_n, 'enrichment')['motif'].tolist()
    top_depleted = sig_motifs.nsmallest(args.top_n, 'enrichment')['motif'].tolist()

    # Add known regulatory motifs
    known_motifs = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'ATTTTT', 'TTTTTA']
    motifs_to_analyze = list(set(top_enriched + top_depleted + known_motifs))
    print(f"  Analyzing {len(motifs_to_analyze)} motifs")

    # Build gene-isoform mappings
    print(f"\nLoading isoform data from: {utr_fasta_path}")
    transcript_to_gene, gene_to_isoforms = build_gene_isoform_map(utr_fasta_path)
    print(f"  {len(transcript_to_gene)} transcripts")
    print(f"  {len(gene_to_isoforms)} genes")

    # Count genes with multiple isoforms
    multi_isoform_genes = {g: isos for g, isos in gene_to_isoforms.items() if len(isos) >= 2}
    print(f"  {len(multi_isoform_genes)} genes with 2+ isoforms")

    # Analyze motifs per isoform
    print(f"\nAnalyzing isoform-specific motifs...")
    gene_motif_isoforms, gene_isoforms_seen = analyze_isoform_motifs(
        blast_path, motifs_to_analyze, transcript_to_gene, gene_to_isoforms
    )

    # Find isoform-specific motifs
    print("\nIdentifying isoform-specific patterns...")
    results = []

    for fbgn, known_isoforms in multi_isoform_genes.items():
        # Get isoforms we actually saw in BLAST results
        seen_isoforms = gene_isoforms_seen.get(fbgn, set())

        if len(seen_isoforms) < 2:
            continue  # Need at least 2 isoforms with TE hits

        motif_data = gene_motif_isoforms.get(fbgn, {})

        for motif in motifs_to_analyze:
            isoforms_with = motif_data.get(motif, set())

            # Only interested if motif is in some but not all isoforms
            if len(isoforms_with) == 0:
                continue

            isoforms_without = seen_isoforms - isoforms_with

            if len(isoforms_without) == 0:
                # Motif in all isoforms - not isoform-specific
                continue

            # Get motif info
            motif_info = motif_df[motif_df['motif'] == motif]
            if len(motif_info) > 0:
                motif_enrichment = motif_info.iloc[0]['enrichment']
                motif_q = motif_info.iloc[0]['q_value']
            else:
                motif_enrichment = 1.0
                motif_q = 1.0

            # Determine specificity type
            if len(isoforms_with) == 1:
                specificity = 'unique_to_one'
            elif len(isoforms_with) < len(seen_isoforms) / 2:
                specificity = 'minority'
            else:
                specificity = 'majority'

            results.append({
                'fbgn': fbgn,
                'n_isoforms_total': len(known_isoforms),
                'n_isoforms_with_te': len(seen_isoforms),
                'motif': motif,
                'motif_enrichment': motif_enrichment,
                'motif_q_value': motif_q,
                'n_isoforms_with_motif': len(isoforms_with),
                'n_isoforms_without_motif': len(isoforms_without),
                'isoforms_with': ','.join(sorted(isoforms_with)),
                'isoforms_without': ','.join(sorted(isoforms_without)),
                'specificity': specificity,
            })

    # Create DataFrame
    df = pd.DataFrame(results)

    if len(df) > 0:
        # Sort by gene, then by number of isoforms with motif (ascending = more specific first)
        df = df.sort_values(['fbgn', 'n_isoforms_with_motif'])

    # Save results
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6g')
    print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if len(df) == 0:
        print("\nNo isoform-specific motifs found")
        return 0

    n_genes = df['fbgn'].nunique()
    n_unique = len(df[df['specificity'] == 'unique_to_one'])

    print(f"\nTotal isoform-specific motif instances: {len(df)}")
    print(f"Genes with isoform-specific motifs: {n_genes}")
    print(f"Unique-to-one-isoform cases: {n_unique}")

    # Summary by motif
    print("\n" + "=" * 70)
    print("MOTIFS BY ISOFORM SPECIFICITY")
    print("=" * 70)

    motif_summary = df.groupby('motif').agg({
        'fbgn': 'nunique',
        'specificity': lambda x: (x == 'unique_to_one').sum(),
        'n_isoforms_with_motif': 'mean',
    }).round(2)
    motif_summary.columns = ['n_genes', 'n_unique', 'mean_isoforms_with']
    motif_summary = motif_summary.sort_values('n_genes', ascending=False)

    print(f"\n{'Motif':<10} {'N genes':>10} {'N unique':>10} {'Mean iso with':>15}")
    print("-" * 50)
    for motif, row in motif_summary.head(20).iterrows():
        print(f"{motif:<10} {row['n_genes']:>10} {row['n_unique']:>10} {row['mean_isoforms_with']:>15.2f}")

    # Top genes with most isoform-specific motifs
    print("\n" + "=" * 70)
    print("GENES WITH MOST ISOFORM-SPECIFIC MOTIFS")
    print("=" * 70)

    gene_summary = df.groupby('fbgn').agg({
        'motif': 'nunique',
        'n_isoforms_total': 'first',
        'n_isoforms_with_te': 'first',
        'specificity': lambda x: (x == 'unique_to_one').sum(),
    })
    gene_summary.columns = ['n_specific_motifs', 'n_isoforms', 'n_iso_with_te', 'n_unique_motifs']
    gene_summary = gene_summary.sort_values('n_specific_motifs', ascending=False)

    print(f"\n{'FBgn':<15} {'N motifs':>10} {'N isoforms':>12} {'N unique':>10}")
    print("-" * 50)
    for fbgn, row in gene_summary.head(20).iterrows():
        print(f"{fbgn:<15} {row['n_specific_motifs']:>10} {row['n_isoforms']:>12} {row['n_unique_motifs']:>10}")

    # Example isoform-specific cases
    print("\n" + "=" * 70)
    print("EXAMPLE ISOFORM-SPECIFIC MOTIFS")
    print("=" * 70)

    unique_cases = df[df['specificity'] == 'unique_to_one'].head(10)
    if len(unique_cases) > 0:
        print("\nMotifs unique to single isoform (first 10):")
        print("-" * 80)
        for _, row in unique_cases.iterrows():
            print(f"Gene: {row['fbgn']}")
            print(f"  Motif: {row['motif']} (enrichment: {row['motif_enrichment']:.2f}x)")
            print(f"  Present in: {row['isoforms_with']}")
            print(f"  Absent from: {row['isoforms_without'][:50]}..." if len(row['isoforms_without']) > 50
                  else f"  Absent from: {row['isoforms_without']}")
            print()

    # Known regulatory motif examples
    print("\n" + "=" * 70)
    print("KNOWN REGULATORY MOTIFS - ISOFORM SPECIFICITY")
    print("=" * 70)

    known_check = ['AATAAA', 'ATTAAA', 'ATTTTT']
    for motif in known_check:
        motif_data = df[df['motif'] == motif]
        if len(motif_data) > 0:
            n_genes = motif_data['fbgn'].nunique()
            n_unique = len(motif_data[motif_data['specificity'] == 'unique_to_one'])
            print(f"\n{motif}:")
            print(f"  Genes with isoform-specific instances: {n_genes}")
            print(f"  Unique-to-one-isoform cases: {n_unique}")
        else:
            print(f"\n{motif}: No isoform-specific instances found")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
