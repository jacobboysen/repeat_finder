#!/usr/bin/env python3
"""
Compare synteny of TE hits that contain specific motifs vs those that don't.

Uses existing pre-computed synteny data from MAF alignments.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--motifs', type=str,
                        default='AATAAA,ATTAAA,TGTATA,TGTAAA,TATTTA,CAGCAG,AGCAGC,ACACAC,CACACA,GGTAAG',
                        help='Comma-separated list of motifs to analyze')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    motif_dir = results_dir / "motif_analysis"
    motif_dir.mkdir(exist_ok=True)

    output_path = args.output or motif_dir / "motif_synteny_by_te_hit.tsv"

    # File paths
    blast_path = results_dir / "genome_wide_all_3utrs.tsv"
    synteny_hq_path = results_dir / "repeatmasker_analysis" / "te_hits_hq_synteny.tsv"
    synteny_all_path = results_dir / "repeatmasker_analysis" / "te_hits_all_synteny_sampled.tsv"

    motifs = [m.strip().upper() for m in args.motifs.split(',')]

    print("=" * 70)
    print("MOTIF SYNTENY ANALYSIS")
    print("=" * 70)
    print(f"Motifs: {', '.join(motifs)}")

    # Load synteny data
    print("\n1. Loading synteny data...")

    # HQ hits synteny
    synteny_hq = pd.read_csv(synteny_hq_path, sep='\t')
    print(f"   HQ synteny data: {len(synteny_hq):,} hits")

    # All hits synteny (sampled)
    synteny_all = pd.read_csv(synteny_all_path, sep='\t')
    print(f"   All hits synteny (sampled): {len(synteny_all):,} hits")

    # Parse synteny names to extract transcript/TE info
    # Format: novel_123|FBtr0123|FBti0456|80.5|50

    def parse_synteny_name(name):
        parts = name.split('|')
        if len(parts) >= 5:
            return {
                'category': parts[0].split('_')[0],
                'transcript': parts[1],
                'te': parts[2],
                'pident': float(parts[3]),
                'length': int(parts[4])
            }
        return None

    # Build lookup for synteny
    print("\n2. Building synteny lookup...")

    synteny_lookup = {}
    for _, row in synteny_all.iterrows():
        parsed = parse_synteny_name(row['name'])
        if parsed:
            key = (parsed['transcript'], parsed['te'],
                   round(parsed['pident'], 1), parsed['length'])
            synteny_lookup[key] = {
                'any_syntenic': row['any_coverage'] >= 0.5,
                'sim_coverage': row['sim_coverage'],
                'yak_coverage': row['yak_coverage'],
                'any_coverage': row['any_coverage'],
            }

    print(f"   Built lookup with {len(synteny_lookup):,} entries")

    # Load BLAST results and check motif presence
    print("\n3. Loading BLAST results and matching synteny...")

    hit_data = []
    n_processed = 0
    n_matched = 0

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        sseqid = hit.get('sseqid', '')
        pident = hit.get('pident', 0)
        length = hit.get('length', 0)
        qseq = hit.get('qseq', '')

        if not qseq:
            continue

        n_processed += 1
        if n_processed % 500000 == 0:
            print(f"   Processed {n_processed:,} hits...")

        key = (qseqid, sseqid, round(pident, 1), length)

        if key not in synteny_lookup:
            continue

        n_matched += 1
        syn = synteny_lookup[key]
        clean_seq = qseq.replace('-', '').upper()

        record = {
            'transcript': qseqid,
            'te': sseqid,
            'pident': pident,
            'length': length,
            'any_syntenic': syn['any_syntenic'],
            'any_coverage': syn['any_coverage'],
            'sim_coverage': syn['sim_coverage'],
        }

        # Check each motif
        for motif in motifs:
            record[f'has_{motif}'] = motif in clean_seq

        hit_data.append(record)

    print(f"   Processed: {n_processed:,}, Matched synteny: {n_matched:,}")

    # Convert to DataFrame
    df = pd.DataFrame(hit_data)
    print(f"\n   Hits with synteny data: {len(df):,}")

    # Analyze synteny by motif presence
    print("\n" + "=" * 70)
    print("RESULTS: Synteny of TE Hits by Motif Presence")
    print("=" * 70)

    overall_syntenic = df['any_syntenic'].mean() * 100
    overall_coverage = df['any_coverage'].mean()
    print(f"\nOverall: {overall_syntenic:.1f}% syntenic, mean coverage: {overall_coverage:.3f}")

    results = []

    for motif in motifs:
        col = f'has_{motif}'
        if col not in df.columns:
            continue

        with_motif = df[df[col] == True]
        without_motif = df[df[col] == False]

        if len(with_motif) < 10:
            continue

        with_syntenic = with_motif['any_syntenic'].mean() * 100
        without_syntenic = without_motif['any_syntenic'].mean() * 100

        with_cov = with_motif['any_coverage'].mean()
        without_cov = without_motif['any_coverage'].mean()

        results.append({
            'motif': motif,
            'n_with': len(with_motif),
            'n_without': len(without_motif),
            'pct_with': len(with_motif) / len(df) * 100,
            'pct_syntenic_with': with_syntenic,
            'pct_syntenic_without': without_syntenic,
            'syntenic_diff': with_syntenic - without_syntenic,
            'mean_cov_with': with_cov,
            'mean_cov_without': without_cov,
            'cov_diff': with_cov - without_cov,
        })

    results_df = pd.DataFrame(results)

    # Print results
    print(f"\n{'Motif':<8} {'N with':>10} {'N w/o':>10} {'%Syn(w)':>9} {'%Syn(w/o)':>10} {'SynDiff':>9} {'Cov(w)':>8} {'Cov(w/o)':>9} {'CovDiff':>9}")
    print("-" * 105)

    for _, row in results_df.iterrows():
        syn_diff = f"{row['syntenic_diff']:+.2f}%"
        cov_diff = f"{row['cov_diff']:+.4f}"
        print(f"{row['motif']:<8} {row['n_with']:>10,} {row['n_without']:>10,} "
              f"{row['pct_syntenic_with']:>8.1f}% {row['pct_syntenic_without']:>9.1f}% "
              f"{syn_diff:>9} {row['mean_cov_with']:>8.4f} {row['mean_cov_without']:>9.4f} {cov_diff:>9}")

    # Save results
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\nResults saved to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    more_syntenic = results_df[results_df['syntenic_diff'] > 0]
    less_syntenic = results_df[results_df['syntenic_diff'] < 0]

    print(f"\nMotifs where TE hits WITH motif are MORE syntenic: {len(more_syntenic)}")
    for _, row in more_syntenic.iterrows():
        print(f"  {row['motif']}: +{row['syntenic_diff']:.2f}% syntenic")

    print(f"\nMotifs where TE hits WITH motif are LESS syntenic: {len(less_syntenic)}")
    for _, row in less_syntenic.iterrows():
        print(f"  {row['motif']}: {row['syntenic_diff']:.2f}% syntenic")

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print("""
Synteny = presence in other Drosophila species (D. simulans, D. yakuba, etc.)
  - Higher synteny = older (pre-speciation) insertion
  - Lower synteny = younger (species-specific) insertion

If regulatory motifs in TE hits are ancient and functional:
  - They should be MORE syntenic (conserved across species)
  - The pattern should match conservation (phyloP) results
""")

    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
