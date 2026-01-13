#!/usr/bin/env python3
"""
Analyze TE family enrichment in germ plasm vs housekeeping genes.

Compares:
1. TE families with hits vs families without hits
2. Germ plasm vs housekeeping enrichment
3. Strand bias per TE family
4. Sense vs antisense query comparison

Outputs comprehensive statistics and visualizations.
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


# BLAST output columns
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'slen', 'qseq', 'sseq'
]


def parse_te_name(description):
    """Extract TE family name from FlyBase TE description."""
    # Look for name= pattern
    match = re.search(r'name=([^;]+)', description)
    if match:
        name = match.group(1)
        # Clean up the name - remove instance identifiers like {}555
        name = re.sub(r'\{[^}]*\}\d*', '', name)
        return name.strip()
    return None


def parse_te_class(name):
    """Classify TE by type based on common naming patterns."""
    name_lower = name.lower()

    # LTR retrotransposons
    ltr_families = ['gypsy', 'copia', 'bel', 'pao', 'mdg', 'roo', '412', '297',
                    'blood', 'accord', 'tirant', 'springer', 'opus', 'diver',
                    'quasimodo', 'idefix', 'invader', 'gtwin', 'tabor', 'stalker']
    for fam in ltr_families:
        if fam in name_lower:
            return 'LTR'

    # Non-LTR retrotransposons (LINEs)
    line_families = ['jockey', 'doc', 'i-element', 'f-element', 'g-element',
                     'x-element', 'het-a', 'tart', 'tahre', 'r1', 'r2', 'cr1',
                     'rt1', 'baggins', 'juan', 'ivk']
    for fam in line_families:
        if fam in name_lower:
            return 'LINE'

    # DNA transposons
    dna_families = ['p-element', 'hobo', 'pogo', 'bari', 's-element',
                    'transib', 'tc1', 'mariner', 'piggybac', 'helitron',
                    'mite', 'protop', 'dnarep']
    for fam in dna_families:
        if fam in name_lower:
            return 'DNA'

    # SINEs
    if 'sine' in name_lower or name_lower.startswith('ine-'):
        return 'SINE'

    # Terminal repeats
    if 'ltr' in name_lower:
        return 'LTR'

    # Default
    return 'Unknown'


def load_te_database(te_fasta):
    """Load all TE families from the database."""
    if SeqIO is None:
        print("Warning: BioPython not available, cannot load TE database")
        return {}

    te_info = {}
    for record in SeqIO.parse(te_fasta, 'fasta'):
        te_id = record.id
        name = parse_te_name(record.description)
        if name:
            te_class = parse_te_class(name)
            te_info[te_id] = {
                'name': name,
                'class': te_class,
                'length': len(record.seq)
            }
        else:
            te_info[te_id] = {
                'name': te_id,
                'class': 'Unknown',
                'length': len(record.seq)
            }

    return te_info


def classify_strand(sstart, send):
    """Classify hit strand based on subject coordinates."""
    return 'plus' if sstart < send else 'minus'


def load_blast_results(results_file):
    """Load BLAST results and add strand classification."""
    if not results_file.exists():
        return pd.DataFrame(columns=BLAST_COLUMNS + ['strand'])

    if results_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS + ['strand'])

    # Detect column count
    first_line = open(results_file).readline()
    num_cols = len(first_line.strip().split('\t'))

    if num_cols == 17:
        df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS + ['strand'])
    else:
        df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS)
        df['strand'] = df.apply(lambda row: classify_strand(row['sstart'], row['send']), axis=1)

    return df


def analyze_te_families(blast_df, te_info, gene_set_name='unknown'):
    """Analyze TE family statistics from BLAST results."""
    results = {
        'gene_set': gene_set_name,
        'total_hits': len(blast_df),
        'families_hit': defaultdict(lambda: {'count': 0, 'plus': 0, 'minus': 0,
                                              'mean_evalue': [], 'mean_pident': []})
    }

    # Count hits per TE family
    for _, row in blast_df.iterrows():
        te_id = row['sseqid']
        strand = row.get('strand', 'plus')

        # Get family name from te_info or use ID
        if te_id in te_info:
            family = te_info[te_id]['name']
            te_class = te_info[te_id]['class']
        else:
            family = te_id
            te_class = 'Unknown'

        results['families_hit'][family]['count'] += 1
        results['families_hit'][family][strand] += 1
        results['families_hit'][family]['mean_evalue'].append(row['evalue'])
        results['families_hit'][family]['mean_pident'].append(row['pident'])
        results['families_hit'][family]['class'] = te_class

    # Compute aggregated statistics
    for family, data in results['families_hit'].items():
        data['mean_evalue'] = np.mean(data['mean_evalue']) if data['mean_evalue'] else 0
        data['mean_pident'] = np.mean(data['mean_pident']) if data['mean_pident'] else 0
        data['strand_ratio'] = data['plus'] / data['minus'] if data['minus'] > 0 else float('inf')

    return results


def compare_gene_sets(germ_plasm_results, housekeeping_results, te_info):
    """Compare TE family enrichment between gene sets."""
    comparison = {
        'germ_plasm_only': [],
        'housekeeping_only': [],
        'both': [],
        'enrichment': []
    }

    gp_families = set(germ_plasm_results['families_hit'].keys())
    hk_families = set(housekeeping_results['families_hit'].keys())

    # Families unique to each set
    comparison['germ_plasm_only'] = list(gp_families - hk_families)
    comparison['housekeeping_only'] = list(hk_families - gp_families)
    comparison['both'] = list(gp_families & hk_families)

    # Calculate enrichment for families in both
    for family in comparison['both']:
        gp_count = germ_plasm_results['families_hit'][family]['count']
        hk_count = housekeeping_results['families_hit'][family]['count']

        # Normalize by total hits
        gp_freq = gp_count / germ_plasm_results['total_hits'] if germ_plasm_results['total_hits'] > 0 else 0
        hk_freq = hk_count / housekeeping_results['total_hits'] if housekeeping_results['total_hits'] > 0 else 0

        enrichment = gp_freq / hk_freq if hk_freq > 0 else float('inf')

        comparison['enrichment'].append({
            'family': family,
            'class': germ_plasm_results['families_hit'][family].get('class', 'Unknown'),
            'germ_plasm_count': gp_count,
            'housekeeping_count': hk_count,
            'germ_plasm_freq': gp_freq,
            'housekeeping_freq': hk_freq,
            'enrichment': enrichment,
            'log2_enrichment': np.log2(enrichment) if enrichment > 0 and enrichment != float('inf') else 0
        })

    return comparison


def compare_strands(sense_results, antisense_results):
    """Compare sense vs antisense query results."""
    comparison = {
        'sense_total': sense_results['total_hits'],
        'antisense_total': antisense_results['total_hits'],
        'sense_families': len(sense_results['families_hit']),
        'antisense_families': len(antisense_results['families_hit']),
        'family_comparison': []
    }

    all_families = set(sense_results['families_hit'].keys()) | set(antisense_results['families_hit'].keys())

    for family in all_families:
        sense_data = sense_results['families_hit'].get(family, {'count': 0, 'plus': 0, 'minus': 0})
        anti_data = antisense_results['families_hit'].get(family, {'count': 0, 'plus': 0, 'minus': 0})

        comparison['family_comparison'].append({
            'family': family,
            'sense_count': sense_data['count'],
            'antisense_count': anti_data['count'],
            'sense_plus': sense_data['plus'],
            'sense_minus': sense_data['minus'],
            'antisense_plus': anti_data['plus'],
            'antisense_minus': anti_data['minus']
        })

    return comparison


def generate_summary_report(analysis_results, output_file):
    """Generate summary report as TSV."""
    rows = []

    for family, data in sorted(analysis_results['families_hit'].items(),
                                key=lambda x: -x[1]['count']):
        rows.append({
            'family': family,
            'class': data.get('class', 'Unknown'),
            'total_hits': data['count'],
            'plus_strand': data['plus'],
            'minus_strand': data['minus'],
            'strand_ratio': data['strand_ratio'] if data['strand_ratio'] != float('inf') else 'inf',
            'mean_evalue': data['mean_evalue'],
            'mean_pident': data['mean_pident']
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)
    return df


def analyze_family_strand_preference(all_results, te_info, output_dir):
    """Analyze per-family strand preference across all datasets.

    For each TE family, calculate strand ratio in each dataset to see if
    certain families consistently show strand preference.
    """
    print("\n" + "=" * 60)
    print("Per-Family Strand Preference Analysis")
    print("=" * 60)

    # Collect all families across all datasets
    all_families = set()
    for dataset_name, results in all_results.items():
        if 'families_hit' in results:
            all_families.update(results['families_hit'].keys())

    # Build comprehensive table
    rows = []
    for family in sorted(all_families):
        row = {'family': family}

        # Get class from any dataset that has it
        te_class = 'Unknown'
        for results in all_results.values():
            if family in results.get('families_hit', {}):
                te_class = results['families_hit'][family].get('class', 'Unknown')
                break
        row['class'] = te_class

        # For each dataset, get strand counts
        for dataset_name, results in all_results.items():
            if family in results.get('families_hit', {}):
                fam_data = results['families_hit'][family]
                row[f'{dataset_name}_plus'] = fam_data['plus']
                row[f'{dataset_name}_minus'] = fam_data['minus']
                row[f'{dataset_name}_total'] = fam_data['count']
                ratio = fam_data['plus'] / fam_data['minus'] if fam_data['minus'] > 0 else float('inf')
                row[f'{dataset_name}_ratio'] = ratio if ratio != float('inf') else 999
            else:
                row[f'{dataset_name}_plus'] = 0
                row[f'{dataset_name}_minus'] = 0
                row[f'{dataset_name}_total'] = 0
                row[f'{dataset_name}_ratio'] = None

        rows.append(row)

    df = pd.DataFrame(rows)

    # Save full table
    output_file = output_dir / 'family_strand_preference_all.tsv'
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved: {output_file.name}")

    # Find families with consistent strand bias across germ plasm
    print("\n" + "-" * 50)
    print("Families with strong strand bias (ratio > 2 or < 0.5):")
    print("-" * 50)

    # Focus on germ_plasm_sense if available
    if 'germ_plasm_sense_ratio' in df.columns:
        biased = df[
            (df['germ_plasm_sense_total'] >= 5) &  # At least 5 hits
            ((df['germ_plasm_sense_ratio'] > 2) | (df['germ_plasm_sense_ratio'] < 0.5))
        ].sort_values('germ_plasm_sense_ratio', ascending=False)

        if len(biased) > 0:
            print(f"\n{'Family':<25} {'Class':<10} {'Total':<8} {'Plus':<8} {'Minus':<8} {'Ratio':<8}")
            print("-" * 70)
            for _, row in biased.head(20).iterrows():
                ratio_str = f"{row['germ_plasm_sense_ratio']:.2f}" if row['germ_plasm_sense_ratio'] < 100 else "999+"
                print(f"{row['family']:<25} {row['class']:<10} {row['germ_plasm_sense_total']:<8} "
                      f"{row['germ_plasm_sense_plus']:<8} {row['germ_plasm_sense_minus']:<8} {ratio_str:<8}")

    # Compare strand bias between control groups
    print("\n" + "-" * 50)
    print("Cross-dataset strand ratio comparison (top families by total hits):")
    print("-" * 50)

    # Get datasets with ratio columns
    ratio_cols = [c for c in df.columns if c.endswith('_ratio')]
    total_cols = [c.replace('_ratio', '_total') for c in ratio_cols]

    if len(ratio_cols) >= 2:
        # Sort by total hits across all datasets
        df['total_all'] = df[[c for c in df.columns if c.endswith('_total')]].sum(axis=1)
        top_families = df.nlargest(15, 'total_all')

        # Build comparison table
        header = f"{'Family':<20} {'Class':<8}"
        for col in ratio_cols:
            dataset = col.replace('_ratio', '').replace('_sense', '').replace('_antisense', '')[:10]
            header += f" {dataset:<10}"
        print(header)
        print("-" * (30 + 11 * len(ratio_cols)))

        for _, row in top_families.iterrows():
            line = f"{row['family'][:19]:<20} {row['class'][:7]:<8}"
            for col in ratio_cols:
                val = row[col]
                if val is None or pd.isna(val):
                    line += f" {'N/A':<10}"
                elif val >= 100:
                    line += f" {'999+':<10}"
                else:
                    line += f" {val:.2f}{'':>6}"
            print(line)

    return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--te-database',
        type=Path,
        default=Path('data/references/dmel_te_flybase.fasta'),
        help='TE database FASTA'
    )
    parser.add_argument(
        '--germ-plasm-sense',
        type=Path,
        help='Germ plasm sense BLAST results'
    )
    parser.add_argument(
        '--germ-plasm-antisense',
        type=Path,
        help='Germ plasm antisense BLAST results'
    )
    parser.add_argument(
        '--housekeeping-sense',
        type=Path,
        help='Housekeeping sense BLAST results'
    )
    parser.add_argument(
        '--housekeeping-antisense',
        type=Path,
        help='Housekeeping antisense BLAST results'
    )
    parser.add_argument(
        '--somatic-sense',
        type=Path,
        help='Somatic sense BLAST results'
    )
    parser.add_argument(
        '--cleared-sense',
        type=Path,
        help='Maternally-cleared sense BLAST results'
    )
    parser.add_argument(
        '--adult-sense',
        type=Path,
        help='Adult-specific sense BLAST results'
    )
    parser.add_argument(
        '--shuffled',
        type=Path,
        help='Shuffled control BLAST results'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('results/te_family_analysis'),
        help='Output directory'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    print("TE Family Analysis")
    print("=" * 60)

    # Load TE database
    print("\nLoading TE database...")
    te_info = load_te_database(args.te_database)
    print(f"  Loaded {len(te_info)} TEs from database")

    # Count TE classes
    class_counts = defaultdict(int)
    for te_id, info in te_info.items():
        class_counts[info['class']] += 1

    print("  TE class breakdown:")
    for te_class, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        print(f"    {te_class}: {count}")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    # Analyze germ plasm sense
    if args.germ_plasm_sense and args.germ_plasm_sense.exists():
        print("\nAnalyzing germ plasm sense...")
        gp_sense_df = load_blast_results(args.germ_plasm_sense)
        results['germ_plasm_sense'] = analyze_te_families(gp_sense_df, te_info, 'germ_plasm_sense')
        print(f"  {results['germ_plasm_sense']['total_hits']} hits, {len(results['germ_plasm_sense']['families_hit'])} families")

        # Save summary
        summary_file = args.output_dir / 'germ_plasm_sense_families.tsv'
        generate_summary_report(results['germ_plasm_sense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze germ plasm antisense
    if args.germ_plasm_antisense and args.germ_plasm_antisense.exists():
        print("\nAnalyzing germ plasm antisense...")
        gp_anti_df = load_blast_results(args.germ_plasm_antisense)
        results['germ_plasm_antisense'] = analyze_te_families(gp_anti_df, te_info, 'germ_plasm_antisense')
        print(f"  {results['germ_plasm_antisense']['total_hits']} hits, {len(results['germ_plasm_antisense']['families_hit'])} families")

        summary_file = args.output_dir / 'germ_plasm_antisense_families.tsv'
        generate_summary_report(results['germ_plasm_antisense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze housekeeping sense
    if args.housekeeping_sense and args.housekeeping_sense.exists():
        print("\nAnalyzing housekeeping sense...")
        hk_sense_df = load_blast_results(args.housekeeping_sense)
        results['housekeeping_sense'] = analyze_te_families(hk_sense_df, te_info, 'housekeeping_sense')
        print(f"  {results['housekeeping_sense']['total_hits']} hits, {len(results['housekeeping_sense']['families_hit'])} families")

        summary_file = args.output_dir / 'housekeeping_sense_families.tsv'
        generate_summary_report(results['housekeeping_sense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze housekeeping antisense
    if args.housekeeping_antisense and args.housekeeping_antisense.exists():
        print("\nAnalyzing housekeeping antisense...")
        hk_anti_df = load_blast_results(args.housekeeping_antisense)
        results['housekeeping_antisense'] = analyze_te_families(hk_anti_df, te_info, 'housekeeping_antisense')
        print(f"  {results['housekeeping_antisense']['total_hits']} hits, {len(results['housekeeping_antisense']['families_hit'])} families")

        summary_file = args.output_dir / 'housekeeping_antisense_families.tsv'
        generate_summary_report(results['housekeeping_antisense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze somatic sense
    if args.somatic_sense and args.somatic_sense.exists():
        print("\nAnalyzing somatic sense...")
        somatic_df = load_blast_results(args.somatic_sense)
        results['somatic_sense'] = analyze_te_families(somatic_df, te_info, 'somatic_sense')
        print(f"  {results['somatic_sense']['total_hits']} hits, {len(results['somatic_sense']['families_hit'])} families")

        summary_file = args.output_dir / 'somatic_sense_families.tsv'
        generate_summary_report(results['somatic_sense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze cleared sense
    if args.cleared_sense and args.cleared_sense.exists():
        print("\nAnalyzing maternally-cleared sense...")
        cleared_df = load_blast_results(args.cleared_sense)
        results['cleared_sense'] = analyze_te_families(cleared_df, te_info, 'cleared_sense')
        print(f"  {results['cleared_sense']['total_hits']} hits, {len(results['cleared_sense']['families_hit'])} families")

        summary_file = args.output_dir / 'cleared_sense_families.tsv'
        generate_summary_report(results['cleared_sense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze adult sense
    if args.adult_sense and args.adult_sense.exists():
        print("\nAnalyzing adult-specific sense...")
        adult_df = load_blast_results(args.adult_sense)
        results['adult_sense'] = analyze_te_families(adult_df, te_info, 'adult_sense')
        print(f"  {results['adult_sense']['total_hits']} hits, {len(results['adult_sense']['families_hit'])} families")

        summary_file = args.output_dir / 'adult_sense_families.tsv'
        generate_summary_report(results['adult_sense'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Analyze shuffled control
    if args.shuffled and args.shuffled.exists():
        print("\nAnalyzing shuffled control...")
        shuffled_df = load_blast_results(args.shuffled)
        results['shuffled'] = analyze_te_families(shuffled_df, te_info, 'shuffled')
        print(f"  {results['shuffled']['total_hits']} hits, {len(results['shuffled']['families_hit'])} families")

        summary_file = args.output_dir / 'shuffled_families.tsv'
        generate_summary_report(results['shuffled'], summary_file)
        print(f"  Saved: {summary_file.name}")

    # Run per-family strand preference analysis across all datasets
    if len(results) > 0:
        analyze_family_strand_preference(results, te_info, args.output_dir)

    # Compare gene sets (germ plasm vs housekeeping)
    if 'germ_plasm_sense' in results and 'housekeeping_sense' in results:
        print("\nComparing germ plasm vs housekeeping...")
        gene_set_comparison = compare_gene_sets(
            results['germ_plasm_sense'],
            results['housekeeping_sense'],
            te_info
        )

        print(f"  Families unique to germ plasm: {len(gene_set_comparison['germ_plasm_only'])}")
        print(f"  Families unique to housekeeping: {len(gene_set_comparison['housekeeping_only'])}")
        print(f"  Families in both: {len(gene_set_comparison['both'])}")

        # Save enrichment data
        if gene_set_comparison['enrichment']:
            enrichment_df = pd.DataFrame(gene_set_comparison['enrichment'])
            enrichment_df = enrichment_df.sort_values('log2_enrichment', ascending=False)
            enrichment_file = args.output_dir / 'germ_plasm_vs_housekeeping_enrichment.tsv'
            enrichment_df.to_csv(enrichment_file, sep='\t', index=False)
            print(f"  Saved: {enrichment_file.name}")

            # Show top enriched in germ plasm
            print("\n  Top 10 families enriched in germ plasm:")
            for _, row in enrichment_df.head(10).iterrows():
                print(f"    {row['family']}: {row['log2_enrichment']:.2f} log2 fold ({row['class']})")

        # Save unique families
        unique_file = args.output_dir / 'unique_families.json'
        with open(unique_file, 'w') as f:
            json.dump({
                'germ_plasm_only': gene_set_comparison['germ_plasm_only'],
                'housekeeping_only': gene_set_comparison['housekeeping_only']
            }, f, indent=2)
        print(f"  Saved: {unique_file.name}")

    # Compare sense vs antisense for germ plasm
    if 'germ_plasm_sense' in results and 'germ_plasm_antisense' in results:
        print("\nComparing sense vs antisense (germ plasm)...")
        strand_comparison = compare_strands(
            results['germ_plasm_sense'],
            results['germ_plasm_antisense']
        )

        print(f"  Sense total hits: {strand_comparison['sense_total']}")
        print(f"  Antisense total hits: {strand_comparison['antisense_total']}")
        print(f"  Sense families: {strand_comparison['sense_families']}")
        print(f"  Antisense families: {strand_comparison['antisense_families']}")

        # Save strand comparison
        strand_df = pd.DataFrame(strand_comparison['family_comparison'])
        strand_file = args.output_dir / 'sense_vs_antisense_comparison.tsv'
        strand_df.to_csv(strand_file, sep='\t', index=False)
        print(f"  Saved: {strand_file.name}")

    # Summary of families NOT hit
    print("\n" + "=" * 60)
    print("TE Families Analysis Summary")
    print("-" * 60)

    all_te_families = set()
    for te_id, info in te_info.items():
        all_te_families.add(info['name'])

    if 'germ_plasm_sense' in results:
        hit_families = set(results['germ_plasm_sense']['families_hit'].keys())
        not_hit = all_te_families - hit_families

        print(f"\nTotal TE families in database: {len(all_te_families)}")
        print(f"Families with hits (germ plasm sense): {len(hit_families)}")
        print(f"Families without hits: {len(not_hit)}")

        # Save families not hit
        not_hit_file = args.output_dir / 'families_not_hit.txt'
        with open(not_hit_file, 'w') as f:
            for family in sorted(not_hit):
                f.write(f"{family}\n")
        print(f"  Saved not-hit families to: {not_hit_file.name}")

    # Save full results as JSON
    results_file = args.output_dir / 'full_analysis.json'
    # Convert defaultdicts to regular dicts for JSON
    json_results = {}
    for key, value in results.items():
        if isinstance(value, dict):
            json_results[key] = {
                'gene_set': value.get('gene_set'),
                'total_hits': value.get('total_hits'),
                'num_families': len(value.get('families_hit', {})),
                'top_families': [
                    {'family': f, 'count': d['count'], 'class': d.get('class', 'Unknown')}
                    for f, d in sorted(value.get('families_hit', {}).items(),
                                      key=lambda x: -x[1]['count'])[:20]
                ]
            }

    with open(results_file, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"\nSaved full results to: {results_file.name}")

    print("\n" + "=" * 60)
    print("TE family analysis complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
