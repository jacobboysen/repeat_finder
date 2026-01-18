#!/usr/bin/env python3
"""
Integrated TE enrichment analysis across exons, 5'UTR, and 3'UTR regions.

This script:
1. Loads TE hit data from all three region types
2. Aggregates by gene to find genes with hits in multiple regions
3. Performs functional enrichment analysis using gene sets
4. Identifies top TE-hit genes and interesting patterns
5. Generates comprehensive summary and visualizations
"""

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

# Add scripts directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_results_dir, get_references_dir
from utils.blast_io import load_blast_results


def load_gene_set(path):
    """Load gene set from FBgn ID file."""
    genes = set()
    with open(path) as f:
        for line in f:
            gene = line.strip()
            if gene.startswith('FBgn'):
                genes.add(gene)
    return genes


def load_all_gene_sets(gene_sets_dir):
    """Load all functional gene sets."""
    gene_sets = {}
    gene_sets_dir = Path(gene_sets_dir)

    for f in gene_sets_dir.glob('*_fbgn_ids.txt'):
        name = f.stem.replace('_fbgn_ids', '')
        gene_sets[name] = load_gene_set(f)

    return gene_sets


def parse_qseqid_for_gene(qseqid):
    """Extract FBgn from various qseqid formats."""
    # Format: FBtr0070000 or FBtr0302344_exon1
    if '_exon' in qseqid:
        # Exon format - need to look up gene
        return None
    return None  # UTR format uses transcript ID


def aggregate_utr_hits_by_gene(blast_file, fbtr_to_fbgn):
    """Aggregate UTR BLAST hits by gene (RAW - use load_deduplicated_utr_stats instead)."""
    gene_stats = defaultdict(lambda: {'hits': 0, 'hit_bp': 0, 'best_pident': 0, 'best_bitscore': 0})

    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            qseqid = parts[0]
            pident = float(parts[2])
            length = int(parts[3])
            bitscore = float(parts[11])

            # Get gene ID
            fbtr = qseqid.split('_')[0] if '_' in qseqid else qseqid
            fbgn = fbtr_to_fbgn.get(fbtr)

            if fbgn:
                gene_stats[fbgn]['hits'] += 1
                gene_stats[fbgn]['hit_bp'] += length
                gene_stats[fbgn]['best_pident'] = max(gene_stats[fbgn]['best_pident'], pident)
                gene_stats[fbgn]['best_bitscore'] = max(gene_stats[fbgn]['best_bitscore'], bitscore)

    return dict(gene_stats)


def load_deduplicated_utr_stats(dedup_file):
    """Load deduplicated UTR stats from *_original_vs_deduplicated.tsv file."""
    gene_stats = {}

    with open(dedup_file) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            fbgn = parts[0]
            unique_hits = int(parts[2])

            gene_stats[fbgn] = {
                'hits': unique_hits,
                'hit_bp': 0,  # Not tracked in dedup summary
                'best_pident': 0,  # Not tracked in dedup summary
            }

    return gene_stats


def load_exon_gene_stats(exon_summary_file, deduplicated=True):
    """Load pre-aggregated exon TE stats by gene.

    Args:
        exon_summary_file: Path to gene summary TSV
        deduplicated: If True, expects deduplicated format; if False, original format
    """
    gene_stats = {}

    with open(exon_summary_file) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            fbgn = parts[1]

            if deduplicated:
                # Deduplicated format: rank, fbgn, gene_symbol, unique_hits, hit_bp, sense_hits, ...
                gene_stats[fbgn] = {
                    'hits': int(parts[3]),
                    'hit_bp': int(parts[4]),
                    'sense_pct': float(parts[7]) if len(parts) > 7 else 50.0,
                }
            else:
                # Original format: rank, fbgn, gene_symbol, exons, exons_with_hits, total_hits, ...
                gene_stats[fbgn] = {
                    'hits': int(parts[5]),
                    'hit_bp': int(parts[6]),
                    'density': float(parts[8]) if len(parts) > 8 else 0,
                    'sense_pct': float(parts[11]) if len(parts) > 11 else 50.0,
                }

    return gene_stats


def load_fbtr_to_fbgn_mapping(gff_file):
    """Build transcript to gene mapping from GFF."""
    import re
    mapping = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'mRNA':
                attrs = parts[8]
                # Extract ID (FBtr) and Parent (FBgn)
                id_match = re.search(r'ID=([^;]+)', attrs)
                parent_match = re.search(r'Parent=([^;]+)', attrs)

                if id_match and parent_match:
                    fbtr = id_match.group(1)
                    fbgn = parent_match.group(1)
                    mapping[fbtr] = fbgn

    return mapping


def load_fbgn_to_symbol(ref_dir):
    """Load FBgn to symbol mapping."""
    mapping = {}
    symbol_file = Path(ref_dir) / 'fbgn_to_symbol.tsv'

    if symbol_file.exists():
        with open(symbol_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2 and parts[0].startswith('FBgn'):
                    mapping[parts[0]] = parts[1]

    return mapping


def fisher_exact_test(a, b, c, d):
    """
    Simple Fisher's exact test approximation.

    Contingency table:
              In Set   Not In Set
    High TE     a          b
    Low TE      c          d
    """
    n = a + b + c + d
    if n == 0:
        return 1.0

    # Calculate odds ratio
    if b == 0 or c == 0:
        odds_ratio = float('inf') if a * d > 0 else 0
    else:
        odds_ratio = (a * d) / (b * c)

    # Use chi-square approximation for p-value
    expected_a = (a + b) * (a + c) / n
    if expected_a == 0:
        return 1.0

    chi2 = ((a - expected_a) ** 2) / expected_a

    # Very rough p-value approximation
    # For chi2 with 1 df, p < 0.05 when chi2 > 3.84
    if chi2 < 0.001:
        return 1.0
    elif chi2 < 1:
        return 0.5
    elif chi2 < 3.84:
        return 0.1
    elif chi2 < 6.63:
        return 0.05
    elif chi2 < 10.83:
        return 0.01
    else:
        return 0.001


def enrichment_analysis(gene_te_data, gene_sets, threshold_pct=10):
    """
    Perform enrichment analysis for gene sets.

    Tests if genes in each set have higher/lower TE content than background.
    """
    # Get all genes with TE data
    all_genes = set(gene_te_data.keys())

    # Calculate threshold for "high TE" (top X%)
    all_hits = sorted([g['total_hits'] for g in gene_te_data.values()], reverse=True)
    if len(all_hits) > 0:
        threshold_idx = int(len(all_hits) * threshold_pct / 100)
        hit_threshold = all_hits[threshold_idx] if threshold_idx < len(all_hits) else 0
    else:
        hit_threshold = 0

    high_te_genes = {g for g, d in gene_te_data.items() if d['total_hits'] >= hit_threshold}
    low_te_genes = all_genes - high_te_genes

    results = []

    for set_name, set_genes in gene_sets.items():
        # Genes in set that we have TE data for
        set_genes_with_data = set_genes & all_genes

        if len(set_genes_with_data) < 10:
            continue  # Skip small sets

        # Contingency table
        a = len(set_genes_with_data & high_te_genes)  # In set, high TE
        b = len(high_te_genes - set_genes_with_data)  # Not in set, high TE
        c = len(set_genes_with_data & low_te_genes)   # In set, low TE
        d = len(low_te_genes - set_genes_with_data)   # Not in set, low TE

        # Calculate enrichment
        set_high_pct = 100 * a / len(set_genes_with_data) if len(set_genes_with_data) > 0 else 0
        bg_high_pct = 100 * len(high_te_genes) / len(all_genes) if len(all_genes) > 0 else 0

        fold_enrichment = set_high_pct / bg_high_pct if bg_high_pct > 0 else 0

        # Calculate mean TE hits for set vs background
        set_mean_hits = sum(gene_te_data[g]['total_hits'] for g in set_genes_with_data) / len(set_genes_with_data)
        bg_mean_hits = sum(d['total_hits'] for d in gene_te_data.values()) / len(gene_te_data)

        p_value = fisher_exact_test(a, b, c, d)

        results.append({
            'set_name': set_name,
            'set_size': len(set_genes_with_data),
            'high_te_in_set': a,
            'set_high_pct': set_high_pct,
            'bg_high_pct': bg_high_pct,
            'fold_enrichment': fold_enrichment,
            'set_mean_hits': set_mean_hits,
            'bg_mean_hits': bg_mean_hits,
            'mean_fold': set_mean_hits / bg_mean_hits if bg_mean_hits > 0 else 0,
            'p_value': p_value,
        })

    return sorted(results, key=lambda x: -x['fold_enrichment'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        default=get_results_dir() / 'integrated_te_analysis')
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Integrated TE Enrichment Analysis")
    print("=" * 70)

    # Load transcript-to-gene mapping
    print("\nLoading reference data...")
    gff_file = get_references_dir() / 'dmel-all-r6.66.gff'
    fbtr_to_fbgn = load_fbtr_to_fbgn_mapping(gff_file)
    print(f"  Transcript-to-gene mapping: {len(fbtr_to_fbgn):,} transcripts")

    fbgn_to_symbol = load_fbgn_to_symbol(get_references_dir())
    print(f"  Gene symbols: {len(fbgn_to_symbol):,} genes")

    # Load functional gene sets
    print("\nLoading functional gene sets...")
    gene_sets_dir = Path('data/gene_lists/functional')
    gene_sets = load_all_gene_sets(gene_sets_dir)
    print(f"  Loaded {len(gene_sets)} gene sets")

    # Load TE data from all three regions
    print("\n" + "=" * 70)
    print("Loading TE hit data from all regions (DEDUPLICATED)...")
    print("=" * 70)

    # 3'UTR data (deduplicated)
    print("\n  Loading 3'UTR data (deduplicated)...")
    utr3_dedup_file = get_results_dir() / '3utr_deduplicated' / '3utr_original_vs_deduplicated.tsv'
    if utr3_dedup_file.exists():
        utr3_stats = load_deduplicated_utr_stats(utr3_dedup_file)
        print(f"    3'UTR: {len(utr3_stats):,} genes with hits (deduplicated)")
    else:
        print("    [Warning: Using non-deduplicated data - run deduplicate_te_hits.py --seq-type 3utr first]")
        utr3_file = get_results_dir() / 'genome_wide_all_3utrs.tsv'
        utr3_stats = aggregate_utr_hits_by_gene(utr3_file, fbtr_to_fbgn)
        print(f"    3'UTR: {len(utr3_stats):,} genes with hits (RAW - inflated)")

    # 5'UTR data (deduplicated)
    print("  Loading 5'UTR data (deduplicated)...")
    utr5_dedup_file = get_results_dir() / '5utr_deduplicated' / '5utr_original_vs_deduplicated.tsv'
    if utr5_dedup_file.exists():
        utr5_stats = load_deduplicated_utr_stats(utr5_dedup_file)
        print(f"    5'UTR: {len(utr5_stats):,} genes with hits (deduplicated)")
    else:
        print("    [Warning: Using non-deduplicated data - run deduplicate_te_hits.py --seq-type 5utr first]")
        utr5_file = get_results_dir() / 'genome_wide_all_5utrs.tsv'
        utr5_stats = aggregate_utr_hits_by_gene(utr5_file, fbtr_to_fbgn)
        print(f"    5'UTR: {len(utr5_stats):,} genes with hits (RAW - inflated)")

    # Exon data (use deduplicated counts)
    print("  Loading exon data (deduplicated)...")
    exon_file = get_results_dir() / 'exon_analysis' / 'deduplicated' / 'gene_te_summary_deduplicated.tsv'
    if not exon_file.exists():
        # Fall back to original if deduplicated doesn't exist
        print("    [Warning: Using non-deduplicated data - run deduplicate_exon_te_hits.py first]")
        exon_file = get_results_dir() / 'exon_analysis' / 'gene_exon_te_summary.tsv'
        exon_stats = load_exon_gene_stats(exon_file, deduplicated=False)
    else:
        exon_stats = load_exon_gene_stats(exon_file, deduplicated=True)
    print(f"    Exons: {len(exon_stats):,} genes with hits")

    # Combine data by gene
    print("\n  Integrating data across regions...")
    all_genes = set(utr3_stats.keys()) | set(utr5_stats.keys()) | set(exon_stats.keys())

    gene_te_data = {}
    for fbgn in all_genes:
        utr3 = utr3_stats.get(fbgn, {'hits': 0, 'hit_bp': 0, 'best_pident': 0})
        utr5 = utr5_stats.get(fbgn, {'hits': 0, 'hit_bp': 0, 'best_pident': 0})
        exon = exon_stats.get(fbgn, {'hits': 0, 'hit_bp': 0, 'density': 0})

        gene_te_data[fbgn] = {
            'symbol': fbgn_to_symbol.get(fbgn, fbgn),
            'utr3_hits': utr3['hits'],
            'utr5_hits': utr5['hits'],
            'exon_hits': exon['hits'],
            'total_hits': utr3['hits'] + utr5['hits'] + exon['hits'],
            'utr3_bp': utr3['hit_bp'],
            'utr5_bp': utr5['hit_bp'],
            'exon_bp': exon['hit_bp'],
            'total_bp': utr3['hit_bp'] + utr5['hit_bp'] + exon['hit_bp'],
            'regions_with_hits': sum([utr3['hits'] > 0, utr5['hits'] > 0, exon['hits'] > 0]),
            'best_3utr_pident': utr3.get('best_pident', 0),
            'best_5utr_pident': utr5.get('best_pident', 0),
        }

    print(f"    Total genes with any TE hits: {len(gene_te_data):,}")

    # Count genes by region combination
    region_combos = defaultdict(int)
    for g, d in gene_te_data.items():
        has_3utr = d['utr3_hits'] > 0
        has_5utr = d['utr5_hits'] > 0
        has_exon = d['exon_hits'] > 0

        combo = []
        if has_3utr: combo.append("3'UTR")
        if has_5utr: combo.append("5'UTR")
        if has_exon: combo.append("Exon")

        if combo:
            region_combos[tuple(combo)] += 1

    print("\n  Genes by region combination:")
    for combo, count in sorted(region_combos.items(), key=lambda x: -x[1]):
        print(f"    {' + '.join(combo)}: {count:,}")

    # ========================================
    # TOP GENES BY TOTAL TE HITS
    # ========================================
    print("\n" + "=" * 70)
    print("TOP 50 GENES BY TOTAL TE HITS (ALL REGIONS)")
    print("=" * 70)

    top_genes = sorted(gene_te_data.items(), key=lambda x: -x[1]['total_hits'])[:50]

    print(f"\n{'Rank':<5} {'Gene':<15} {'Total':>8} {'3UTR':>8} {'5UTR':>8} {'Exon':>8} {'Regions':>8}")
    print("-" * 70)
    for i, (fbgn, d) in enumerate(top_genes, 1):
        print(f"{i:<5} {d['symbol']:<15} {d['total_hits']:>8,} {d['utr3_hits']:>8,} "
              f"{d['utr5_hits']:>8,} {d['exon_hits']:>8,} {d['regions_with_hits']:>8}")

    # ========================================
    # GENES WITH HITS IN ALL 3 REGIONS
    # ========================================
    print("\n" + "=" * 70)
    print("TOP 30 GENES WITH HITS IN ALL 3 REGIONS")
    print("=" * 70)

    all_region_genes = [(g, d) for g, d in gene_te_data.items() if d['regions_with_hits'] == 3]
    all_region_genes = sorted(all_region_genes, key=lambda x: -x[1]['total_hits'])[:30]

    print(f"\n{'Rank':<5} {'Gene':<15} {'Total':>8} {'3UTR':>8} {'5UTR':>8} {'Exon':>8}")
    print("-" * 60)
    for i, (fbgn, d) in enumerate(all_region_genes, 1):
        print(f"{i:<5} {d['symbol']:<15} {d['total_hits']:>8,} {d['utr3_hits']:>8,} "
              f"{d['utr5_hits']:>8,} {d['exon_hits']:>8,}")

    # ========================================
    # FUNCTIONAL ENRICHMENT ANALYSIS
    # ========================================
    print("\n" + "=" * 70)
    print("FUNCTIONAL ENRICHMENT ANALYSIS")
    print("=" * 70)

    enrichment_results = enrichment_analysis(gene_te_data, gene_sets, threshold_pct=10)

    # Top enriched sets
    print("\nGene sets ENRICHED for high TE content (top 20):")
    print(f"{'Set Name':<40} {'Size':>6} {'High%':>7} {'BG%':>6} {'Fold':>6} {'Mean Hits':>10}")
    print("-" * 85)

    enriched = [r for r in enrichment_results if r['fold_enrichment'] > 1.0][:20]
    for r in enriched:
        print(f"{r['set_name']:<40} {r['set_size']:>6} {r['set_high_pct']:>6.1f}% "
              f"{r['bg_high_pct']:>5.1f}% {r['fold_enrichment']:>5.2f}x {r['set_mean_hits']:>10.1f}")

    # Depleted sets
    print("\nGene sets DEPLETED for high TE content (bottom 20):")
    print(f"{'Set Name':<40} {'Size':>6} {'High%':>7} {'BG%':>6} {'Fold':>6} {'Mean Hits':>10}")
    print("-" * 85)

    depleted = sorted([r for r in enrichment_results if r['fold_enrichment'] < 1.0 and r['fold_enrichment'] > 0],
                      key=lambda x: x['fold_enrichment'])[:20]
    for r in depleted:
        print(f"{r['set_name']:<40} {r['set_size']:>6} {r['set_high_pct']:>6.1f}% "
              f"{r['bg_high_pct']:>5.1f}% {r['fold_enrichment']:>5.2f}x {r['set_mean_hits']:>10.1f}")

    # ========================================
    # REGION-SPECIFIC PATTERNS
    # ========================================
    print("\n" + "=" * 70)
    print("REGION-SPECIFIC PATTERNS")
    print("=" * 70)

    # Genes with high 3'UTR but low exon hits (classic UTR TE insertions)
    print("\nGenes with HIGH 3'UTR hits but LOW exon hits (UTR-specific):")
    utr_specific = [(g, d) for g, d in gene_te_data.items()
                   if d['utr3_hits'] > 100 and d['exon_hits'] < 10]
    utr_specific = sorted(utr_specific, key=lambda x: -x[1]['utr3_hits'])[:15]

    print(f"{'Gene':<15} {'3UTR':>8} {'5UTR':>8} {'Exon':>8}")
    print("-" * 45)
    for fbgn, d in utr_specific:
        print(f"{d['symbol']:<15} {d['utr3_hits']:>8,} {d['utr5_hits']:>8,} {d['exon_hits']:>8,}")

    # Genes with high exon hits but low UTR hits (exon-specific)
    print("\nGenes with HIGH exon hits but LOW UTR hits (exon-specific):")
    exon_specific = [(g, d) for g, d in gene_te_data.items()
                    if d['exon_hits'] > 100 and d['utr3_hits'] < 10 and d['utr5_hits'] < 10]
    exon_specific = sorted(exon_specific, key=lambda x: -x[1]['exon_hits'])[:15]

    print(f"{'Gene':<15} {'Exon':>8} {'3UTR':>8} {'5UTR':>8}")
    print("-" * 45)
    for fbgn, d in exon_specific:
        print(f"{d['symbol']:<15} {d['exon_hits']:>8,} {d['utr3_hits']:>8,} {d['utr5_hits']:>8,}")

    # Genes with balanced hits across all regions
    print("\nGenes with BALANCED hits across all regions:")
    balanced = []
    for g, d in gene_te_data.items():
        if d['total_hits'] > 50 and d['regions_with_hits'] == 3:
            min_hits = min(d['utr3_hits'], d['utr5_hits'], d['exon_hits'])
            max_hits = max(d['utr3_hits'], d['utr5_hits'], d['exon_hits'])
            if max_hits > 0 and min_hits / max_hits > 0.2:  # Ratio > 0.2
                balanced.append((g, d, min_hits / max_hits))

    balanced = sorted(balanced, key=lambda x: -x[2])[:15]
    print(f"{'Gene':<15} {'3UTR':>8} {'5UTR':>8} {'Exon':>8} {'Balance':>8}")
    print("-" * 55)
    for fbgn, d, ratio in balanced:
        print(f"{d['symbol']:<15} {d['utr3_hits']:>8,} {d['utr5_hits']:>8,} "
              f"{d['exon_hits']:>8,} {ratio:>7.2f}")

    # ========================================
    # SAVE RESULTS
    # ========================================
    print("\n" + "=" * 70)
    print("Saving results...")
    print("=" * 70)

    # Save full gene TE data
    gene_file = args.output_dir / 'gene_te_all_regions.tsv'
    with open(gene_file, 'w') as f:
        f.write('fbgn\tsymbol\ttotal_hits\tutr3_hits\tutr5_hits\texon_hits\t'
                'total_bp\tutr3_bp\tutr5_bp\texon_bp\tregions_with_hits\n')
        for fbgn, d in sorted(gene_te_data.items(), key=lambda x: -x[1]['total_hits']):
            f.write(f"{fbgn}\t{d['symbol']}\t{d['total_hits']}\t{d['utr3_hits']}\t"
                    f"{d['utr5_hits']}\t{d['exon_hits']}\t{d['total_bp']}\t"
                    f"{d['utr3_bp']}\t{d['utr5_bp']}\t{d['exon_bp']}\t{d['regions_with_hits']}\n")
    print(f"  Saved gene data: {gene_file}")

    # Save enrichment results
    enrich_file = args.output_dir / 'functional_enrichment_results.tsv'
    with open(enrich_file, 'w') as f:
        f.write('set_name\tset_size\thigh_te_in_set\tset_high_pct\tbg_high_pct\t'
                'fold_enrichment\tset_mean_hits\tbg_mean_hits\tmean_fold\tp_value\n')
        for r in enrichment_results:
            f.write(f"{r['set_name']}\t{r['set_size']}\t{r['high_te_in_set']}\t"
                    f"{r['set_high_pct']:.2f}\t{r['bg_high_pct']:.2f}\t"
                    f"{r['fold_enrichment']:.3f}\t{r['set_mean_hits']:.1f}\t"
                    f"{r['bg_mean_hits']:.1f}\t{r['mean_fold']:.3f}\t{r['p_value']}\n")
    print(f"  Saved enrichment: {enrich_file}")

    # Save summary JSON
    summary = {
        'total_genes': len(gene_te_data),
        'genes_3utr_only': len([g for g, d in gene_te_data.items() if d['utr3_hits'] > 0 and d['utr5_hits'] == 0 and d['exon_hits'] == 0]),
        'genes_5utr_only': len([g for g, d in gene_te_data.items() if d['utr5_hits'] > 0 and d['utr3_hits'] == 0 and d['exon_hits'] == 0]),
        'genes_exon_only': len([g for g, d in gene_te_data.items() if d['exon_hits'] > 0 and d['utr3_hits'] == 0 and d['utr5_hits'] == 0]),
        'genes_all_regions': len([g for g, d in gene_te_data.items() if d['regions_with_hits'] == 3]),
        'region_combinations': {' + '.join(k): v for k, v in region_combos.items()},
        'top_enriched': [r['set_name'] for r in enriched[:5]],
        'top_depleted': [r['set_name'] for r in depleted[:5]],
    }

    summary_file = args.output_dir / 'analysis_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"  Saved summary: {summary_file}")

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
