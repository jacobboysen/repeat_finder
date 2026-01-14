#!/usr/bin/env python3
"""
Compare 5'UTR and 3'UTR TE analysis results.
"""

from pathlib import Path
from collections import defaultdict


def load_te_data(path):
    """Load gene-level TE density data."""
    gene_data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                rank = int(parts[0])
                fbgn = parts[1]
                density = float(parts[2])
                hits = int(parts[3])
                hit_bp = int(parts[4])
                utr_len = int(parts[5])
                gene_data[fbgn] = {
                    'rank': rank,
                    'density': density,
                    'hits': hits,
                    'hit_bp': hit_bp,
                    'utr_len': utr_len
                }
    return gene_data


def load_strand_bias(path):
    """Load strand bias data."""
    utr_data = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                fbtr = parts[0]
                fbgn = parts[1]
                utr_len = int(parts[2])
                total_hits = int(parts[3])
                sense_hits = int(parts[4])
                anti_hits = int(parts[5])
                sense_pct = float(parts[6])

                if fbgn not in utr_data:
                    utr_data[fbgn] = {
                        'total_hits': 0,
                        'sense_hits': 0,
                        'anti_hits': 0,
                        'total_utr_len': 0
                    }
                utr_data[fbgn]['total_hits'] += total_hits
                utr_data[fbgn]['sense_hits'] += sense_hits
                utr_data[fbgn]['anti_hits'] += anti_hits
                utr_data[fbgn]['total_utr_len'] += utr_len

    # Calculate percentages
    for fbgn, data in utr_data.items():
        total = data['total_hits']
        if total > 0:
            data['sense_pct'] = 100 * data['sense_hits'] / total
        else:
            data['sense_pct'] = 50

    return utr_data


def load_germ_plasm_genes(path):
    """Load germ plasm gene list from TSV file."""
    genes = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                symbol = parts[0]
                fbgn = parts[1]
                genes[fbgn] = symbol
    return genes


def main():
    base_dir = Path('/Users/jacobboysen/git_repos/repeat_finder')

    # Load data
    print("Loading data...")
    utr3_top = load_te_data(base_dir / 'results/top_100_te_genes_FIXED.tsv')
    utr3_bottom = load_te_data(base_dir / 'results/bottom_100_te_genes_FIXED.tsv')
    utr3_all = {**utr3_top, **utr3_bottom}

    utr5_top = load_te_data(base_dir / 'results/5utr_analysis/top_100_te_genes_FIXED.tsv')
    utr5_bottom = load_te_data(base_dir / 'results/5utr_analysis/bottom_100_te_genes_FIXED.tsv')
    utr5_all = {**utr5_top, **utr5_bottom}

    utr3_strand = load_strand_bias(base_dir / 'results/strand_bias_by_utr.tsv')
    utr5_strand = load_strand_bias(base_dir / 'results/5utr_analysis/strand_bias_by_utr.tsv')

    # Load germ plasm genes from source of truth file
    germ_plasm_genes = load_germ_plasm_genes(
        base_dir / 'data/gene_lists/germ_plasm_genes_consolidated.tsv'
    )

    # Summary statistics
    print("\n" + "="*80)
    print("5'UTR vs 3'UTR TE CONTENT COMPARISON")
    print("="*80)

    print("\n## Overall Statistics")
    utr5_label = "5'UTR"
    utr3_label = "3'UTR"
    print(f"\n{'Metric':<30} {utr5_label:>15} {utr3_label:>15} {'Ratio':>10}")
    print("-" * 70)

    # Total genes
    total_5utr = 12314  # From analysis output
    total_3utr = 13298
    print(f"{'Genes with hits':<30} {total_5utr:>15,} {total_3utr:>15,} {total_5utr/total_3utr:>10.2f}x")

    # Total hits (from strand bias files)
    total_5utr_hits = sum(d['total_hits'] for d in utr5_strand.values())
    total_3utr_hits = sum(d['total_hits'] for d in utr3_strand.values())
    print(f"{'Total BLAST hits':<30} {total_5utr_hits:>15,} {total_3utr_hits:>15,} {total_5utr_hits/total_3utr_hits:>10.2f}x")

    # Strand bias
    sense_5utr = sum(d['sense_hits'] for d in utr5_strand.values())
    sense_3utr = sum(d['sense_hits'] for d in utr3_strand.values())
    sense_pct_5utr = 100 * sense_5utr / total_5utr_hits if total_5utr_hits > 0 else 0
    sense_pct_3utr = 100 * sense_3utr / total_3utr_hits if total_3utr_hits > 0 else 0
    print(f"{'Sense strand bias':<30} {sense_pct_5utr:>14.1f}% {sense_pct_3utr:>14.1f}% {sense_pct_5utr/sense_pct_3utr:>10.2f}x")

    # Top gene densities
    if utr5_top and utr3_top:
        max_5utr = max(d['density'] for d in utr5_top.values())
        max_3utr = max(d['density'] for d in utr3_top.values())
        print(f"{'Max TE density (bp/kb)':<30} {max_5utr:>15,.0f} {max_3utr:>15,.0f} {max_5utr/max_3utr:>10.2f}x")

    # Average UTR lengths
    avg_5utr_len = sum(d['total_utr_len'] for d in utr5_strand.values()) / len(utr5_strand) if utr5_strand else 0
    avg_3utr_len = sum(d['total_utr_len'] for d in utr3_strand.values()) / len(utr3_strand) if utr3_strand else 0
    print(f"{'Avg UTR length (bp)':<30} {avg_5utr_len:>15,.0f} {avg_3utr_len:>15,.0f} {avg_5utr_len/avg_3utr_len:>10.2f}x")

    # Top 10 comparison
    print("\n" + "="*80)
    print("TOP 10 GENES BY TE DENSITY")
    print("="*80)

    col_5utr = utr5_label + " Gene"
    col_3utr = utr3_label + " Gene"
    print(f"\n{'Rank':<6} {col_5utr:<15} {'Density':>12} {col_3utr:<15} {'Density':>12}")
    print("-" * 70)

    utr5_sorted = sorted(utr5_top.items(), key=lambda x: x[1]['rank'])[:10]
    utr3_sorted = sorted(utr3_top.items(), key=lambda x: x[1]['rank'])[:10]

    for i in range(10):
        g5, d5 = utr5_sorted[i] if i < len(utr5_sorted) else ('—', {'density': 0})
        g3, d3 = utr3_sorted[i] if i < len(utr3_sorted) else ('—', {'density': 0})
        print(f"{i+1:<6} {g5:<15} {d5['density']:>12,.0f} {g3:<15} {d3['density']:>12,.0f}")

    # Overlap analysis
    print("\n" + "="*80)
    print("OVERLAP ANALYSIS")
    print("="*80)

    top100_5utr = set(utr5_top.keys())
    top100_3utr = set(utr3_top.keys())
    overlap_top = top100_5utr & top100_3utr

    print(f"\nTop 100 genes overlap: {len(overlap_top)} genes")
    if overlap_top:
        print("Genes in BOTH top 100 lists:")
        for fbgn in sorted(overlap_top):
            r5 = utr5_top[fbgn]['rank']
            r3 = utr3_top[fbgn]['rank']
            print(f"  {fbgn}: 5'UTR rank #{r5}, 3'UTR rank #{r3}")

    # Bottom overlap
    bot100_5utr = set(utr5_bottom.keys())
    bot100_3utr = set(utr3_bottom.keys())
    overlap_bot = bot100_5utr & bot100_3utr

    print(f"\nBottom 100 genes overlap: {len(overlap_bot)} genes")

    # Genes in one top and other bottom
    top5_bot3 = top100_5utr & bot100_3utr
    top3_bot5 = top100_3utr & bot100_5utr

    if top5_bot3:
        print(f"\nGenes HIGH in 5'UTR but LOW in 3'UTR: {len(top5_bot3)}")
        for fbgn in sorted(top5_bot3)[:5]:
            print(f"  {fbgn}")

    if top3_bot5:
        print(f"\nGenes HIGH in 3'UTR but LOW in 5'UTR: {len(top3_bot5)}")
        for fbgn in sorted(top3_bot5)[:5]:
            print(f"  {fbgn}")

    # Germ plasm gene comparison
    print("\n" + "="*80)
    print("GERM PLASM GENES COMPARISON")
    print("="*80)

    # germ_plasm_genes loaded from file at start of main()

    h1 = utr5_label + " Hits"
    h2 = utr5_label + " %Sense"
    h3 = utr3_label + " Hits"
    h4 = utr3_label + " %Sense"
    print(f"\n{'Gene':<8} {h1:>12} {h2:>14} {h3:>12} {h4:>14}")
    print("-" * 70)

    for fbgn, symbol in germ_plasm_genes.items():
        d5 = utr5_strand.get(fbgn, {'total_hits': 0, 'sense_pct': 0})
        d3 = utr3_strand.get(fbgn, {'total_hits': 0, 'sense_pct': 0})
        print(f"{symbol:<8} {d5['total_hits']:>12,} {d5['sense_pct']:>13.1f}% {d3['total_hits']:>12,} {d3['sense_pct']:>13.1f}%")

    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    print("""
1. 5'UTRs have FEWER total TE hits (1.9M vs 2.6M for 3'UTRs)
2. 5'UTRs show STRONGER sense strand bias (65.8% vs 60.2%)
3. 5'UTR top genes have HIGHER TE density (~3x higher than 3'UTR top genes)
4. Top gene lists are almost completely DIFFERENT between 5'UTR and 3'UTR
5. 5'UTRs are shorter on average than 3'UTRs
""")


if __name__ == '__main__':
    main()
