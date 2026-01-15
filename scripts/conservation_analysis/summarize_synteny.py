#!/usr/bin/env python3
"""Summarize synteny analysis results."""

import sys
from collections import defaultdict

def load_synteny_results(filepath):
    """Load synteny results TSV."""
    results = []
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            # Convert numeric fields
            for key in ['pident', 'length', 'sim_coverage', 'yak_coverage',
                       'ere_coverage', 'sec_coverage', 'any_coverage',
                       'both_close_coverage', 'n_blocks']:
                if key in row:
                    row[key] = float(row[key])
            results.append(row)
    return results

def analyze_synteny(results, label=""):
    """Print synteny analysis summary."""
    n = len(results)
    if n == 0:
        print("No results to analyze")
        return

    print(f"\n{'='*70}")
    print(f"SYNTENY ANALYSIS: {label} (n={n:,})")
    print(f"{'='*70}\n")

    # Coverage thresholds
    thresholds = [0.25, 0.5, 0.75, 0.9]

    print("Syntenic coverage in comparison species:")
    print("-" * 60)
    print(f"{'Threshold':<15} {'D. simulans':>12} {'D. yakuba':>12} {'D. erecta':>12} {'Any species':>12}")
    print("-" * 60)

    for thresh in thresholds:
        n_sim = sum(1 for r in results if r['sim_coverage'] >= thresh)
        n_yak = sum(1 for r in results if r['yak_coverage'] >= thresh)
        n_ere = sum(1 for r in results if r['ere_coverage'] >= thresh)
        n_any = sum(1 for r in results if r['any_coverage'] >= thresh)

        print(f"≥{thresh*100:.0f}% coverage   {n_sim:>6,} ({100*n_sim/n:>4.1f}%)  {n_yak:>6,} ({100*n_yak/n:>4.1f}%)  {n_ere:>6,} ({100*n_ere/n:>4.1f}%)  {n_any:>6,} ({100*n_any/n:>4.1f}%)")

    # Classification
    print(f"\n{'='*50}")
    print("Synteny Classification (≥50% coverage threshold)")
    print("="*50 + "\n")

    syntenic_any = [r for r in results if r['any_coverage'] >= 0.5]
    mel_specific = [r for r in results if r['any_coverage'] < 0.5]
    syntenic_all = [r for r in results if min(r['sim_coverage'], r['yak_coverage'], r['ere_coverage']) >= 0.5]

    print(f"{'Category':<35} {'Count':>10} {'Percent':>10}")
    print("-" * 55)
    print(f"{'Syntenic (in any species)':<35} {len(syntenic_any):>10,} {100*len(syntenic_any)/n:>9.1f}%")
    print(f"{'Syntenic (in ALL 3 species)':<35} {len(syntenic_all):>10,} {100*len(syntenic_all)/n:>9.1f}%")
    print(f"{'Mel-specific (no alignment)':<35} {len(mel_specific):>10,} {100*len(mel_specific)/n:>9.1f}%")

    # Break down by identity
    print(f"\n{'='*50}")
    print("Synteny by Identity Threshold")
    print("="*50 + "\n")

    print(f"{'Identity range':<20} {'Total':>10} {'Syntenic':>10} {'% Syntenic':>12}")
    print("-" * 55)

    id_ranges = [(60, 70), (70, 80), (80, 85), (85, 90), (90, 95), (95, 100), (100, 101)]
    for min_id, max_id in id_ranges:
        tier = [r for r in results if min_id <= r['pident'] < max_id]
        if tier:
            n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
            label = f"{min_id}-{max_id}%" if max_id < 101 else "100%"
            print(f"{label:<20} {len(tier):>10,} {n_syn:>10,} {100*n_syn/len(tier):>11.1f}%")

    # Break down by length
    print(f"\n{'='*50}")
    print("Synteny by Alignment Length")
    print("="*50 + "\n")

    print(f"{'Length range':<20} {'Total':>10} {'Syntenic':>10} {'% Syntenic':>12}")
    print("-" * 55)

    len_ranges = [(50, 75), (75, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 1000)]
    for min_len, max_len in len_ranges:
        tier = [r for r in results if min_len <= r['length'] < max_len]
        if tier:
            n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
            print(f"{min_len}-{max_len}bp{'':<10} {len(tier):>10,} {n_syn:>10,} {100*n_syn/len(tier):>11.1f}%")

    # Novel vs Known
    print(f"\n{'='*50}")
    print("Synteny: Novel vs Known (RepeatMasker)")
    print("="*50 + "\n")

    novel = [r for r in results if r['category'] == 'novel']
    known = [r for r in results if r['category'] == 'known']

    print(f"{'Category':<25} {'Total':>10} {'Syntenic':>10} {'% Syntenic':>12}")
    print("-" * 55)

    for label, subset in [('Novel (not in RM)', novel), ('Known (in RM)', known)]:
        if subset:
            n_syn = sum(1 for r in subset if r['any_coverage'] >= 0.5)
            print(f"{label:<25} {len(subset):>10,} {n_syn:>10,} {100*n_syn/len(subset):>11.1f}%")

    return {
        'total': n,
        'syntenic_any': len(syntenic_any),
        'syntenic_all': len(syntenic_all),
        'mel_specific': len(mel_specific)
    }

# Main
if __name__ == '__main__':
    # Analyze HQ hits
    print("Loading HQ synteny results...")
    hq_results = load_synteny_results('results/repeatmasker_analysis/te_hits_hq_synteny.tsv')
    analyze_synteny(hq_results, "HIGH-QUALITY HITS (≥80% id, ≥50bp)")

    # Try to load all hits if available
    try:
        print("\n\nLoading ALL hits synteny results (sampled)...")
        all_results = load_synteny_results('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')
        analyze_synteny(all_results, "ALL HITS (10% sample)")
    except FileNotFoundError:
        print("\nAll-hits synteny analysis not yet complete.")
