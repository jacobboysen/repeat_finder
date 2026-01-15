#!/usr/bin/env python3
"""Run fast synteny analysis on HQ hits only for comparison."""

import gzip
import sys
import os
from collections import defaultdict

def extract_species_coverage_from_maf(maf_file, chrom):
    """Extract coverage intervals for each species from MAF file."""
    species_intervals = defaultdict(list)
    species_of_interest = ['droSim1', 'droYak2', 'droYak3', 'droEre2', 'droSec1']

    if maf_file.endswith('.gz'):
        fh = gzip.open(maf_file, 'rt')
    else:
        fh = open(maf_file)

    current_ref_start = None
    current_ref_end = None
    current_species = set()

    for line in fh:
        if line.startswith('a '):
            if current_ref_start is not None:
                for sp in current_species:
                    species_intervals[sp].append((current_ref_start, current_ref_end))
            current_ref_start = None
            current_ref_end = None
            current_species = set()

        elif line.startswith('s '):
            parts = line.strip().split()
            species_chrom = parts[1]
            start = int(parts[2])
            length = int(parts[3])

            if '.' in species_chrom:
                species = species_chrom.split('.')[0]
            else:
                species = species_chrom

            if species == 'dm6':
                current_ref_start = start
                current_ref_end = start + length
            elif species in species_of_interest:
                current_species.add(species)

    if current_ref_start is not None:
        for sp in current_species:
            species_intervals[sp].append((current_ref_start, current_ref_end))

    fh.close()

    for sp in species_intervals:
        species_intervals[sp] = merge_intervals(sorted(species_intervals[sp]))

    return dict(species_intervals)


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def calculate_coverage(query_start, query_end, intervals):
    if not intervals:
        return 0.0

    query_len = query_end - query_start
    if query_len <= 0:
        return 0.0

    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= query_start:
            lo = mid + 1
        else:
            hi = mid

    covered = 0
    for i in range(lo, len(intervals)):
        int_start, int_end = intervals[i]
        if int_start >= query_end:
            break
        overlap_start = max(query_start, int_start)
        overlap_end = min(query_end, int_end)
        if overlap_end > overlap_start:
            covered += overlap_end - overlap_start

    return covered / query_len


def main():
    maf_dir = '/Users/jacobboysen/git_repos/repeat_finder/data/references/maf'
    results_dir = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis'

    # Load HQ hits
    print("Loading HQ hits...", file=sys.stderr)
    hits_by_chrom = defaultdict(list)

    with open(os.path.join(results_dir, 'te_hits_hq_for_synteny.bed')) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]
            name_parts = name.split('|')
            category = name_parts[0].split('_')[0]
            pident = float(name_parts[3])
            length = int(name_parts[4])
            hits_by_chrom[chrom].append((start, end, name, pident, length, category))

    total = sum(len(h) for h in hits_by_chrom.values())
    print(f"Loaded {total:,} HQ hits", file=sys.stderr)

    # Process each chromosome
    results = []
    for chrom in ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']:
        if chrom not in hits_by_chrom:
            continue
        maf_file = os.path.join(maf_dir, f'{chrom}.maf.gz')
        if not os.path.exists(maf_file):
            continue

        print(f"Processing {chrom}...", file=sys.stderr)
        coverage = extract_species_coverage_from_maf(maf_file, chrom)

        for start, end, name, pident, length, category in hits_by_chrom[chrom]:
            sim_cov = calculate_coverage(start, end, coverage.get('droSim1', []))
            sec_cov = calculate_coverage(start, end, coverage.get('droSec1', []))
            yak_cov = max(
                calculate_coverage(start, end, coverage.get('droYak2', [])),
                calculate_coverage(start, end, coverage.get('droYak3', []))
            )
            ere_cov = calculate_coverage(start, end, coverage.get('droEre2', []))
            any_cov = max(sim_cov, sec_cov, yak_cov, ere_cov)

            results.append({
                'name': name, 'pident': pident, 'length': length, 'category': category,
                'sim': sim_cov, 'yak': yak_cov, 'ere': ere_cov, 'sec': sec_cov, 'any': any_cov
            })

    # Summary
    print("\n" + "="*70)
    print("FAST SYNTENY: HQ HITS SUMMARY")
    print("="*70)

    n = len(results)
    print(f"\nTotal HQ hits: {n:,}")

    for thresh in [0.25, 0.5, 0.75]:
        n_any = sum(1 for r in results if r['any'] >= thresh)
        print(f"Syntenic (any species ≥{100*thresh:.0f}%): {n_any:,} ({100*n_any/n:.1f}%)")

    # By identity
    print("\nBy identity tier (≥50% coverage):")
    for min_id, max_id in [(80, 85), (85, 90), (90, 95), (95, 100), (100, 101)]:
        tier = [r for r in results if min_id <= r['pident'] < max_id]
        if tier:
            n_syn = sum(1 for r in tier if r['any'] >= 0.5)
            label = f"{min_id}-{max_id}%" if max_id < 101 else "100%"
            print(f"  {label}: {n_syn}/{len(tier)} ({100*n_syn/len(tier):.1f}%)")

    # By length
    print("\nBy length (≥50% coverage):")
    for min_len, max_len in [(50, 100), (100, 200), (200, 500), (500, 1000)]:
        tier = [r for r in results if min_len <= r['length'] < max_len]
        if tier:
            n_syn = sum(1 for r in tier if r['any'] >= 0.5)
            print(f"  {min_len}-{max_len}bp: {n_syn}/{len(tier)} ({100*n_syn/len(tier):.1f}%)")


if __name__ == '__main__':
    main()
