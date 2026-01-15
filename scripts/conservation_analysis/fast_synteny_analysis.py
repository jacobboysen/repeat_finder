#!/usr/bin/env python3
"""
Fast synteny analysis using sorted merge approach.

Strategy:
1. Extract species coverage from MAF files into simple BED format
2. Sort both TE hits and species coverage by position
3. Use merge-sort style overlap detection for O(n+m) complexity
"""

import gzip
import sys
import os
from collections import defaultdict

def extract_species_coverage_from_maf(maf_file, chrom):
    """Extract coverage intervals for each species from MAF file."""
    print(f"  Extracting coverage from {maf_file}...", file=sys.stderr)

    # Species coverage: list of (start, end) intervals
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
            # Save previous block
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

    # Don't forget last block
    if current_ref_start is not None:
        for sp in current_species:
            species_intervals[sp].append((current_ref_start, current_ref_end))

    fh.close()

    # Merge overlapping intervals and sort
    for sp in species_intervals:
        species_intervals[sp] = merge_intervals(sorted(species_intervals[sp]))
        print(f"    {sp}: {len(species_intervals[sp])} merged intervals", file=sys.stderr)

    return dict(species_intervals)


def merge_intervals(intervals):
    """Merge overlapping intervals."""
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
    """Calculate what fraction of query is covered by intervals (binary search + scan)."""
    if not intervals:
        return 0.0

    query_len = query_end - query_start
    if query_len <= 0:
        return 0.0

    # Binary search for first interval that could overlap
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= query_start:
            lo = mid + 1
        else:
            hi = mid

    # Count covered bases
    covered = 0
    for i in range(lo, len(intervals)):
        int_start, int_end = intervals[i]
        if int_start >= query_end:
            break

        # Calculate overlap
        overlap_start = max(query_start, int_start)
        overlap_end = min(query_end, int_end)
        if overlap_end > overlap_start:
            covered += overlap_end - overlap_start

    return covered / query_len


def analyze_chromosome(chrom, maf_file, te_hits, output_results):
    """Analyze synteny for all TE hits on one chromosome."""
    # Extract species coverage from MAF
    coverage_intervals = extract_species_coverage_from_maf(maf_file, chrom)

    print(f"  Processing {len(te_hits)} hits on {chrom}...", file=sys.stderr)

    for hit in te_hits:
        start, end, name, pident, length, category = hit

        # Calculate coverage for each species
        sim_cov = calculate_coverage(start, end,
                    coverage_intervals.get('droSim1', []))
        sec_cov = calculate_coverage(start, end,
                    coverage_intervals.get('droSec1', []))
        yak_cov = max(
            calculate_coverage(start, end, coverage_intervals.get('droYak2', [])),
            calculate_coverage(start, end, coverage_intervals.get('droYak3', []))
        )
        ere_cov = calculate_coverage(start, end,
                    coverage_intervals.get('droEre2', []))

        any_cov = max(sim_cov, sec_cov, yak_cov, ere_cov)
        both_close = min(sim_cov, sec_cov)

        output_results.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'name': name,
            'category': category,
            'pident': pident,
            'length': length,
            'sim_coverage': sim_cov,
            'yak_coverage': yak_cov,
            'ere_coverage': ere_cov,
            'sec_coverage': sec_cov,
            'any_coverage': any_cov,
            'both_close_coverage': both_close
        })


def load_te_hits(bed_file, sample_rate=1):
    """Load TE hits from BED file, grouped by chromosome."""
    hits_by_chrom = defaultdict(list)

    with open(bed_file) as f:
        for i, line in enumerate(f):
            if sample_rate > 1 and i % sample_rate != 0:
                continue

            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]

            # Parse metadata
            name_parts = name.split('|')
            category = name_parts[0].split('_')[0] if name_parts else 'unknown'
            pident = float(name_parts[3]) if len(name_parts) > 3 else 0
            length = int(name_parts[4]) if len(name_parts) > 4 else 0

            hits_by_chrom[chrom].append((start, end, name, pident, length, category))

    # Sort by position
    for chrom in hits_by_chrom:
        hits_by_chrom[chrom].sort()

    return dict(hits_by_chrom)


def main():
    maf_dir = '/Users/jacobboysen/git_repos/repeat_finder/data/references/maf'
    results_dir = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis'

    chromosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']

    # Analyze ALL hits with sampling
    print("="*70, file=sys.stderr)
    print("FAST SYNTENY ANALYSIS - ALL HITS (10% sample)", file=sys.stderr)
    print("="*70, file=sys.stderr)

    print("\nLoading TE hits (sampling every 10th)...", file=sys.stderr)
    all_hits = load_te_hits(
        os.path.join(results_dir, 'te_hits_all_genomic.bed'),
        sample_rate=10
    )

    total_hits = sum(len(hits) for hits in all_hits.values())
    print(f"Loaded {total_hits:,} sampled hits", file=sys.stderr)

    results = []
    for chrom in chromosomes:
        maf_file = os.path.join(maf_dir, f'{chrom}.maf.gz')
        if chrom in all_hits and os.path.exists(maf_file):
            analyze_chromosome(chrom, maf_file, all_hits[chrom], results)

    # Write results
    output_file = os.path.join(results_dir, 'te_hits_all_synteny_sampled.tsv')
    print(f"\nWriting results to {output_file}...", file=sys.stderr)

    with open(output_file, 'w') as out:
        headers = ['chrom', 'start', 'end', 'name', 'category', 'pident', 'length',
                   'sim_coverage', 'yak_coverage', 'ere_coverage', 'sec_coverage',
                   'any_coverage', 'both_close_coverage']
        out.write('\t'.join(headers) + '\n')
        for r in results:
            out.write('\t'.join(str(r[h]) for h in headers) + '\n')

    print(f"Wrote {len(results):,} results", file=sys.stderr)

    # Print summary
    print("\n" + "="*70)
    print("ALL HITS SYNTENY SUMMARY (10% sample)")
    print("="*70)

    n = len(results)
    if n > 0:
        n_sim = sum(1 for r in results if r['sim_coverage'] >= 0.5)
        n_yak = sum(1 for r in results if r['yak_coverage'] >= 0.5)
        n_any = sum(1 for r in results if r['any_coverage'] >= 0.5)

        print(f"\nTotal sampled hits: {n:,}")
        print(f"Syntenic in D. simulans (≥50%): {n_sim:,} ({100*n_sim/n:.1f}%)")
        print(f"Syntenic in D. yakuba (≥50%): {n_yak:,} ({100*n_yak/n:.1f}%)")
        print(f"Syntenic in any species (≥50%): {n_any:,} ({100*n_any/n:.1f}%)")

        # By quality tier
        print("\nSynteny by identity tier (any species ≥50%):")
        for min_id, max_id in [(60, 70), (70, 80), (80, 90), (90, 100)]:
            tier = [r for r in results if min_id <= r['pident'] < max_id]
            if tier:
                n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
                print(f"  {min_id}-{max_id}%: {n_syn}/{len(tier)} ({100*n_syn/len(tier):.1f}%)")


if __name__ == '__main__':
    main()
