#!/usr/bin/env python3
"""
Analyze synteny of TE hits using MAF alignments.

For each TE hit, check if the sequence is present (aligned) in D. simulans and D. yakuba.
"""

import gzip
import sys
import os
from collections import defaultdict
import bisect

# Species of interest (besides dm6 reference)
COMPARISON_SPECIES = ['droSim1', 'droYak2', 'droYak3', 'droEre2', 'droAna3', 'droSec1']

def parse_maf_blocks(maf_file):
    """Parse MAF file and yield alignment blocks with their genomic coordinates."""
    if maf_file.endswith('.gz'):
        fh = gzip.open(maf_file, 'rt')
    else:
        fh = open(maf_file)

    current_block = None
    ref_info = None

    for line in fh:
        if line.startswith('a '):
            # New alignment block
            if current_block and ref_info:
                yield ref_info, current_block
            current_block = {}
            ref_info = None

        elif line.startswith('s '):
            # Sequence line
            parts = line.strip().split()
            species_chrom = parts[1]
            start = int(parts[2])
            length = int(parts[3])
            strand = parts[4]
            src_len = int(parts[5])
            seq = parts[6] if len(parts) > 6 else ''

            # Parse species name
            if '.' in species_chrom:
                species = species_chrom.split('.')[0]
            else:
                species = species_chrom

            if species == 'dm6':
                # This is the reference
                ref_info = {
                    'chrom': species_chrom.replace('dm6.', ''),
                    'start': start,
                    'end': start + length,
                    'strand': strand,
                    'seq': seq
                }
                current_block['dm6'] = {
                    'start': start,
                    'end': start + length,
                    'seq': seq,
                    'aligned_length': length
                }
            elif species in COMPARISON_SPECIES or species.startswith('dro'):
                # Store comparison species
                current_block[species] = {
                    'start': start,
                    'end': start + length,
                    'seq': seq,
                    'aligned_length': length
                }

    # Don't forget the last block
    if current_block and ref_info:
        yield ref_info, current_block

    fh.close()


def build_maf_index(maf_file):
    """Build an index of MAF blocks by genomic position."""
    print(f"  Indexing {maf_file}...", file=sys.stderr)

    # Store (start, end, block_data) tuples sorted by start
    blocks = []
    block_count = 0

    for ref_info, block in parse_maf_blocks(maf_file):
        blocks.append((ref_info['start'], ref_info['end'], block))
        block_count += 1
        if block_count % 100000 == 0:
            print(f"    {block_count} blocks...", file=sys.stderr)

    # Sort by start position
    blocks.sort(key=lambda x: x[0])
    print(f"  Indexed {block_count} blocks", file=sys.stderr)

    return blocks


def find_overlapping_blocks(blocks, query_start, query_end):
    """Find all MAF blocks overlapping a query region using binary search."""
    # Find the first block that could overlap (start <= query_end)
    # Binary search for blocks where start <= query_end
    overlapping = []

    # Find insertion point for query_start
    starts = [b[0] for b in blocks]
    idx = bisect.bisect_right(starts, query_end)

    # Check blocks from idx backwards that might overlap
    for i in range(idx - 1, -1, -1):
        block_start, block_end, block_data = blocks[i]
        if block_end <= query_start:
            # No more overlaps possible (blocks are sorted by start)
            # But we need to keep going since a block starting earlier could extend further
            continue
        if block_start > query_end:
            continue
        if block_start < query_end and block_end > query_start:
            overlapping.append((block_start, block_end, block_data))

    # Also check forward from the insertion point
    for i in range(idx, len(blocks)):
        block_start, block_end, block_data = blocks[i]
        if block_start >= query_end:
            break
        if block_start < query_end and block_end > query_start:
            overlapping.append((block_start, block_end, block_data))

    return overlapping


def calculate_synteny_coverage(blocks, query_start, query_end, species):
    """Calculate what fraction of the query region has alignment to the given species."""
    if not blocks:
        return 0.0, 0

    query_len = query_end - query_start
    covered_bases = 0

    # Track covered positions to avoid double-counting overlapping blocks
    covered = set()

    for block_start, block_end, block_data in blocks:
        if species not in block_data:
            continue

        # Calculate overlap
        overlap_start = max(query_start, block_start)
        overlap_end = min(query_end, block_end)

        if overlap_end > overlap_start:
            for pos in range(overlap_start, overlap_end):
                covered.add(pos)

    coverage = len(covered) / query_len if query_len > 0 else 0
    return coverage, len(covered)


def analyze_te_hits(bed_file, maf_indices, output_file, sample_rate=1):
    """Analyze synteny for TE hits."""
    print(f"\nAnalyzing {bed_file}...", file=sys.stderr)

    results = []
    processed = 0
    skipped_no_maf = 0

    with open(bed_file) as f:
        for i, line in enumerate(f):
            if sample_rate > 1 and i % sample_rate != 0:
                continue

            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3]

            # Parse metadata from name
            name_parts = name.split('|')
            category = name_parts[0].split('_')[0] if name_parts else 'unknown'
            pident = float(name_parts[3]) if len(name_parts) > 3 else 0
            length = int(name_parts[4]) if len(name_parts) > 4 else 0

            # Get MAF index for this chromosome
            if chrom not in maf_indices:
                skipped_no_maf += 1
                continue

            blocks = maf_indices[chrom]

            # Find overlapping MAF blocks
            overlapping = find_overlapping_blocks(blocks, start, end)

            # Calculate coverage for each comparison species
            sim_cov, sim_bases = calculate_synteny_coverage(overlapping, start, end, 'droSim1')
            yak_cov, yak_bases = calculate_synteny_coverage(overlapping, start, end, 'droYak2')
            # Try droYak3 if droYak2 not found
            if yak_cov == 0:
                yak_cov, yak_bases = calculate_synteny_coverage(overlapping, start, end, 'droYak3')
            ere_cov, ere_bases = calculate_synteny_coverage(overlapping, start, end, 'droEre2')
            sec_cov, sec_bases = calculate_synteny_coverage(overlapping, start, end, 'droSec1')

            # Classify synteny status
            any_species = max(sim_cov, yak_cov, ere_cov, sec_cov)
            both_close = min(sim_cov, sec_cov)  # Both closest relatives

            results.append({
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
                'any_coverage': any_species,
                'both_close_coverage': both_close,
                'n_blocks': len(overlapping)
            })

            processed += 1
            if processed % 10000 == 0:
                print(f"  Processed {processed} hits...", file=sys.stderr)

    print(f"  Done: {processed} analyzed, {skipped_no_maf} skipped (no MAF)", file=sys.stderr)

    # Write results
    with open(output_file, 'w') as out:
        headers = ['chrom', 'start', 'end', 'name', 'category', 'pident', 'length',
                   'sim_coverage', 'yak_coverage', 'ere_coverage', 'sec_coverage',
                   'any_coverage', 'both_close_coverage', 'n_blocks']
        out.write('\t'.join(headers) + '\n')
        for r in results:
            out.write('\t'.join(str(r[h]) for h in headers) + '\n')

    return results


def main():
    # Paths
    maf_dir = '/Users/jacobboysen/git_repos/repeat_finder/data/references/maf'
    results_dir = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis'

    # Build MAF indices for each chromosome
    print("Building MAF indices...", file=sys.stderr)
    maf_indices = {}

    for chrom in ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX']:
        maf_file = os.path.join(maf_dir, f'{chrom}.maf.gz')
        if os.path.exists(maf_file):
            maf_indices[chrom] = build_maf_index(maf_file)

    # Analyze high-quality hits first (smaller set)
    print("\n" + "="*70, file=sys.stderr)
    print("Analyzing HIGH-QUALITY hits (≥80% id, ≥50bp)", file=sys.stderr)
    print("="*70, file=sys.stderr)

    # Filter HQ hits from the all hits file
    hq_bed = os.path.join(results_dir, 'te_hits_hq_for_synteny.bed')
    with open(os.path.join(results_dir, 'te_hits_all_genomic.bed')) as fin:
        with open(hq_bed, 'w') as fout:
            for line in fin:
                parts = line.strip().split('\t')
                name = parts[3]
                name_parts = name.split('|')
                if len(name_parts) >= 5:
                    pident = float(name_parts[3])
                    length = int(name_parts[4])
                    if pident >= 80 and length >= 50:
                        fout.write(line)

    hq_results = analyze_te_hits(
        hq_bed,
        maf_indices,
        os.path.join(results_dir, 'te_hits_hq_synteny.tsv')
    )

    # Print HQ summary
    print("\n" + "="*70)
    print("HIGH-QUALITY HITS SYNTENY SUMMARY")
    print("="*70)

    n_total = len(hq_results)
    if n_total > 0:
        n_sim = sum(1 for r in hq_results if r['sim_coverage'] >= 0.5)
        n_yak = sum(1 for r in hq_results if r['yak_coverage'] >= 0.5)
        n_any = sum(1 for r in hq_results if r['any_coverage'] >= 0.5)
        n_both = sum(1 for r in hq_results if r['both_close_coverage'] >= 0.5)

        print(f"\nTotal HQ hits: {n_total:,}")
        print(f"Syntenic in D. simulans (≥50% coverage): {n_sim:,} ({100*n_sim/n_total:.1f}%)")
        print(f"Syntenic in D. yakuba (≥50% coverage): {n_yak:,} ({100*n_yak/n_total:.1f}%)")
        print(f"Syntenic in any species (≥50% coverage): {n_any:,} ({100*n_any/n_total:.1f}%)")
        print(f"Syntenic in BOTH sim+sec (≥50% coverage): {n_both:,} ({100*n_both/n_total:.1f}%)")

    # Now analyze ALL hits (will take longer)
    print("\n" + "="*70, file=sys.stderr)
    print("Analyzing ALL hits (sampling every 10th hit for speed)", file=sys.stderr)
    print("="*70, file=sys.stderr)

    all_results = analyze_te_hits(
        os.path.join(results_dir, 'te_hits_all_genomic.bed'),
        maf_indices,
        os.path.join(results_dir, 'te_hits_all_synteny_sampled.tsv'),
        sample_rate=10  # Sample every 10th hit for speed
    )

    # Print ALL hits summary
    print("\n" + "="*70)
    print("ALL HITS SYNTENY SUMMARY (10% sample)")
    print("="*70)

    n_total = len(all_results)
    if n_total > 0:
        n_sim = sum(1 for r in all_results if r['sim_coverage'] >= 0.5)
        n_yak = sum(1 for r in all_results if r['yak_coverage'] >= 0.5)
        n_any = sum(1 for r in all_results if r['any_coverage'] >= 0.5)

        print(f"\nTotal sampled hits: {n_total:,} (est. {n_total*10:,} total)")
        print(f"Syntenic in D. simulans (≥50% coverage): {n_sim:,} ({100*n_sim/n_total:.1f}%)")
        print(f"Syntenic in D. yakuba (≥50% coverage): {n_yak:,} ({100*n_yak/n_total:.1f}%)")
        print(f"Syntenic in any species (≥50% coverage): {n_any:,} ({100*n_any/n_total:.1f}%)")

        # Break down by quality tier
        print("\n" + "-"*50)
        print("Synteny by quality tier (any species ≥50%):")
        print("-"*50)

        for min_id, max_id in [(60, 70), (70, 80), (80, 90), (90, 100)]:
            tier = [r for r in all_results if min_id <= r['pident'] < max_id]
            if tier:
                n_syn = sum(1 for r in tier if r['any_coverage'] >= 0.5)
                print(f"  {min_id}-{max_id}% identity: {n_syn}/{len(tier)} ({100*n_syn/len(tier):.1f}%) syntenic")


if __name__ == '__main__':
    main()
