#!/usr/bin/env python3
"""
Filter TE hits for ancient functional candidates using both conservation and synteny.

Criteria:
- Syntenic in 2+ species (any_coverage >= 0.5 for 2+ species)
- Conserved (phyloP > 1)
- Any identity (no filter)
"""

import sys
import os
from collections import defaultdict

def load_synteny_data(filepath):
    """Load synteny data indexed by (transcript, te_id)."""
    data = {}
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))

            # Parse name to get key
            name_parts = row['name'].split('|')
            if len(name_parts) >= 3:
                transcript = name_parts[1]
                te_id = name_parts[2]
                key = (transcript, te_id, row['chrom'], row['start'], row['end'])

                # Convert numeric fields
                row['pident'] = float(row['pident'])
                row['length'] = int(row['length'])
                row['sim_coverage'] = float(row['sim_coverage'])
                row['yak_coverage'] = float(row['yak_coverage'])
                row['ere_coverage'] = float(row['ere_coverage'])
                row['sec_coverage'] = float(row['sec_coverage'])
                row['any_coverage'] = float(row['any_coverage'])

                # Count species with synteny
                species_count = sum([
                    row['sim_coverage'] >= 0.5,
                    row['yak_coverage'] >= 0.5,
                    row['ere_coverage'] >= 0.5,
                    row['sec_coverage'] >= 0.5
                ])
                row['syntenic_species'] = species_count

                data[key] = row
    return data

def load_conservation_data(filepath):
    """Load conservation scores indexed by (transcript, te_id)."""
    data = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            phyloP = float(parts[5])  # mean conservation

            # Parse name
            name_parts = name.split('|')
            if len(name_parts) >= 3:
                transcript = name_parts[1]
                te_id = name_parts[2]
                # Use a simpler key for matching
                simple_key = (transcript, te_id)

                # Store with full name for later matching
                if simple_key not in data:
                    data[simple_key] = []
                data[simple_key].append({
                    'name': name,
                    'phyloP': phyloP,
                    'pident': float(name_parts[3]) if len(name_parts) > 3 else 0,
                    'length': int(name_parts[4]) if len(name_parts) > 4 else 0
                })
    return data

def load_blast_sequences(filepath, target_hits):
    """Load sequence data for specific hits from BLAST results."""
    # Build lookup set
    target_set = set()
    for hit in target_hits:
        target_set.add((hit['transcript'], hit['te_id']))

    sequences = {}
    with open(filepath) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            transcript = parts[0]
            te_id = parts[3]

            if (transcript, te_id) in target_set:
                key = (transcript, te_id)
                if key not in sequences:
                    sequences[key] = {
                        'qseqid': transcript,
                        'sseqid': te_id,
                        'qstart': int(parts[1]),
                        'qend': int(parts[2]),
                        'pident': float(parts[4]),
                        'length': int(parts[5]),
                        'evalue': parts[6]
                    }
    return sequences

def main():
    results_dir = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis'

    print("Loading synteny data...", file=sys.stderr)
    synteny = load_synteny_data(os.path.join(results_dir, 'te_hits_all_synteny_sampled.tsv'))
    print(f"  Loaded {len(synteny):,} hits with synteny data", file=sys.stderr)

    print("Loading conservation data...", file=sys.stderr)
    conservation = load_conservation_data(os.path.join(results_dir, 'te_hits_all_conservation.tab'))
    print(f"  Loaded conservation for {len(conservation):,} transcript-TE pairs", file=sys.stderr)

    # Match synteny with conservation
    print("\nMatching synteny with conservation...", file=sys.stderr)
    matched_hits = []

    for key, syn_row in synteny.items():
        transcript, te_id, chrom, start, end = key
        simple_key = (transcript, te_id)

        if simple_key in conservation:
            # Find best matching conservation entry (by pident/length)
            cons_entries = conservation[simple_key]
            best_cons = None
            best_match = 0

            for cons in cons_entries:
                # Match by pident and length similarity
                pident_match = 1 - abs(cons['pident'] - syn_row['pident']) / 100
                len_match = 1 - abs(cons['length'] - syn_row['length']) / max(cons['length'], syn_row['length'], 1)
                match_score = pident_match + len_match

                if match_score > best_match:
                    best_match = match_score
                    best_cons = cons

            if best_cons and best_match > 1.5:  # Reasonable match
                hit = {
                    'transcript': transcript,
                    'te_id': te_id,
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'pident': syn_row['pident'],
                    'length': syn_row['length'],
                    'category': syn_row['category'],
                    'phyloP': best_cons['phyloP'],
                    'sim_cov': syn_row['sim_coverage'],
                    'yak_cov': syn_row['yak_coverage'],
                    'ere_cov': syn_row['ere_coverage'],
                    'any_cov': syn_row['any_coverage'],
                    'syntenic_species': syn_row['syntenic_species']
                }
                matched_hits.append(hit)

    print(f"  Matched {len(matched_hits):,} hits", file=sys.stderr)

    # Apply filters for ancient functional candidates
    print("\nFiltering for ancient functional candidates...", file=sys.stderr)
    print("  Criteria: syntenic in ≥2 species AND phyloP > 1", file=sys.stderr)

    ancient_candidates = [
        h for h in matched_hits
        if h['syntenic_species'] >= 2 and h['phyloP'] > 1
    ]

    print(f"  Found {len(ancient_candidates):,} ancient functional candidates", file=sys.stderr)

    # Summary statistics
    print("\n" + "="*70)
    print("ANCIENT FUNCTIONAL TE CANDIDATE ANALYSIS")
    print("="*70)

    print(f"\nTotal matched hits: {len(matched_hits):,}")
    print(f"Ancient candidates (≥2 species, phyloP>1): {len(ancient_candidates):,}")
    print(f"  ({100*len(ancient_candidates)/len(matched_hits):.1f}% of matched hits)")

    # Distribution by filter
    print("\n" + "-"*50)
    print("Filter breakdown:")
    print("-"*50)

    syn_only = [h for h in matched_hits if h['syntenic_species'] >= 2]
    cons_only = [h for h in matched_hits if h['phyloP'] > 1]

    print(f"Syntenic in ≥2 species only: {len(syn_only):,}")
    print(f"Conserved (phyloP>1) only: {len(cons_only):,}")
    print(f"Both criteria: {len(ancient_candidates):,}")

    # Identity distribution
    print("\n" + "-"*50)
    print("Ancient candidates by identity:")
    print("-"*50)

    id_bins = [(60, 70), (70, 80), (80, 90), (90, 100)]
    for lo, hi in id_bins:
        subset = [h for h in ancient_candidates if lo <= h['pident'] < hi]
        pct = 100 * len(subset) / len(ancient_candidates) if ancient_candidates else 0
        print(f"  {lo}-{hi}%: {len(subset):,} ({pct:.1f}%)")

    # Length distribution
    print("\n" + "-"*50)
    print("Ancient candidates by length:")
    print("-"*50)

    len_bins = [(0, 30), (30, 50), (50, 100), (100, 200), (200, 500)]
    for lo, hi in len_bins:
        subset = [h for h in ancient_candidates if lo <= h['length'] < hi]
        pct = 100 * len(subset) / len(ancient_candidates) if ancient_candidates else 0
        print(f"  {lo}-{hi}bp: {len(subset):,} ({pct:.1f}%)")

    # TE families
    print("\n" + "-"*50)
    print("Top 15 TE instances in ancient candidates:")
    print("-"*50)

    te_counts = defaultdict(int)
    for h in ancient_candidates:
        te_counts[h['te_id']] += 1

    for te_id, count in sorted(te_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {te_id}: {count:,}")

    # Unique genes
    genes = set(h['transcript'] for h in ancient_candidates)
    print(f"\nUnique transcripts with ancient TE candidates: {len(genes):,}")

    # Top hits
    print("\n" + "="*70)
    print("TOP 30 ANCIENT FUNCTIONAL TE CANDIDATES")
    print("(Ranked by phyloP conservation)")
    print("="*70)

    # Sort by conservation
    top_hits = sorted(ancient_candidates, key=lambda x: -x['phyloP'])[:30]

    print(f"\n{'Rank':<5} {'Transcript':<15} {'TE':<15} {'phyloP':>8} {'Syn#':>5} {'ID%':>6} {'Len':>5} {'Chrom':<8}")
    print("-" * 75)

    for i, h in enumerate(top_hits, 1):
        print(f"{i:<5} {h['transcript']:<15} {h['te_id']:<15} {h['phyloP']:>8.2f} {h['syntenic_species']:>5} {h['pident']:>5.1f}% {h['length']:>5} {h['chrom']:<8}")

    # Save results
    output_file = os.path.join(results_dir, 'ancient_te_candidates.tsv')
    print(f"\nSaving {len(ancient_candidates):,} candidates to {output_file}...", file=sys.stderr)

    with open(output_file, 'w') as out:
        headers = ['transcript', 'te_id', 'chrom', 'start', 'end', 'pident', 'length',
                   'phyloP', 'syntenic_species', 'sim_cov', 'yak_cov', 'ere_cov', 'category']
        out.write('\t'.join(headers) + '\n')
        for h in sorted(ancient_candidates, key=lambda x: -x['phyloP']):
            out.write('\t'.join(str(h[k]) for k in headers) + '\n')

    # Save top 100 for detailed analysis
    top100_file = os.path.join(results_dir, 'ancient_te_candidates_top100.tsv')
    with open(top100_file, 'w') as out:
        headers = ['rank', 'transcript', 'te_id', 'chrom', 'start', 'end', 'pident', 'length',
                   'phyloP', 'syntenic_species', 'sim_cov', 'yak_cov', 'ere_cov', 'category']
        out.write('\t'.join(headers) + '\n')
        for i, h in enumerate(sorted(ancient_candidates, key=lambda x: -x['phyloP'])[:100], 1):
            out.write(f"{i}\t" + '\t'.join(str(h[k]) for k in headers[1:]) + '\n')

    print(f"Saved top 100 to {top100_file}", file=sys.stderr)

if __name__ == '__main__':
    main()
