#!/usr/bin/env python3
"""
Filter TE hits for ancient functional candidates with improved quality controls.

Version 2 adds critical filters missing from v1:
1. CDS overlap exclusion - removes hits overlapping coding sequences in any isoform
2. Hit deduplication - collapses overlapping hits from different TE families into loci
3. Sequence complexity filter - removes low-complexity/AT-rich spurious matches

These filters address systematic false positives identified in the v1 analysis where:
- Top candidates were actually CDS regions with coding-level conservation
- Single AT-rich regions generated hits to multiple unrelated TE families
- HeT-A (telomeric) and other family matches were spurious low-complexity matches

Usage:
    python scripts/filter_ancient_te_candidates_v2.py [--output-dir DIR]
"""

import sys
import os
import math
import argparse
from collections import defaultdict
from pathlib import Path

# Add scripts directory to path for utils
sys.path.insert(0, str(Path(__file__).parent))

from utils.paths import get_project_root, get_results_dir, get_references_dir


def calculate_shannon_entropy(sequence: str) -> float:
    """
    Calculate Shannon entropy of a DNA sequence.

    Higher entropy = more complex sequence.
    Lower entropy = more repetitive/biased sequence.

    Max entropy for DNA = 2.0 (equal base frequencies)
    Poly-A has entropy = 0.0

    Args:
        sequence: DNA sequence string

    Returns:
        Shannon entropy (0.0 to 2.0)
    """
    if not sequence:
        return 0.0

    sequence = sequence.upper()
    length = len(sequence)

    # Count nucleotide frequencies
    counts = defaultdict(int)
    for base in sequence:
        if base in 'ATGC':
            counts[base] += 1

    # Calculate entropy
    entropy = 0.0
    for base in 'ATGC':
        if counts[base] > 0:
            freq = counts[base] / length
            entropy -= freq * math.log2(freq)

    return entropy


def calculate_at_content(sequence: str) -> float:
    """Calculate AT content of a sequence."""
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    at = sum(1 for b in sequence if b in 'AT')
    total = sum(1 for b in sequence if b in 'ATGC')
    return at / total if total > 0 else 0.0


def count_homopolymer_runs(sequence: str, min_length: int = 4) -> int:
    """Count homopolymer runs of at least min_length bases."""
    if not sequence:
        return 0

    sequence = sequence.upper()
    runs = 0
    i = 0
    while i < len(sequence):
        if sequence[i] in 'ATGC':
            run_len = 1
            while i + run_len < len(sequence) and sequence[i + run_len] == sequence[i]:
                run_len += 1
            if run_len >= min_length:
                runs += 1
            i += run_len
        else:
            i += 1
    return runs


def load_cds_regions(gff_path: str) -> dict:
    """
    Load CDS regions from GFF file.

    Returns:
        Dictionary mapping chromosome to list of (start, end) tuples
    """
    cds_regions = defaultdict(list)

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'CDS':
                continue

            chrom = 'chr' + parts[0]  # Add chr prefix for consistency
            start = int(parts[3])
            end = int(parts[4])
            cds_regions[chrom].append((start, end))

    # Sort each chromosome's CDS regions for efficient lookup
    for chrom in cds_regions:
        cds_regions[chrom].sort()

    return cds_regions


def overlaps_cds(chrom: str, start: int, end: int, cds_regions: dict) -> bool:
    """Check if a region overlaps any CDS."""
    if chrom not in cds_regions:
        return False

    for cds_start, cds_end in cds_regions[chrom]:
        if max(start, cds_start) < min(end, cds_end):
            return True
        # Early exit if we've passed the region
        if cds_start > end:
            break

    return False


def load_te_families(te_fasta_path: str) -> dict:
    """Load TE ID to family mapping from FASTA headers."""
    families = {}
    with open(te_fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split()
                te_id = parts[0][1:]  # Remove >
                for part in parts:
                    if part.startswith('name='):
                        name = part.split('=')[1].split('{')[0]
                        families[te_id] = name
                        break
    return families


def collapse_overlapping_hits(hits: list, max_gap: int = 10) -> list:
    """
    Collapse overlapping or adjacent hits into single loci.

    When multiple TE families hit the same region, this likely indicates
    spurious matches to low-complexity sequence rather than independent
    TE insertions.

    Args:
        hits: List of hit dictionaries with chrom, start, end
        max_gap: Maximum gap between hits to merge (default 10bp)

    Returns:
        List of collapsed loci with metadata about constituent hits
    """
    if not hits:
        return []

    # Group hits by transcript and chromosome
    grouped = defaultdict(list)
    for hit in hits:
        key = (hit['transcript'], hit['chrom'])
        grouped[key].append(hit)

    collapsed = []

    for (transcript, chrom), group_hits in grouped.items():
        # Sort by start position
        sorted_hits = sorted(group_hits, key=lambda x: x['start'])

        # Merge overlapping/adjacent hits
        loci = []
        current_locus = None

        for hit in sorted_hits:
            if current_locus is None:
                current_locus = {
                    'transcript': transcript,
                    'chrom': chrom,
                    'start': hit['start'],
                    'end': hit['end'],
                    'constituent_hits': [hit],
                    'te_ids': {hit['te_id']},
                    'te_families': {hit.get('te_family', 'unknown')},
                    'best_phyloP': hit['phyloP'],
                    'best_pident': hit['pident'],
                    'best_length': hit['length']
                }
            elif hit['start'] <= current_locus['end'] + max_gap:
                # Merge into current locus
                current_locus['end'] = max(current_locus['end'], hit['end'])
                current_locus['constituent_hits'].append(hit)
                current_locus['te_ids'].add(hit['te_id'])
                current_locus['te_families'].add(hit.get('te_family', 'unknown'))
                if hit['phyloP'] > current_locus['best_phyloP']:
                    current_locus['best_phyloP'] = hit['phyloP']
                    current_locus['best_pident'] = hit['pident']
                    current_locus['best_length'] = hit['length']
            else:
                # Start new locus
                loci.append(current_locus)
                current_locus = {
                    'transcript': transcript,
                    'chrom': chrom,
                    'start': hit['start'],
                    'end': hit['end'],
                    'constituent_hits': [hit],
                    'te_ids': {hit['te_id']},
                    'te_families': {hit.get('te_family', 'unknown')},
                    'best_phyloP': hit['phyloP'],
                    'best_pident': hit['pident'],
                    'best_length': hit['length']
                }

        if current_locus:
            loci.append(current_locus)

        collapsed.extend(loci)

    return collapsed


def load_blast_sequences(blast_tsv_path: str, target_keys: set) -> dict:
    """Load aligned sequences from BLAST results for complexity analysis."""
    sequences = {}

    with open(blast_tsv_path) as f:
        header = next(f).strip().split('\t')
        col_idx = {col: i for i, col in enumerate(header)}

        # Check for qseq column
        if 'qseq' not in col_idx:
            print(f"  Warning: No qseq column in {blast_tsv_path}", file=sys.stderr)
            return sequences

        qseq_idx = col_idx['qseq']
        qseqid_idx = col_idx.get('qseqid', 0)
        sseqid_idx = col_idx.get('sseqid', 1)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(qseq_idx, qseqid_idx, sseqid_idx):
                continue

            key = (parts[qseqid_idx], parts[sseqid_idx])
            if key in target_keys:
                # Remove gaps from sequence
                qseq = parts[qseq_idx].replace('-', '')
                sequences[key] = qseq

    return sequences


def load_synteny_data(filepath: str) -> dict:
    """Load synteny data indexed by hit key."""
    data = {}
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))

            name_parts = row['name'].split('|')
            if len(name_parts) >= 3:
                transcript = name_parts[1]
                te_id = name_parts[2]
                key = (transcript, te_id, row['chrom'], row['start'], row['end'])

                row['pident'] = float(row['pident'])
                row['length'] = int(row['length'])
                row['sim_coverage'] = float(row['sim_coverage'])
                row['yak_coverage'] = float(row['yak_coverage'])
                row['ere_coverage'] = float(row['ere_coverage'])
                row['sec_coverage'] = float(row['sec_coverage'])

                species_count = sum([
                    row['sim_coverage'] >= 0.5,
                    row['yak_coverage'] >= 0.5,
                    row['ere_coverage'] >= 0.5,
                    row['sec_coverage'] >= 0.5
                ])
                row['syntenic_species'] = species_count

                data[key] = row
    return data


def load_conservation_data(filepath: str) -> dict:
    """Load conservation scores indexed by (transcript, te_id)."""
    data = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            phyloP = float(parts[5])

            name_parts = name.split('|')
            if len(name_parts) >= 3:
                transcript = name_parts[1]
                te_id = name_parts[2]
                simple_key = (transcript, te_id)

                if simple_key not in data:
                    data[simple_key] = []
                data[simple_key].append({
                    'name': name,
                    'phyloP': phyloP,
                    'pident': float(name_parts[3]) if len(name_parts) > 3 else 0,
                    'length': int(name_parts[4]) if len(name_parts) > 4 else 0
                })
    return data


def main():
    parser = argparse.ArgumentParser(description='Filter ancient TE candidates with improved QC')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: results/repeatmasker_analysis_v2)')
    parser.add_argument('--min-entropy', type=float, default=1.5,
                        help='Minimum Shannon entropy (default: 1.5, max=2.0)')
    parser.add_argument('--max-at-content', type=float, default=0.75,
                        help='Maximum AT content (default: 0.75)')
    parser.add_argument('--max-homopolymers', type=int, default=3,
                        help='Maximum number of 4+ bp homopolymer runs (default: 3)')
    args = parser.parse_args()

    # Set up paths
    project_root = get_project_root()
    results_dir = get_results_dir() / 'repeatmasker_analysis'

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = get_results_dir() / 'repeatmasker_analysis_v2'

    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70, file=sys.stderr)
    print("ANCIENT TE CANDIDATE FILTERING v2", file=sys.stderr)
    print("="*70, file=sys.stderr)
    print(f"\nFilters applied:", file=sys.stderr)
    print(f"  - CDS overlap exclusion", file=sys.stderr)
    print(f"  - Hit deduplication (collapse overlapping loci)", file=sys.stderr)
    print(f"  - Sequence complexity: entropy >= {args.min_entropy}", file=sys.stderr)
    print(f"  - AT content <= {args.max_at_content*100:.0f}%", file=sys.stderr)
    print(f"  - Homopolymer runs (4+bp) <= {args.max_homopolymers}", file=sys.stderr)

    # Load reference data
    print("\n[1/6] Loading CDS regions from GFF...", file=sys.stderr)
    gff_path = get_references_dir() / 'dmel-all-r6.66.gff'
    cds_regions = load_cds_regions(str(gff_path))
    total_cds = sum(len(v) for v in cds_regions.values())
    print(f"  Loaded {total_cds:,} CDS regions across {len(cds_regions)} chromosomes", file=sys.stderr)

    print("\n[2/6] Loading TE family annotations...", file=sys.stderr)
    te_fasta_path = get_references_dir() / 'dmel_te_flybase.fasta'
    te_families = load_te_families(str(te_fasta_path))
    print(f"  Loaded {len(te_families):,} TE family mappings", file=sys.stderr)

    print("\n[3/6] Loading synteny data...", file=sys.stderr)
    synteny = load_synteny_data(str(results_dir / 'te_hits_all_synteny_sampled.tsv'))
    print(f"  Loaded {len(synteny):,} hits with synteny data", file=sys.stderr)

    print("\n[4/6] Loading conservation data...", file=sys.stderr)
    conservation = load_conservation_data(str(results_dir / 'te_hits_all_conservation.tab'))
    print(f"  Loaded conservation for {len(conservation):,} transcript-TE pairs", file=sys.stderr)

    # Match synteny with conservation
    print("\n[5/6] Matching and filtering hits...", file=sys.stderr)
    matched_hits = []

    for key, syn_row in synteny.items():
        transcript, te_id, chrom, start, end = key
        simple_key = (transcript, te_id)

        if simple_key in conservation:
            cons_entries = conservation[simple_key]
            best_cons = None
            best_match = 0

            for cons in cons_entries:
                pident_match = 1 - abs(cons['pident'] - syn_row['pident']) / 100
                len_match = 1 - abs(cons['length'] - syn_row['length']) / max(cons['length'], syn_row['length'], 1)
                match_score = pident_match + len_match

                if match_score > best_match:
                    best_match = match_score
                    best_cons = cons

            if best_cons and best_match > 1.5:
                hit = {
                    'transcript': transcript,
                    'te_id': te_id,
                    'te_family': te_families.get(te_id, 'unknown'),
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
                    'syntenic_species': syn_row['syntenic_species']
                }
                matched_hits.append(hit)

    print(f"  Matched {len(matched_hits):,} hits", file=sys.stderr)

    # Apply base filters (synteny + conservation)
    base_filtered = [
        h for h in matched_hits
        if h['syntenic_species'] >= 2 and h['phyloP'] > 1
    ]
    print(f"  After synteny+conservation filter: {len(base_filtered):,}", file=sys.stderr)

    # Filter 1: Remove CDS overlaps
    pre_cds = len(base_filtered)
    cds_filtered = [
        h for h in base_filtered
        if not overlaps_cds(h['chrom'], h['start'], h['end'], cds_regions)
    ]
    cds_removed = pre_cds - len(cds_filtered)
    print(f"  After CDS exclusion: {len(cds_filtered):,} ({cds_removed:,} removed)", file=sys.stderr)

    # Load sequences for complexity analysis
    print("\n[6/6] Loading sequences for complexity analysis...", file=sys.stderr)
    target_keys = {(h['transcript'], h['te_id']) for h in cds_filtered}

    # Try to load from genome-wide BLAST results
    blast_path = get_results_dir() / 'genome_wide_all_3utrs.tsv'
    sequences = {}
    if blast_path.exists():
        sequences = load_blast_sequences(str(blast_path), target_keys)
        print(f"  Loaded {len(sequences):,} sequences for complexity analysis", file=sys.stderr)
    else:
        print(f"  Warning: BLAST file not found, skipping complexity filter", file=sys.stderr)

    # Filter 2: Sequence complexity
    if sequences:
        complexity_filtered = []
        complexity_stats = {'low_entropy': 0, 'high_at': 0, 'homopolymer': 0, 'passed': 0}

        for hit in cds_filtered:
            key = (hit['transcript'], hit['te_id'])
            if key in sequences:
                seq = sequences[key]
                entropy = calculate_shannon_entropy(seq)
                at_content = calculate_at_content(seq)
                homopolymers = count_homopolymer_runs(seq)

                hit['entropy'] = entropy
                hit['at_content'] = at_content
                hit['homopolymers'] = homopolymers

                # Apply filters
                if entropy < args.min_entropy:
                    complexity_stats['low_entropy'] += 1
                    continue
                if at_content > args.max_at_content:
                    complexity_stats['high_at'] += 1
                    continue
                if homopolymers > args.max_homopolymers:
                    complexity_stats['homopolymer'] += 1
                    continue

                complexity_stats['passed'] += 1
                complexity_filtered.append(hit)
            else:
                # No sequence data - keep hit but flag it
                hit['entropy'] = None
                hit['at_content'] = None
                hit['homopolymers'] = None
                complexity_filtered.append(hit)

        print(f"  Complexity filter results:", file=sys.stderr)
        print(f"    Low entropy (<{args.min_entropy}): {complexity_stats['low_entropy']:,} removed", file=sys.stderr)
        print(f"    High AT (>{args.max_at_content*100:.0f}%): {complexity_stats['high_at']:,} removed", file=sys.stderr)
        print(f"    Homopolymer-rich: {complexity_stats['homopolymer']:,} removed", file=sys.stderr)
        print(f"    Passed: {complexity_stats['passed']:,}", file=sys.stderr)
    else:
        complexity_filtered = cds_filtered

    # Filter 3: Collapse overlapping hits
    print("\n  Collapsing overlapping hits into loci...", file=sys.stderr)
    collapsed_loci = collapse_overlapping_hits(complexity_filtered)

    # Flag suspicious loci (multiple TE families)
    suspicious_loci = [l for l in collapsed_loci if len(l['te_families']) > 1]
    clean_loci = [l for l in collapsed_loci if len(l['te_families']) == 1]

    print(f"  Collapsed into {len(collapsed_loci):,} unique loci", file=sys.stderr)
    print(f"    Clean (single TE family): {len(clean_loci):,}", file=sys.stderr)
    print(f"    Suspicious (multiple families): {len(suspicious_loci):,}", file=sys.stderr)

    # Summary
    print("\n" + "="*70, file=sys.stderr)
    print("FILTERING SUMMARY", file=sys.stderr)
    print("="*70, file=sys.stderr)
    print(f"Starting hits (matched): {len(matched_hits):,}", file=sys.stderr)
    print(f"After synteny+conservation: {len(base_filtered):,}", file=sys.stderr)
    print(f"After CDS exclusion: {len(cds_filtered):,}", file=sys.stderr)
    print(f"After complexity filter: {len(complexity_filtered):,}", file=sys.stderr)
    print(f"Unique loci: {len(collapsed_loci):,}", file=sys.stderr)
    print(f"Clean loci (final): {len(clean_loci):,}", file=sys.stderr)

    # Save results
    print(f"\nSaving results to {output_dir}/...", file=sys.stderr)

    # Save clean loci (primary output)
    clean_output = output_dir / 'ancient_te_candidates_clean.tsv'
    with open(clean_output, 'w') as out:
        headers = ['rank', 'transcript', 'chrom', 'start', 'end', 'te_family', 'te_id',
                   'phyloP', 'pident', 'length', 'syntenic_species', 'entropy', 'at_content']
        out.write('\t'.join(headers) + '\n')

        for i, locus in enumerate(sorted(clean_loci, key=lambda x: -x['best_phyloP']), 1):
            # Get representative hit
            best_hit = locus['constituent_hits'][0]
            for h in locus['constituent_hits']:
                if h['phyloP'] == locus['best_phyloP']:
                    best_hit = h
                    break

            row = [
                str(i),
                locus['transcript'],
                locus['chrom'],
                str(locus['start']),
                str(locus['end']),
                list(locus['te_families'])[0],
                list(locus['te_ids'])[0],
                f"{locus['best_phyloP']:.4f}",
                f"{locus['best_pident']:.1f}",
                str(locus['best_length']),
                str(best_hit.get('syntenic_species', 'NA')),
                f"{best_hit.get('entropy', 'NA'):.3f}" if best_hit.get('entropy') else 'NA',
                f"{best_hit.get('at_content', 'NA'):.3f}" if best_hit.get('at_content') else 'NA'
            ]
            out.write('\t'.join(row) + '\n')

    print(f"  Saved {len(clean_loci):,} clean candidates to {clean_output.name}", file=sys.stderr)

    # Save suspicious loci (for review)
    suspicious_output = output_dir / 'ancient_te_candidates_suspicious.tsv'
    with open(suspicious_output, 'w') as out:
        headers = ['rank', 'transcript', 'chrom', 'start', 'end', 'n_te_families',
                   'te_families', 'n_hits', 'best_phyloP', 'issue']
        out.write('\t'.join(headers) + '\n')

        for i, locus in enumerate(sorted(suspicious_loci, key=lambda x: -x['best_phyloP']), 1):
            families_str = ','.join(sorted(locus['te_families']))
            row = [
                str(i),
                locus['transcript'],
                locus['chrom'],
                str(locus['start']),
                str(locus['end']),
                str(len(locus['te_families'])),
                families_str,
                str(len(locus['constituent_hits'])),
                f"{locus['best_phyloP']:.4f}",
                'multi_family_overlap'
            ]
            out.write('\t'.join(row) + '\n')

    print(f"  Saved {len(suspicious_loci):,} suspicious loci to {suspicious_output.name}", file=sys.stderr)

    # Save filtering report
    report_output = output_dir / 'FILTERING_REPORT_v2.md'
    with open(report_output, 'w') as out:
        out.write("# Ancient TE Candidate Filtering Report v2\n\n")
        out.write(f"**Generated**: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")

        out.write("## Issues Addressed\n\n")
        out.write("This filtering addresses systematic false positives identified in v1:\n\n")
        out.write("1. **CDS Overlap**: Top candidates were 3'UTR regions that overlap CDS in alternative\n")
        out.write("   transcript isoforms, explaining their coding-level conservation (phyloP ~5)\n\n")
        out.write("2. **Multi-family Hits**: Single low-complexity regions generating spurious matches\n")
        out.write("   to multiple unrelated TE families (e.g., nAChRalpha5 hitting 1360, antonia,\n")
        out.write("   Stalker, and HeT-A in the same 60bp window)\n\n")
        out.write("3. **AT-rich Sequences**: Poly-A regions bypassing DUST filter but generating\n")
        out.write("   spurious alignments to AT-rich TEs\n\n")

        out.write("## Filter Parameters\n\n")
        out.write(f"- Minimum Shannon entropy: {args.min_entropy} (max=2.0)\n")
        out.write(f"- Maximum AT content: {args.max_at_content*100:.0f}%\n")
        out.write(f"- Maximum homopolymer runs (4+bp): {args.max_homopolymers}\n\n")

        out.write("## Results Summary\n\n")
        out.write("| Stage | Count | Removed |\n")
        out.write("|-------|-------|--------|\n")
        out.write(f"| Matched hits | {len(matched_hits):,} | - |\n")
        out.write(f"| Synteny + conservation | {len(base_filtered):,} | {len(matched_hits)-len(base_filtered):,} |\n")
        out.write(f"| CDS exclusion | {len(cds_filtered):,} | {len(base_filtered)-len(cds_filtered):,} |\n")
        out.write(f"| Complexity filter | {len(complexity_filtered):,} | {len(cds_filtered)-len(complexity_filtered):,} |\n")
        out.write(f"| Unique loci | {len(collapsed_loci):,} | - |\n")
        out.write(f"| **Clean loci (final)** | **{len(clean_loci):,}** | - |\n")
        out.write(f"| Suspicious (multi-family) | {len(suspicious_loci):,} | - |\n\n")

        out.write("## Output Files\n\n")
        out.write(f"- `ancient_te_candidates_clean.tsv`: {len(clean_loci):,} high-confidence candidates\n")
        out.write(f"- `ancient_te_candidates_suspicious.tsv`: {len(suspicious_loci):,} multi-family loci for review\n")

    print(f"  Saved filtering report to {report_output.name}", file=sys.stderr)
    print("\nDone!", file=sys.stderr)


if __name__ == '__main__':
    main()
