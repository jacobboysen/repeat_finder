#!/usr/bin/env python3
"""
Build TE loci per UTR isoform with annotations (no filtering).

Each transcript's 3'UTR is a distinct model. TE hits are collapsed into loci
within each transcript, then annotated with:
- CDS overlap status (in any isoform)
- TE families involved
- Conservation (phyloP)
- Synteny

Naming convention: {FBtr}|{chrom}|{start}-{end}
  Example: FBtr0080594|chr2L|14089270-14089329

Usage:
    python scripts/build_utr_te_loci.py
"""

import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir, get_references_dir


def load_cds_regions(gff_path: str) -> dict:
    """Load all CDS regions indexed by chromosome."""
    cds_regions = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            chrom = 'chr' + parts[0]
            start = int(parts[3])
            end = int(parts[4])
            # Extract parent transcript
            parent = None
            for attr in parts[8].split(';'):
                if attr.startswith('Parent='):
                    parent = attr.split('=')[1]
                    break
            cds_regions[chrom].append((start, end, parent))

    for chrom in cds_regions:
        cds_regions[chrom].sort()
    return cds_regions


def get_cds_overlaps(chrom: str, start: int, end: int, cds_regions: dict) -> list:
    """Return list of transcripts whose CDS overlaps this region."""
    overlapping = []
    if chrom not in cds_regions:
        return overlapping

    for cds_start, cds_end, parent in cds_regions[chrom]:
        if cds_start > end:
            break
        if max(start, cds_start) < min(end, cds_end):
            if parent and parent not in overlapping:
                overlapping.append(parent)
    return overlapping


def load_te_families(te_fasta_path: str) -> dict:
    """Load TE ID to family mapping."""
    families = {}
    with open(te_fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split()
                te_id = parts[0][1:]
                for part in parts:
                    if part.startswith('name='):
                        name = part.split('=')[1].split('{')[0]
                        families[te_id] = name
                        break
    return families


def load_synteny_data(filepath: str) -> dict:
    """Load synteny data."""
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
    """Load conservation scores."""
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
                    'phyloP': phyloP,
                    'pident': float(name_parts[3]) if len(name_parts) > 3 else 0,
                    'length': int(name_parts[4]) if len(name_parts) > 4 else 0
                })
    return data


def collapse_hits_to_loci(hits: list, max_gap: int = 10) -> list:
    """
    Collapse overlapping/adjacent hits into loci within each transcript.
    Preserves all TE families that hit each locus.
    """
    if not hits:
        return []

    # Group by transcript
    by_transcript = defaultdict(list)
    for hit in hits:
        by_transcript[hit['transcript']].append(hit)

    all_loci = []

    for transcript, transcript_hits in by_transcript.items():
        # Sort by start position
        sorted_hits = sorted(transcript_hits, key=lambda x: x['start'])

        loci = []
        current = None

        for hit in sorted_hits:
            if current is None:
                current = {
                    'transcript': transcript,
                    'chrom': hit['chrom'],
                    'start': hit['start'],
                    'end': hit['end'],
                    'hits': [hit],
                    'te_ids': {hit['te_id']},
                    'te_families': {hit['te_family']},
                    'best_phyloP': hit['phyloP'],
                    'best_pident': hit['pident'],
                    'max_length': hit['length'],
                    'syntenic_species': hit['syntenic_species']
                }
            elif hit['start'] <= current['end'] + max_gap:
                # Merge
                current['end'] = max(current['end'], hit['end'])
                current['hits'].append(hit)
                current['te_ids'].add(hit['te_id'])
                current['te_families'].add(hit['te_family'])
                if hit['phyloP'] > current['best_phyloP']:
                    current['best_phyloP'] = hit['phyloP']
                    current['best_pident'] = hit['pident']
                current['max_length'] = max(current['max_length'], hit['length'])
                current['syntenic_species'] = max(current['syntenic_species'], hit['syntenic_species'])
            else:
                loci.append(current)
                current = {
                    'transcript': transcript,
                    'chrom': hit['chrom'],
                    'start': hit['start'],
                    'end': hit['end'],
                    'hits': [hit],
                    'te_ids': {hit['te_id']},
                    'te_families': {hit['te_family']},
                    'best_phyloP': hit['phyloP'],
                    'best_pident': hit['pident'],
                    'max_length': hit['length'],
                    'syntenic_species': hit['syntenic_species']
                }

        if current:
            loci.append(current)

        all_loci.extend(loci)

    return all_loci


def main():
    project_root = get_project_root()
    results_dir = get_results_dir() / 'repeatmasker_analysis'
    output_dir = get_results_dir() / 'utr_te_loci'
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*70, file=sys.stderr)
    print("BUILD UTR TE LOCI (per isoform, with annotations)", file=sys.stderr)
    print("="*70, file=sys.stderr)

    # Load reference data
    print("\n[1/5] Loading CDS regions...", file=sys.stderr)
    gff_path = get_references_dir() / 'dmel-all-r6.66.gff'
    cds_regions = load_cds_regions(str(gff_path))
    total_cds = sum(len(v) for v in cds_regions.values())
    print(f"  Loaded {total_cds:,} CDS regions", file=sys.stderr)

    print("\n[2/5] Loading TE families...", file=sys.stderr)
    te_fasta = get_references_dir() / 'dmel_te_flybase.fasta'
    te_families = load_te_families(str(te_fasta))
    print(f"  Loaded {len(te_families):,} TE family mappings", file=sys.stderr)

    print("\n[3/5] Loading synteny data...", file=sys.stderr)
    synteny = load_synteny_data(str(results_dir / 'te_hits_all_synteny_sampled.tsv'))
    print(f"  Loaded {len(synteny):,} hits", file=sys.stderr)

    print("\n[4/5] Loading conservation data...", file=sys.stderr)
    conservation = load_conservation_data(str(results_dir / 'te_hits_all_conservation.tab'))
    print(f"  Loaded conservation for {len(conservation):,} pairs", file=sys.stderr)

    # Match and build hits
    print("\n[5/5] Building matched hits...", file=sys.stderr)
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
                    'phyloP': best_cons['phyloP'],
                    'syntenic_species': syn_row['syntenic_species'],
                    'category': syn_row['category']
                }
                matched_hits.append(hit)

    print(f"  Matched {len(matched_hits):,} hits", file=sys.stderr)

    # Filter for ancient candidates (syntenic + conserved)
    ancient_hits = [h for h in matched_hits if h['syntenic_species'] >= 2 and h['phyloP'] > 1]
    print(f"  Ancient candidates (syntenic + conserved): {len(ancient_hits):,}", file=sys.stderr)

    # Collapse to loci per transcript
    print("\n  Collapsing to loci per transcript...", file=sys.stderr)
    loci = collapse_hits_to_loci(ancient_hits)
    print(f"  Total loci: {len(loci):,}", file=sys.stderr)
    print(f"  Unique transcripts: {len(set(l['transcript'] for l in loci)):,}", file=sys.stderr)

    # Annotate with CDS overlap
    print("\n  Annotating CDS overlaps...", file=sys.stderr)
    cds_overlap_count = 0
    for locus in loci:
        overlaps = get_cds_overlaps(locus['chrom'], locus['start'], locus['end'], cds_regions)
        locus['cds_overlap'] = len(overlaps) > 0
        locus['cds_overlap_transcripts'] = overlaps
        if locus['cds_overlap']:
            cds_overlap_count += 1

    print(f"  Loci with CDS overlap: {cds_overlap_count:,} ({100*cds_overlap_count/len(loci):.1f}%)", file=sys.stderr)

    # Build locus ID: {FBtr}|{chrom}|{start}-{end}
    for locus in loci:
        locus['locus_id'] = f"{locus['transcript']}|{locus['chrom']}|{locus['start']}-{locus['end']}"

    # Sort by phyloP
    loci.sort(key=lambda x: -x['best_phyloP'])

    # Save results
    output_file = output_dir / 'ancient_te_loci_by_isoform.tsv'
    print(f"\n  Saving to {output_file}...", file=sys.stderr)

    with open(output_file, 'w') as out:
        headers = [
            'locus_id', 'transcript', 'chrom', 'start', 'end', 'locus_length',
            'n_te_families', 'te_families', 'n_te_hits', 'te_ids',
            'best_phyloP', 'best_pident', 'max_hit_length', 'syntenic_species',
            'cds_overlap', 'cds_overlap_transcripts'
        ]
        out.write('\t'.join(headers) + '\n')

        for locus in loci:
            row = [
                locus['locus_id'],
                locus['transcript'],
                locus['chrom'],
                str(locus['start']),
                str(locus['end']),
                str(locus['end'] - locus['start']),
                str(len(locus['te_families'])),
                ','.join(sorted(locus['te_families'])),
                str(len(locus['hits'])),
                ','.join(sorted(locus['te_ids'])),
                f"{locus['best_phyloP']:.4f}",
                f"{locus['best_pident']:.1f}",
                str(locus['max_length']),
                str(locus['syntenic_species']),
                'yes' if locus['cds_overlap'] else 'no',
                ','.join(locus['cds_overlap_transcripts']) if locus['cds_overlap'] else ''
            ]
            out.write('\t'.join(row) + '\n')

    # Summary stats
    print("\n" + "="*70, file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print("="*70, file=sys.stderr)
    print(f"Total loci: {len(loci):,}", file=sys.stderr)
    print(f"Unique transcripts (UTR isoforms): {len(set(l['transcript'] for l in loci)):,}", file=sys.stderr)
    print(f"Loci with CDS overlap: {cds_overlap_count:,} ({100*cds_overlap_count/len(loci):.1f}%)", file=sys.stderr)

    single_family = sum(1 for l in loci if len(l['te_families']) == 1)
    multi_family = sum(1 for l in loci if len(l['te_families']) > 1)
    print(f"Single TE family: {single_family:,} ({100*single_family/len(loci):.1f}%)", file=sys.stderr)
    print(f"Multiple TE families: {multi_family:,} ({100*multi_family/len(loci):.1f}%)", file=sys.stderr)

    print(f"\nphyloP range: {min(l['best_phyloP'] for l in loci):.2f} - {max(l['best_phyloP'] for l in loci):.2f}", file=sys.stderr)

    print(f"\nOutput: {output_file}", file=sys.stderr)
    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
