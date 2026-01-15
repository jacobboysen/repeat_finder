#!/usr/bin/env python3
"""
Analyze overlap between RepeatMasker annotations and UTR regions.

This script:
1. Extracts 3'UTR coordinates from FlyBase GFF
2. Parses RepeatMasker .out file
3. Finds RepeatMasker hits that overlap UTR regions
4. Cross-references with BLAST TE hits to classify as "known" vs "novel"
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

# Add scripts directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_references_dir, get_results_dir
from utils.blast_io import load_blast_results


def parse_repeatmasker_out(rm_file):
    """
    Parse RepeatMasker .out file.

    Returns list of dicts with repeat annotations.
    """
    repeats = []

    with open(rm_file) as f:
        # Skip header lines (first 3 lines)
        for _ in range(3):
            next(f)

        for line in f:
            line = line.strip()
            if not line:
                continue

            # RepeatMasker output is space-delimited with variable spacing
            parts = line.split()
            if len(parts) < 14:
                continue

            try:
                # Handle the (left) field which may have parentheses
                score = int(parts[0])
                pct_div = float(parts[1])
                pct_del = float(parts[2])
                pct_ins = float(parts[3])
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                # parts[7] is (remaining) - skip
                strand = '+' if parts[8] == '+' else '-'
                repeat_name = parts[9]
                repeat_class = parts[10]

                # Normalize chromosome names (chr2L -> 2L)
                if chrom.startswith('chr'):
                    chrom = chrom[3:]

                repeats.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'repeat_name': repeat_name,
                    'repeat_class': repeat_class,
                    'score': score,
                    'pct_divergence': pct_div,
                    'pct_deletion': pct_del,
                    'pct_insertion': pct_ins,
                })
            except (ValueError, IndexError):
                continue

    return repeats


def parse_utr_coordinates(gff_file, utr_type='three_prime_UTR'):
    """
    Extract UTR coordinates from FlyBase GFF.

    Returns dict mapping (chrom, start, end) -> list of transcript IDs
    """
    utrs = defaultdict(list)
    utr_by_transcript = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != utr_type:
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Extract transcript ID(s) from Parent attribute
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            if parent_match:
                transcript_ids = parent_match.group(1).split(',')
                for tr_id in transcript_ids:
                    utrs[(chrom, start, end, strand)].append(tr_id)
                    utr_by_transcript[tr_id] = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }

    return utrs, utr_by_transcript


def build_interval_index(utrs):
    """
    Build chromosome-indexed interval structure for fast overlap queries.
    """
    by_chrom = defaultdict(list)

    for (chrom, start, end, strand), transcripts in utrs.items():
        by_chrom[chrom].append({
            'start': start,
            'end': end,
            'strand': strand,
            'transcripts': transcripts
        })

    # Sort by start position for binary search
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x['start'])

    return by_chrom


def find_overlaps(repeats, utr_index):
    """
    Find RepeatMasker hits that overlap UTR regions.
    """
    overlaps = []

    for repeat in repeats:
        chrom = repeat['chrom']
        r_start = repeat['start']
        r_end = repeat['end']

        if chrom not in utr_index:
            continue

        # Check all UTRs on this chromosome
        for utr in utr_index[chrom]:
            u_start = utr['start']
            u_end = utr['end']

            # Check for overlap
            if r_start <= u_end and r_end >= u_start:
                # Calculate overlap region
                overlap_start = max(r_start, u_start)
                overlap_end = min(r_end, u_end)
                overlap_len = overlap_end - overlap_start + 1

                # Calculate position within UTR
                utr_len = u_end - u_start + 1
                utr_position = overlap_start - u_start  # 0-based position in UTR

                overlaps.append({
                    'chrom': chrom,
                    'repeat_start': r_start,
                    'repeat_end': r_end,
                    'utr_start': u_start,
                    'utr_end': u_end,
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_len': overlap_len,
                    'utr_position': utr_position,
                    'utr_len': utr_len,
                    'utr_strand': utr['strand'],
                    'repeat_name': repeat['repeat_name'],
                    'repeat_class': repeat['repeat_class'],
                    'repeat_score': repeat['score'],
                    'repeat_divergence': repeat['pct_divergence'],
                    'transcripts': utr['transcripts'],
                })

    return overlaps


def summarize_repeat_classes(overlaps):
    """
    Summarize repeat classes found in UTRs.
    """
    class_counts = defaultdict(int)
    class_bp = defaultdict(int)

    for overlap in overlaps:
        class_counts[overlap['repeat_class']] += 1
        class_bp[overlap['repeat_class']] += overlap['overlap_len']

    return class_counts, class_bp


def compare_with_blast_hits(overlaps, blast_df, utr_by_transcript):
    """
    Compare RepeatMasker overlaps with BLAST TE hits.

    Identifies which BLAST hits correspond to known RepeatMasker regions.
    """
    # Index RepeatMasker overlaps by transcript
    rm_by_transcript = defaultdict(list)
    for overlap in overlaps:
        for tr_id in overlap['transcripts']:
            rm_by_transcript[tr_id].append(overlap)

    # Compare each BLAST hit
    known_hits = []
    novel_hits = []

    for _, hit in blast_df.iterrows():
        qseqid = hit['qseqid']
        qstart = hit['qstart']
        qend = hit['qend']

        # Check if this BLAST hit overlaps any RepeatMasker region
        is_known = False
        matching_rm = None

        if qseqid in rm_by_transcript:
            for rm_overlap in rm_by_transcript[qseqid]:
                # BLAST positions are relative to UTR, so compare directly
                rm_utr_start = rm_overlap['utr_position']
                rm_utr_end = rm_utr_start + rm_overlap['overlap_len'] - 1

                # Check for overlap
                if qstart <= rm_utr_end and qend >= rm_utr_start:
                    is_known = True
                    matching_rm = rm_overlap
                    break

        hit_info = {
            'qseqid': qseqid,
            'qstart': qstart,
            'qend': qend,
            'sseqid': hit['sseqid'],
            'pident': hit['pident'],
            'length': hit['length'],
            'evalue': hit['evalue'],
        }

        if is_known:
            hit_info['rm_repeat_name'] = matching_rm['repeat_name']
            hit_info['rm_repeat_class'] = matching_rm['repeat_class']
            hit_info['rm_divergence'] = matching_rm['repeat_divergence']
            known_hits.append(hit_info)
        else:
            novel_hits.append(hit_info)

    return known_hits, novel_hits


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--rm-file',
        type=Path,
        default=get_references_dir() / 'dm6.fa.out',
        help='RepeatMasker .out file'
    )
    parser.add_argument(
        '--gff-file',
        type=Path,
        default=get_references_dir() / 'dmel-all-r6.66.gff',
        help='FlyBase GFF annotation file'
    )
    parser.add_argument(
        '--blast-file',
        type=Path,
        help='BLAST results TSV to compare (optional)'
    )
    parser.add_argument(
        '--utr-type',
        choices=['three_prime_UTR', 'five_prime_UTR'],
        default='three_prime_UTR',
        help='UTR type to analyze'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=get_results_dir() / 'repeatmasker_analysis',
        help='Output directory'
    )
    parser.add_argument(
        '--min-overlap',
        type=int,
        default=20,
        help='Minimum overlap length to report'
    )

    args = parser.parse_args()

    print("=" * 70)
    print("RepeatMasker / UTR Overlap Analysis")
    print("=" * 70)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Parse RepeatMasker file
    print(f"\nParsing RepeatMasker file: {args.rm_file}")
    repeats = parse_repeatmasker_out(args.rm_file)
    print(f"  Loaded {len(repeats):,} repeat annotations")

    # Count repeat classes
    rm_classes = defaultdict(int)
    for r in repeats:
        rm_classes[r['repeat_class']] += 1
    print(f"\n  Top repeat classes in genome:")
    for cls, count in sorted(rm_classes.items(), key=lambda x: -x[1])[:10]:
        print(f"    {cls}: {count:,}")

    # Parse UTR coordinates
    print(f"\nParsing UTR coordinates from: {args.gff_file}")
    utrs, utr_by_transcript = parse_utr_coordinates(args.gff_file, args.utr_type)
    print(f"  Found {len(utrs):,} unique {args.utr_type} regions")
    print(f"  Covering {len(utr_by_transcript):,} transcripts")

    # Build interval index
    print("\nBuilding interval index...")
    utr_index = build_interval_index(utrs)

    # Find overlaps
    print("\nFinding RepeatMasker/UTR overlaps...")
    overlaps = find_overlaps(repeats, utr_index)

    # Filter by minimum overlap
    overlaps = [o for o in overlaps if o['overlap_len'] >= args.min_overlap]
    print(f"  Found {len(overlaps):,} overlaps (>= {args.min_overlap}bp)")

    # Summarize repeat classes in UTRs
    class_counts, class_bp = summarize_repeat_classes(overlaps)

    print(f"\n  Repeat classes found in {args.utr_type}s:")
    print(f"  {'Class':<25} {'Count':>10} {'Total BP':>12}")
    print(f"  {'-'*25} {'-'*10} {'-'*12}")
    for cls, count in sorted(class_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {cls:<25} {count:>10,} {class_bp[cls]:>12,}")

    # Count unique transcripts with repeats
    transcripts_with_repeats = set()
    for o in overlaps:
        transcripts_with_repeats.update(o['transcripts'])
    print(f"\n  Transcripts with RepeatMasker hits: {len(transcripts_with_repeats):,}")
    print(f"  Percentage of UTRs with repeats: {100*len(transcripts_with_repeats)/len(utr_by_transcript):.1f}%")

    # Save overlaps to TSV
    overlap_file = args.output_dir / f'{args.utr_type}_repeatmasker_overlaps.tsv'
    with open(overlap_file, 'w') as f:
        f.write('chrom\tutr_start\tutr_end\trepeat_start\trepeat_end\toverlap_len\t'
                'repeat_name\trepeat_class\tdivergence\ttranscripts\n')
        for o in sorted(overlaps, key=lambda x: -x['overlap_len'])[:10000]:
            transcripts = ','.join(o['transcripts'][:3])  # Limit transcript list
            f.write(f"{o['chrom']}\t{o['utr_start']}\t{o['utr_end']}\t"
                    f"{o['repeat_start']}\t{o['repeat_end']}\t{o['overlap_len']}\t"
                    f"{o['repeat_name']}\t{o['repeat_class']}\t{o['repeat_divergence']:.1f}\t"
                    f"{transcripts}\n")
    print(f"\n  Saved overlaps to: {overlap_file}")

    # Compare with BLAST hits if provided
    if args.blast_file and args.blast_file.exists():
        print(f"\n{'='*70}")
        print("Comparing with BLAST TE hits")
        print("=" * 70)

        print(f"\nLoading BLAST results: {args.blast_file}")
        blast_df = load_blast_results(args.blast_file)
        print(f"  Loaded {len(blast_df):,} BLAST hits")

        known_hits, novel_hits = compare_with_blast_hits(
            overlaps, blast_df, utr_by_transcript
        )

        print(f"\n  BLAST hits in known RepeatMasker regions: {len(known_hits):,} ({100*len(known_hits)/len(blast_df):.1f}%)")
        print(f"  BLAST hits NOT in RepeatMasker regions: {len(novel_hits):,} ({100*len(novel_hits)/len(blast_df):.1f}%)")

        # Breakdown of known hits by RM class
        if known_hits:
            known_by_class = defaultdict(int)
            for h in known_hits:
                known_by_class[h['rm_repeat_class']] += 1

            print(f"\n  Known hits by RepeatMasker class:")
            for cls, count in sorted(known_by_class.items(), key=lambda x: -x[1])[:10]:
                print(f"    {cls}: {count:,}")

        # Save known/novel hits
        known_file = args.output_dir / 'blast_hits_known_repeatmasker.tsv'
        novel_file = args.output_dir / 'blast_hits_novel.tsv'

        with open(known_file, 'w') as f:
            f.write('qseqid\tqstart\tqend\tsseqid\tpident\tlength\tevalue\t'
                    'rm_repeat_name\trm_repeat_class\trm_divergence\n')
            for h in known_hits:
                f.write(f"{h['qseqid']}\t{h['qstart']}\t{h['qend']}\t{h['sseqid']}\t"
                        f"{h['pident']:.1f}\t{h['length']}\t{h['evalue']:.2e}\t"
                        f"{h['rm_repeat_name']}\t{h['rm_repeat_class']}\t{h['rm_divergence']:.1f}\n")

        with open(novel_file, 'w') as f:
            f.write('qseqid\tqstart\tqend\tsseqid\tpident\tlength\tevalue\n')
            for h in novel_hits:
                f.write(f"{h['qseqid']}\t{h['qstart']}\t{h['qend']}\t{h['sseqid']}\t"
                        f"{h['pident']:.1f}\t{h['length']}\t{h['evalue']:.2e}\n")

        print(f"\n  Saved known hits to: {known_file}")
        print(f"  Saved novel hits to: {novel_file}")

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
