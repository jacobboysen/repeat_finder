#!/usr/bin/env python3
"""
Analyze which TE structural regions (LTR, coding, etc.) are matched by UTR sequences.

This script:
1. Parses TE structural annotations from GFF
2. BLASTs UTR sequences against TE consensus database
3. Maps hits to TE structural regions
4. Reports enrichment of hits in different regions
"""

import argparse
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import csv


def parse_te_gff(gff_path):
    """
    Parse TE structural annotations from GFF file.

    Returns dict: te_id -> {region_type -> [(start, end), ...]}
    """
    te_annotations = defaultdict(lambda: defaultdict(list))
    te_lengths = {}

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('##sequence-region'):
                parts = line.split()
                te_id = parts[1]
                te_lengths[te_id] = int(parts[3])
            elif not line.startswith('#') and line:
                fields = line.split('\t')
                if len(fields) >= 9:
                    te_id = fields[0]
                    feature_type = fields[2]
                    start = int(fields[3])
                    end = int(fields[4])

                    # Normalize feature types
                    if 'LTR' in feature_type:
                        region = 'LTR'
                    elif feature_type == 'CDS':
                        region = 'CDS'
                    elif 'transposon' in feature_type.lower():
                        region = 'full_element'
                        continue  # Don't count full element as a region
                    else:
                        region = 'other'  # TATA, primer_binding_site, etc.

                    te_annotations[te_id][region].append((start, end))

    return te_annotations, te_lengths


def parse_consensus_header(header):
    """
    Parse consensus sequence header to get family and class.
    Format: family#CLASS/Subclass
    """
    if '#' in header:
        family, class_info = header.split('#', 1)
        if '/' in class_info:
            te_class, subclass = class_info.split('/', 1)
        else:
            te_class, subclass = class_info, ''
    else:
        family = header
        te_class, subclass = 'Unknown', ''

    return family, te_class, subclass


def map_flybase_to_consensus(flybase_name):
    """
    Map FlyBase TE instance name to consensus family.
    FlyBase format: family{}location or family{}gene[location]
    """
    if '{}' in flybase_name:
        family = flybase_name.split('{}')[0]
    else:
        family = flybase_name
    return family


def run_blast(query_fasta, db_path, output_path, num_threads=4):
    """Run BLAST against TE consensus database."""
    # Find blastn in PATH
    blastn = shutil.which('blastn') or 'blastn'
    cmd = [
        blastn,
        '-query', query_fasta,
        '-db', db_path,
        '-word_size', '7',
        '-gapopen', '2',
        '-gapextend', '1',
        '-penalty', '-1',
        '-reward', '1',
        '-dust', 'yes',
        '-evalue', '10',
        '-max_target_seqs', '50',
        '-max_hsps', '10',
        '-num_threads', str(num_threads),
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
        '-out', output_path
    ]

    print(f"Running BLAST: {' '.join(cmd[:6])}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"BLAST error: {result.stderr}", file=sys.stderr)
        return False
    return True


def classify_hit_region(hit_start, hit_end, te_annotations, te_id):
    """
    Classify which structural region(s) a hit overlaps.
    Returns list of (region_type, overlap_bp).
    """
    overlaps = []

    # Map FlyBase ID pattern to consensus ID pattern
    # GFF uses FBte format, consensus uses family#class format
    # For now, we'll handle this at the analysis stage

    for region_type, regions in te_annotations.get(te_id, {}).items():
        for reg_start, reg_end in regions:
            # Calculate overlap
            overlap_start = max(hit_start, reg_start)
            overlap_end = min(hit_end, reg_end)

            if overlap_start < overlap_end:
                overlap_bp = overlap_end - overlap_start
                overlaps.append((region_type, overlap_bp))

    return overlaps if overlaps else [('intergenic', hit_end - hit_start)]


def analyze_blast_results(blast_output, te_annotations, min_length=50):
    """
    Analyze BLAST results to count hits by TE region.
    """
    region_counts = defaultdict(int)
    region_bp = defaultdict(int)
    family_counts = defaultdict(int)
    class_counts = defaultdict(int)

    hits_by_gene = defaultdict(list)

    with open(blast_output) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 14:
                continue

            qseqid = row[0]
            sseqid = row[1]  # e.g., "mdg1#LTR/Gypsy"
            pident = float(row[2])
            length = int(row[3])
            sstart = int(row[8])
            send = int(row[9])

            if length < min_length:
                continue

            # Parse TE family and class from consensus header
            family, te_class, subclass = parse_consensus_header(sseqid)

            # Record hit
            family_counts[family] += 1
            class_counts[te_class] += 1

            # Determine hit position on TE
            hit_start = min(sstart, send)
            hit_end = max(sstart, send)

            # Map to GFF annotations (need to convert family name to FBte ID)
            # For now, categorize based on hit position relative to TE length
            hits_by_gene[qseqid].append({
                'family': family,
                'class': te_class,
                'subclass': subclass,
                'length': length,
                'pident': pident,
                'hit_start': hit_start,
                'hit_end': hit_end
            })

    return {
        'region_counts': dict(region_counts),
        'region_bp': dict(region_bp),
        'family_counts': dict(family_counts),
        'class_counts': dict(class_counts),
        'hits_by_gene': dict(hits_by_gene)
    }


def analyze_with_gff(blast_output, gff_path, min_length=50):
    """
    Analyze BLAST results using GFF structural annotations.
    Maps hit positions to LTR/CDS/other regions.
    """
    # Parse GFF
    te_annotations, te_lengths = parse_te_gff(gff_path)

    # Build mapping from family name to FBte ID
    # GFF uses FBte IDs, need to extract name from attributes
    family_to_fbte = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) >= 9 and 'name=' in fields[8]:
                te_id = fields[0]
                attrs = fields[8]
                for attr in attrs.split(';'):
                    if attr.startswith('name='):
                        name = attr.split('=')[1].replace('Dmel\\', '')
                        family_to_fbte[name.lower()] = te_id

    # Analyze hits
    region_stats = defaultdict(lambda: {'count': 0, 'bp': 0})
    family_region_stats = defaultdict(lambda: defaultdict(lambda: {'count': 0, 'bp': 0}))

    total_hits = 0
    hits_with_annotation = 0

    with open(blast_output) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 14:
                continue

            sseqid = row[1]
            length = int(row[3])
            sstart = int(row[8])
            send = int(row[9])

            if length < min_length:
                continue

            total_hits += 1

            # Parse family from consensus header
            family, te_class, _ = parse_consensus_header(sseqid)
            hit_start = min(sstart, send)
            hit_end = max(sstart, send)

            # Try to find corresponding FBte annotation
            fbte_id = family_to_fbte.get(family.lower())

            if fbte_id and fbte_id in te_annotations:
                hits_with_annotation += 1
                overlaps = classify_hit_region(hit_start, hit_end, te_annotations, fbte_id)

                for region_type, overlap_bp in overlaps:
                    region_stats[region_type]['count'] += 1
                    region_stats[region_type]['bp'] += overlap_bp
                    family_region_stats[family][region_type]['count'] += 1
                    family_region_stats[family][region_type]['bp'] += overlap_bp
            else:
                # No annotation available - categorize by position
                region_stats['unannotated']['count'] += 1
                region_stats['unannotated']['bp'] += length

    return {
        'total_hits': total_hits,
        'hits_with_annotation': hits_with_annotation,
        'region_stats': dict(region_stats),
        'family_region_stats': {k: dict(v) for k, v in family_region_stats.items()}
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query', '-q', required=True, help='Query FASTA (UTR sequences)')
    parser.add_argument('--blast-output', '-b', help='Pre-existing BLAST output (skip BLAST)')
    parser.add_argument('--db', default='data/blastdb/dmel_te_consensus', help='BLAST database')
    parser.add_argument('--gff', default='data/references/te_annotations.gff', help='TE structural annotations GFF')
    parser.add_argument('--output', '-o', default='results/te_region_analysis.tsv', help='Output file')
    parser.add_argument('--min-length', type=int, default=50, help='Minimum hit length')
    parser.add_argument('--threads', '-t', type=int, default=4, help='BLAST threads')

    args = parser.parse_args()

    # Run BLAST if needed
    if args.blast_output:
        blast_output = args.blast_output
    else:
        blast_output = Path(args.output).with_suffix('.blast.tsv')
        print(f"Running BLAST against TE consensus database...")
        if not run_blast(args.query, args.db, str(blast_output), args.threads):
            sys.exit(1)

    print(f"\nAnalyzing BLAST results...")

    # Basic analysis without GFF
    results = analyze_blast_results(blast_output, {}, args.min_length)

    print(f"\n=== TE Class Distribution ===")
    for te_class, count in sorted(results['class_counts'].items(), key=lambda x: -x[1]):
        print(f"  {te_class}: {count} hits")

    print(f"\n=== Top 20 TE Families ===")
    for family, count in sorted(results['family_counts'].items(), key=lambda x: -x[1])[:20]:
        print(f"  {family}: {count} hits")

    # Detailed analysis with GFF
    if Path(args.gff).exists():
        print(f"\nAnalyzing TE structural regions from GFF...")
        gff_results = analyze_with_gff(blast_output, args.gff, args.min_length)

        print(f"\n=== TE Region Distribution ===")
        print(f"Total hits: {gff_results['total_hits']}")
        print(f"Hits with structural annotation: {gff_results['hits_with_annotation']}")

        for region, stats in sorted(gff_results['region_stats'].items(), key=lambda x: -x[1]['count']):
            pct = 100 * stats['count'] / gff_results['total_hits'] if gff_results['total_hits'] > 0 else 0
            print(f"  {region}: {stats['count']} hits ({pct:.1f}%), {stats['bp']} bp")

    # Write detailed output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['metric', 'category', 'count', 'bp'])

        for te_class, count in results['class_counts'].items():
            writer.writerow(['te_class', te_class, count, ''])

        for family, count in results['family_counts'].items():
            writer.writerow(['te_family', family, count, ''])

        if Path(args.gff).exists():
            for region, stats in gff_results['region_stats'].items():
                writer.writerow(['te_region', region, stats['count'], stats['bp']])

    print(f"\nResults written to: {output_path}")


if __name__ == '__main__':
    main()
