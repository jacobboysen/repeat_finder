#!/usr/bin/env python3
"""Convert ALL BLAST TE hits to genomic BED format for conservation analysis."""

import sys

# Load UTR coordinates from GFF
print("Loading UTR coordinates...", file=sys.stderr)
utr_coords = {}
with open('/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel-all-r6.66.gff') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'three_prime_UTR':
            continue

        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]

        # Extract transcript ID
        attrs = parts[8]
        parent = None
        for attr in attrs.split(';'):
            if attr.startswith('Parent='):
                parent = attr.split('=')[1]
                break

        if parent:
            # Store by transcript ID (may have multiple UTR segments)
            if parent not in utr_coords:
                utr_coords[parent] = []
            utr_coords[parent].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })

print(f"Loaded UTR coords for {len(utr_coords)} transcripts", file=sys.stderr)

# Process BLAST hits in chunks to manage memory
def process_hits_file(filepath, output_file, prefix, sample_rate=1):
    """Process a BLAST hits file and write BED output.

    sample_rate: 1 = all hits, 10 = every 10th hit, etc.
    """
    count = 0
    written = 0
    skipped_no_utr = 0

    with open(filepath) as f:
        header = next(f)  # skip header
        for i, line in enumerate(f):
            if sample_rate > 1 and i % sample_rate != 0:
                continue

            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue

            # Format: qseqid, qstart, qend, sseqid, pident, length, evalue
            transcript_id = parts[0]
            qstart = int(parts[1])
            qend = int(parts[2])
            te_id = parts[3]
            pident = float(parts[4])
            length = int(parts[5])

            # Get UTR coordinates
            if transcript_id not in utr_coords:
                skipped_no_utr += 1
                continue

            # Use first UTR segment (simplified - most transcripts have one)
            utr = utr_coords[transcript_id][0]

            # Convert UTR-relative to genomic coordinates
            if utr['strand'] == '+':
                genomic_start = utr['start'] + qstart - 1
                genomic_end = utr['start'] + qend - 1
            else:
                # Minus strand
                genomic_start = utr['end'] - qend + 1
                genomic_end = utr['end'] - qstart + 1

            # BED is 0-based, half-open
            bed_start = genomic_start - 1
            bed_end = genomic_end

            # Add chr prefix
            chrom = utr['chrom']
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom

            # Create unique name with key info
            name = f"{prefix}_{written}|{transcript_id}|{te_id}|{pident:.1f}|{length}"

            output_file.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t0\t{utr['strand']}\n")
            written += 1

            if written % 100000 == 0:
                print(f"  Processed {written} hits...", file=sys.stderr)

    return written, skipped_no_utr

# Process both files
print("\nProcessing novel hits...", file=sys.stderr)
with open('/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/te_hits_all_genomic.bed', 'w') as out:
    novel_count, novel_skip = process_hits_file(
        '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/blast_hits_novel.tsv',
        out, 'novel'
    )
    print(f"Novel: {novel_count} written, {novel_skip} skipped (no UTR)", file=sys.stderr)

    print("\nProcessing known hits...", file=sys.stderr)
    known_count, known_skip = process_hits_file(
        '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/blast_hits_known_repeatmasker.tsv',
        out, 'known'
    )
    print(f"Known: {known_count} written, {known_skip} skipped (no UTR)", file=sys.stderr)

print(f"\nTotal: {novel_count + known_count} hits written to BED file", file=sys.stderr)
