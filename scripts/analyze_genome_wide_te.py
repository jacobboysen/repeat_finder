#!/usr/bin/env python3
"""
Genome-wide TE analysis: gene-level aggregation, HTML visualizations, and strand bias.
"""

import argparse
from pathlib import Path
from collections import defaultdict
import re


def parse_3utr_fasta(fasta_path):
    """Parse 3'UTR FASTA to get transcript lengths and gene mappings."""
    transcript_info = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                # >FBtr0070000 type=three_prime_UTR; loc=...; length=1019; parent=FBgn0031081; ...
                parts = line.strip().split()
                fbtr = parts[0][1:]  # Remove >
                length = None
                parent = None
                for part in parts:
                    if part.startswith('length='):
                        length = int(part.split('=')[1].rstrip(';'))
                    if part.startswith('parent='):
                        parent = part.split('=')[1].rstrip(';')
                if fbtr and length and parent:
                    transcript_info[fbtr] = {'length': length, 'fbgn': parent}
    return transcript_info


def parse_te_database(fasta_path):
    """Parse TE database FASTA to get TE info (family, class)."""
    te_info = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                # >FBti0019046 roo{}412
                # or from dfam: >Gypsy-2_DM#LTR/Gypsy
                parts = line.strip().split()
                te_id = parts[0][1:]
                te_name = parts[1] if len(parts) > 1 else te_id

                # Extract TE class from name
                te_class = 'Unknown'
                if '#' in te_name:
                    class_part = te_name.split('#')[1]
                    if '/' in class_part:
                        te_class = class_part.split('/')[0]
                    else:
                        te_class = class_part
                elif '{}' in te_name:
                    # FlyBase format like roo{}412
                    te_class = 'LTR'  # Most FlyBase TEs are LTR

                te_info[te_id] = {'name': te_name, 'class': te_class}
    return te_info


def parse_blast_results(blast_path, transcript_info):
    """Parse genome-wide BLAST results and aggregate by transcript."""
    transcript_hits = defaultdict(list)

    print(f"Parsing BLAST results from {blast_path}...")
    with open(blast_path) as f:
        for i, line in enumerate(f):
            if i % 500000 == 0 and i > 0:
                print(f"  Processed {i:,} lines...")

            parts = line.strip().split('\t')
            if len(parts) < 14:  # Minimum: qseqid through slen (14 columns)
                continue

            fbtr = parts[0]
            te_id = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            qstart = int(parts[6])
            qend = int(parts[7])
            sstart = int(parts[8])
            send = int(parts[9])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qseq = parts[14] if len(parts) > 14 else ''
            sseq = parts[15] if len(parts) > 15 else ''

            # Determine strand: if sstart > send, it's antisense
            strand = '+' if sstart < send else '-'

            transcript_hits[fbtr].append({
                'te_id': te_id,
                'pident': pident,
                'length': length,
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'evalue': evalue,
                'bitscore': bitscore,
                'strand': strand,
                'qseq': qseq,
                'sseq': sseq
            })

    print(f"  Total transcripts with hits: {len(transcript_hits):,}")
    return transcript_hits


def aggregate_by_gene(transcript_hits, transcript_info):
    """Aggregate transcript-level data to gene level."""
    gene_data = defaultdict(lambda: {
        'transcripts': [],
        'total_utr_len': 0,
        'hits': [],
        'hit_count': 0,
        'hit_bp': 0
    })

    for fbtr, hits in transcript_hits.items():
        if fbtr not in transcript_info:
            continue

        fbgn = transcript_info[fbtr]['fbgn']
        utr_len = transcript_info[fbtr]['length']

        gene_data[fbgn]['transcripts'].append(fbtr)
        gene_data[fbgn]['total_utr_len'] += utr_len
        gene_data[fbgn]['hits'].extend(hits)
        gene_data[fbgn]['hit_count'] += len(hits)
        gene_data[fbgn]['hit_bp'] += sum(h['length'] for h in hits)

    # Calculate density (hit_bp per kb of total UTR)
    for fbgn, data in gene_data.items():
        if data['total_utr_len'] > 0:
            data['density'] = (data['hit_bp'] / data['total_utr_len']) * 1000
        else:
            data['density'] = 0

    return gene_data


def get_gene_symbols(fbgn_list):
    """Get gene symbols from FlyBase IDs using preloaded mapping."""
    # Load gene symbol mapping from FlyBase annotation
    symbol_map = {}

    # Try to load from common annotation sources
    gene_map_paths = [
        Path('/Users/jacobboysen/git_repos/repeat_finder/data/references/fbgn_annotation_ID.tsv'),
        Path('/Users/jacobboysen/git_repos/repeat_finder/data/references/gene_symbols.tsv'),
    ]

    for path in gene_map_paths:
        if path.exists():
            with open(path) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('##'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        # Try different column arrangements
                        for i, part in enumerate(parts):
                            if part.startswith('FBgn'):
                                # Look for symbol in adjacent columns
                                for j, p in enumerate(parts):
                                    if j != i and not p.startswith('FBgn') and not p.startswith('FBtr'):
                                        symbol_map[part] = p
                                        break

    return symbol_map


def write_top_bottom_100(gene_data, output_dir):
    """Write top and bottom 100 genes by TE density."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort genes by density
    sorted_genes = sorted(gene_data.items(), key=lambda x: x[1]['density'], reverse=True)

    # Get top 100 and bottom 100
    top_100 = sorted_genes[:100]
    bottom_100 = sorted_genes[-100:]

    # Get TE family counts for each gene
    def get_te_families(hits):
        te_counts = defaultdict(int)
        for h in hits:
            te_counts[h['te_id']] += 1
        # Return top 3
        sorted_tes = sorted(te_counts.items(), key=lambda x: x[1], reverse=True)[:3]
        return sorted_tes

    # Write top 100
    with open(output_dir / 'top_100_te_genes_FIXED.tsv', 'w') as f:
        f.write('rank\tfbgn\tdensity\thits\thit_bp\ttotal_utr_len\tnum_transcripts\ttop_te_1\ttop_te_2\ttop_te_3\n')
        for i, (fbgn, data) in enumerate(top_100, 1):
            top_tes = get_te_families(data['hits'])
            te_strs = [f"{te}:{count}" for te, count in top_tes]
            te_strs.extend([''] * (3 - len(te_strs)))
            f.write(f"{i}\t{fbgn}\t{data['density']:.1f}\t{data['hit_count']}\t{data['hit_bp']}\t"
                   f"{data['total_utr_len']}\t{len(data['transcripts'])}\t{te_strs[0]}\t{te_strs[1]}\t{te_strs[2]}\n")

    # Write bottom 100
    with open(output_dir / 'bottom_100_te_genes_FIXED.tsv', 'w') as f:
        f.write('rank\tfbgn\tdensity\thits\thit_bp\ttotal_utr_len\tnum_transcripts\ttop_te_1\ttop_te_2\ttop_te_3\n')
        for i, (fbgn, data) in enumerate(bottom_100, 1):
            top_tes = get_te_families(data['hits'])
            te_strs = [f"{te}:{count}" for te, count in top_tes]
            te_strs.extend([''] * (3 - len(te_strs)))
            rank = len(sorted_genes) - 99 + (i - 1)
            f.write(f"{rank}\t{fbgn}\t{data['density']:.1f}\t{data['hit_count']}\t{data['hit_bp']}\t"
                   f"{data['total_utr_len']}\t{len(data['transcripts'])}\t{te_strs[0]}\t{te_strs[1]}\t{te_strs[2]}\n")

    print(f"Wrote top/bottom 100 to {output_dir}")
    return top_100, bottom_100, sorted_genes


def generate_html_visualization(fbgn, gene_data, transcript_info, utr_sequences, te_info, output_path):
    """Generate HTML visualization for a gene's 3'UTR TE content."""

    data = gene_data[fbgn]
    hits = data['hits']

    # Get the longest transcript for visualization
    transcripts = data['transcripts']
    longest_tr = max(transcripts, key=lambda t: transcript_info.get(t, {}).get('length', 0))
    utr_len = transcript_info[longest_tr]['length']
    utr_seq = utr_sequences.get(longest_tr, 'N' * utr_len)

    # Filter hits for this transcript
    tr_hits = [h for h in hits if True]  # For gene-level, show all hits

    # Build position coverage map
    sense_pos = set()
    anti_pos = set()
    for h in tr_hits:
        for p in range(h['qstart'] - 1, h['qend']):
            if h['strand'] == '+':
                sense_pos.add(p)
            else:
                anti_pos.add(p)

    # Count strands
    sense_count = sum(1 for h in tr_hits if h['strand'] == '+')
    anti_count = sum(1 for h in tr_hits if h['strand'] == '-')

    # Get TE class distribution
    class_counts = defaultdict(int)
    for h in tr_hits:
        te_id = h['te_id']
        if te_id in te_info:
            class_counts[te_info[te_id]['class']] += 1
        else:
            class_counts['Unknown'] += 1

    # Generate HTML
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{fbgn} 3'UTR TE Annotation</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 14px;
    line-height: 1.6;
    margin: 40px;
    background: #fafafa;
    color: #333;
}}
h1, h2, h3 {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-weight: 600;
}}
h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
.sequence-container {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin: 20px 0;
    overflow-x: auto;
}}
.seq-line {{
    white-space: pre;
    font-family: "SF Mono", "Monaco", "Inconsolata", "Fira Mono", monospace;
    font-size: 13px;
}}
.pos-num {{ color: #888; margin-right: 10px; }}
.te-sense {{ background-color: #cce5ff; color: #004085; font-weight: 600; }}
.te-anti {{ background-color: #fff3cd; color: #856404; font-weight: 600; }}
.te-both {{ background-color: #d4edda; color: #155724; font-weight: 600; }}
.normal {{ color: #333; }}
.alignment-box {{
    background: white;
    padding: 15px;
    margin: 15px 0;
    border-radius: 4px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}}
.alignment-box.sense {{ border-left: 4px solid #004085; }}
.alignment-box.antisense {{ border-left: 4px solid #856404; }}
.align-header {{ font-weight: 600; margin-bottom: 8px; }}
.align-meta {{ color: #666; font-size: 12px; margin-bottom: 10px; }}
.align-seq {{
    white-space: pre;
    font-family: "SF Mono", "Monaco", "Inconsolata", "Fira Mono", monospace;
    font-size: 12px;
    line-height: 1.4;
}}
table {{ border-collapse: collapse; margin: 20px 0; width: 100%; max-width: 900px; }}
th, td {{ padding: 8px 12px; border: 1px solid #ddd; text-align: left; }}
th {{ background: #f5f5f5; font-weight: 600; }}
tr:nth-child(even) {{ background: #fafafa; }}
.summary-box {{
    background: #e8f4e8;
    padding: 15px 20px;
    border-radius: 8px;
    margin: 20px 0;
}}
.strand-box {{
    display: inline-block;
    padding: 15px 20px;
    border-radius: 8px;
    margin: 5px;
    text-align: center;
}}
.strand-sense {{ background: #cce5ff; }}
.strand-anti {{ background: #fff3cd; }}
.legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 15px;
    margin: 15px 0;
    padding: 10px;
    background: white;
    border-radius: 4px;
}}
.legend-item {{ display: flex; align-items: center; gap: 6px; font-size: 13px; }}
.legend-color {{ width: 16px; height: 16px; display: inline-block; border-radius: 3px; }}
.class-tag {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 3px;
    font-size: 11px;
    font-weight: 600;
    color: white;
}}
.class-LTR {{ background: #0066cc; }}
.class-LINE {{ background: #cc6600; }}
.class-DNA {{ background: #009900; }}
.class-Helitron {{ background: #990099; }}
.class-Unknown {{ background: #666666; }}
</style>
</head>
<body>

<h1>{fbgn} 3'UTR TE Annotation</h1>

<div class="summary-box">
<strong>Summary:</strong> {len(tr_hits)} TE hits | {data['hit_bp']:,} bp aligned | Density: {data['density']:.1f} bp/kb | UTR length: {data['total_utr_len']:,} bp ({len(transcripts)} transcripts)
</div>

<h2>Strand Distribution</h2>
<div class="strand-box strand-sense">
    <strong>Sense (+)</strong><br>
    {sense_count} hits ({100*sense_count/max(len(tr_hits),1):.1f}%)
</div>
<div class="strand-box strand-anti">
    <strong>Antisense (-)</strong><br>
    {anti_count} hits ({100*anti_count/max(len(tr_hits),1):.1f}%)
</div>

<h2>TE Class Distribution</h2>
<table>
<tr><th>TE Class</th><th>Hits</th><th>%</th></tr>
'''

    for te_class, count in sorted(class_counts.items(), key=lambda x: x[1], reverse=True):
        pct = 100 * count / max(len(tr_hits), 1)
        html += f'<tr><td><span class="class-tag class-{te_class}">{te_class}</span></td><td>{count}</td><td>{pct:.1f}%</td></tr>\n'

    html += '</table>\n'

    # Sequence visualization (first 2000 bp)
    html += '<h2>Sequence Annotation</h2>\n'
    html += '''<div class="legend">
<div class="legend-item"><div class="legend-color" style="background:#cce5ff"></div> TE sense (+)</div>
<div class="legend-item"><div class="legend-color" style="background:#fff3cd"></div> TE antisense (-)</div>
<div class="legend-item"><div class="legend-color" style="background:#d4edda"></div> Both strands</div>
</div>
'''

    html += '<div class="sequence-container">\n'

    # Show sequence in 80-char lines
    display_len = min(len(utr_seq), 2000)
    for line_start in range(0, display_len, 80):
        line_end = min(line_start + 80, display_len)
        html += f'<div class="seq-line"><span class="pos-num">{line_start+1:5d}</span>'

        for i in range(line_start, line_end):
            in_sense = i in sense_pos
            in_anti = i in anti_pos

            if in_sense and in_anti:
                html += f'<span class="te-both">{utr_seq[i]}</span>'
            elif in_sense:
                html += f'<span class="te-sense">{utr_seq[i]}</span>'
            elif in_anti:
                html += f'<span class="te-anti">{utr_seq[i]}</span>'
            else:
                html += f'<span class="normal">{utr_seq[i]}</span>'

        html += '</div>\n'

    if len(utr_seq) > 2000:
        html += f'<div class="seq-line"><em>... ({len(utr_seq) - 2000} more bp not shown)</em></div>\n'

    html += '</div>\n'

    # Top hits table
    html += '<h2>Top TE Hits</h2>\n'
    html += '<table>\n'
    html += '<tr><th>TE</th><th>Class</th><th>UTR pos</th><th>Length</th><th>Identity</th><th>E-value</th><th>Strand</th></tr>\n'

    sorted_hits = sorted(tr_hits, key=lambda x: x['bitscore'], reverse=True)[:50]
    for h in sorted_hits:
        te_id = h['te_id']
        te_class = te_info.get(te_id, {}).get('class', 'Unknown')
        strand_sym = '+' if h['strand'] == '+' else '-'
        html += f'<tr><td>{te_id}</td><td><span class="class-tag class-{te_class}">{te_class}</span></td>'
        html += f'<td>{h["qstart"]}-{h["qend"]}</td><td>{h["length"]}</td>'
        html += f'<td>{h["pident"]:.1f}%</td><td>{h["evalue"]:.2e}</td><td>{strand_sym}</td></tr>\n'

    html += '</table>\n'

    # Top alignments
    html += '<h2>Top Alignments</h2>\n'
    for h in sorted_hits[:10]:
        strand_class = 'sense' if h['strand'] == '+' else 'antisense'
        strand_sym = '+' if h['strand'] == '+' else '-'
        te_id = h['te_id']
        te_class = te_info.get(te_id, {}).get('class', 'Unknown')

        html += f'<div class="alignment-box {strand_class}">\n'
        html += f'<div class="align-header">{te_id} <span class="class-tag class-{te_class}">{te_class}</span> ({strand_sym})</div>\n'
        html += f'<div class="align-meta">UTR {h["qstart"]}-{h["qend"]} | TE {h["sstart"]}-{h["send"]} | {h["pident"]:.1f}% identity | {h["length"]}bp | E={h["evalue"]:.2e}</div>\n'

        if h['qseq'] and h['sseq']:
            html += '<div class="align-seq">'
            html += f'UTR: {h["qseq"][:100]}{"..." if len(h["qseq"]) > 100 else ""}\n'
            html += f'TE:  {h["sseq"][:100]}{"..." if len(h["sseq"]) > 100 else ""}'
            html += '</div>\n'

        html += '</div>\n'

    html += '</body>\n</html>'

    with open(output_path, 'w') as f:
        f.write(html)


def load_utr_sequences(fasta_path):
    """Load UTR sequences from FASTA file."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def analyze_strand_bias(transcript_hits, transcript_info, te_info, output_dir):
    """Analyze strand bias at UTR and TE family levels."""
    output_dir = Path(output_dir)

    # UTR-level strand bias
    utr_strand_data = []
    for fbtr, hits in transcript_hits.items():
        if fbtr not in transcript_info:
            continue

        sense_count = sum(1 for h in hits if h['strand'] == '+')
        anti_count = sum(1 for h in hits if h['strand'] == '-')
        total = sense_count + anti_count

        if total > 0:
            utr_strand_data.append({
                'fbtr': fbtr,
                'fbgn': transcript_info[fbtr]['fbgn'],
                'utr_len': transcript_info[fbtr]['length'],
                'total_hits': total,
                'sense_hits': sense_count,
                'anti_hits': anti_count,
                'sense_pct': 100 * sense_count / total,
                'anti_pct': 100 * anti_count / total
            })

    # Write UTR-level strand bias
    with open(output_dir / 'strand_bias_by_utr.tsv', 'w') as f:
        f.write('fbtr\tfbgn\tutr_len\ttotal_hits\tsense_hits\tanti_hits\tsense_pct\tanti_pct\n')
        for d in sorted(utr_strand_data, key=lambda x: x['total_hits'], reverse=True):
            f.write(f"{d['fbtr']}\t{d['fbgn']}\t{d['utr_len']}\t{d['total_hits']}\t"
                   f"{d['sense_hits']}\t{d['anti_hits']}\t{d['sense_pct']:.1f}\t{d['anti_pct']:.1f}\n")

    # TE family-level strand bias
    te_strand_data = defaultdict(lambda: {'sense': 0, 'anti': 0, 'total_bp': 0})

    for fbtr, hits in transcript_hits.items():
        for h in hits:
            te_id = h['te_id']
            # Get TE family name (strip instance number)
            # FBti IDs map to families - for now use te_id directly
            if h['strand'] == '+':
                te_strand_data[te_id]['sense'] += 1
            else:
                te_strand_data[te_id]['anti'] += 1
            te_strand_data[te_id]['total_bp'] += h['length']

    # Write TE-level strand bias
    te_strand_list = []
    for te_id, data in te_strand_data.items():
        total = data['sense'] + data['anti']
        if total >= 10:  # Minimum hits for meaningful stats
            te_strand_list.append({
                'te_id': te_id,
                'total_hits': total,
                'sense_hits': data['sense'],
                'anti_hits': data['anti'],
                'sense_pct': 100 * data['sense'] / total,
                'anti_pct': 100 * data['anti'] / total,
                'total_bp': data['total_bp']
            })

    with open(output_dir / 'strand_bias_by_te.tsv', 'w') as f:
        f.write('te_id\ttotal_hits\tsense_hits\tanti_hits\tsense_pct\tanti_pct\ttotal_bp\n')
        for d in sorted(te_strand_list, key=lambda x: x['total_hits'], reverse=True):
            f.write(f"{d['te_id']}\t{d['total_hits']}\t{d['sense_hits']}\t{d['anti_hits']}\t"
                   f"{d['sense_pct']:.1f}\t{d['anti_pct']:.1f}\t{d['total_bp']}\n")

    # Calculate overall stats
    total_sense = sum(d['sense_hits'] for d in utr_strand_data)
    total_anti = sum(d['anti_hits'] for d in utr_strand_data)
    total_all = total_sense + total_anti

    # Find extreme biased UTRs and TEs
    strongly_sense_utrs = [d for d in utr_strand_data if d['total_hits'] >= 50 and d['sense_pct'] >= 70]
    strongly_anti_utrs = [d for d in utr_strand_data if d['total_hits'] >= 50 and d['anti_pct'] >= 70]
    strongly_sense_tes = [d for d in te_strand_list if d['sense_pct'] >= 60]
    strongly_anti_tes = [d for d in te_strand_list if d['anti_pct'] >= 60]

    return {
        'total_sense': total_sense,
        'total_anti': total_anti,
        'sense_pct': 100 * total_sense / total_all if total_all > 0 else 0,
        'anti_pct': 100 * total_anti / total_all if total_all > 0 else 0,
        'utrs_analyzed': len(utr_strand_data),
        'tes_analyzed': len(te_strand_list),
        'strongly_sense_utrs': len(strongly_sense_utrs),
        'strongly_anti_utrs': len(strongly_anti_utrs),
        'strongly_sense_tes': len(strongly_sense_tes),
        'strongly_anti_tes': len(strongly_anti_tes),
        'utr_strand_data': utr_strand_data,
        'te_strand_list': te_strand_list
    }


def main():
    parser = argparse.ArgumentParser(description='Genome-wide TE analysis')
    parser.add_argument('--blast', default='results/genome_wide_all_3utrs.tsv',
                       help='BLAST results file')
    parser.add_argument('--utr-fasta', default='data/references/dmel_3utr.fasta',
                       help='3\'UTR FASTA file')
    parser.add_argument('--te-fasta', default='data/references/dmel_te_flybase.fasta',
                       help='TE database FASTA')
    parser.add_argument('--output', default='results',
                       help='Output directory')
    parser.add_argument('--html-top', type=int, default=10,
                       help='Number of top genes to visualize')
    parser.add_argument('--html-bottom', type=int, default=10,
                       help='Number of bottom genes to visualize')

    args = parser.parse_args()

    # Parse reference files
    print("Loading reference data...")
    transcript_info = parse_3utr_fasta(args.utr_fasta)
    print(f"  Loaded {len(transcript_info):,} transcripts")

    te_info = parse_te_database(args.te_fasta)
    print(f"  Loaded {len(te_info):,} TEs")

    # Parse BLAST results
    transcript_hits = parse_blast_results(args.blast, transcript_info)

    # Aggregate by gene
    print("\nAggregating by gene...")
    gene_data = aggregate_by_gene(transcript_hits, transcript_info)
    print(f"  {len(gene_data):,} genes with hits")

    # Write top/bottom 100
    print("\nWriting top/bottom 100 genes...")
    top_100, bottom_100, sorted_genes = write_top_bottom_100(gene_data, args.output)

    # Load UTR sequences for HTML visualization
    print("\nLoading UTR sequences for visualization...")
    utr_sequences = load_utr_sequences(args.utr_fasta)

    # Create HTML visualizations
    output_dir = Path(args.output) / 'te_annotations_genomewide'
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nGenerating HTML for top {args.html_top} genes...")
    for i, (fbgn, data) in enumerate(top_100[:args.html_top], 1):
        html_path = output_dir / f'top_{i:02d}_{fbgn}_te_annotation.html'
        generate_html_visualization(fbgn, gene_data, transcript_info, utr_sequences, te_info, html_path)
        print(f"  Created {html_path.name}")

    print(f"\nGenerating HTML for bottom {args.html_bottom} genes...")
    for i, (fbgn, data) in enumerate(bottom_100[:args.html_bottom], 1):
        html_path = output_dir / f'bottom_{i:02d}_{fbgn}_te_annotation.html'
        generate_html_visualization(fbgn, gene_data, transcript_info, utr_sequences, te_info, html_path)
        print(f"  Created {html_path.name}")

    # Create index page
    print("\nCreating index page...")
    with open(output_dir / 'index.html', 'w') as f:
        f.write('''<!DOCTYPE html>
<html>
<head>
<title>Genome-wide TE Analysis - Top/Bottom Genes</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, sans-serif; margin: 40px; }
h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
table { border-collapse: collapse; margin: 20px 0; }
th, td { padding: 8px 16px; border: 1px solid #ddd; text-align: left; }
th { background: #f5f5f5; }
a { color: #0066cc; }
.high { background: #ffe6e6; }
.low { background: #e6ffe6; }
</style>
</head>
<body>
<h1>Genome-wide 3'UTR TE Content Analysis</h1>
<h2>Top 10 Genes (Highest TE Density)</h2>
<table>
<tr><th>Rank</th><th>Gene (FBgn)</th><th>Density (bp/kb)</th><th>Hits</th><th>UTR Length</th><th>View</th></tr>
''')
        for i, (fbgn, data) in enumerate(top_100[:args.html_top], 1):
            f.write(f'<tr class="high"><td>{i}</td><td>{fbgn}</td><td>{data["density"]:.1f}</td>'
                   f'<td>{data["hit_count"]:,}</td><td>{data["total_utr_len"]:,}</td>'
                   f'<td><a href="top_{i:02d}_{fbgn}_te_annotation.html">View</a></td></tr>\n')

        f.write('''</table>
<h2>Bottom 10 Genes (Lowest TE Density)</h2>
<table>
<tr><th>Rank</th><th>Gene (FBgn)</th><th>Density (bp/kb)</th><th>Hits</th><th>UTR Length</th><th>View</th></tr>
''')
        for i, (fbgn, data) in enumerate(bottom_100[:args.html_bottom], 1):
            rank = len(sorted_genes) - 99 + (i - 1)
            f.write(f'<tr class="low"><td>{rank}</td><td>{fbgn}</td><td>{data["density"]:.1f}</td>'
                   f'<td>{data["hit_count"]:,}</td><td>{data["total_utr_len"]:,}</td>'
                   f'<td><a href="bottom_{i:02d}_{fbgn}_te_annotation.html">View</a></td></tr>\n')

        f.write('</table>\n</body>\n</html>')

    # Strand bias analysis
    print("\nAnalyzing strand bias genome-wide...")
    strand_stats = analyze_strand_bias(transcript_hits, transcript_info, te_info, args.output)

    print(f"\n=== Strand Bias Summary ===")
    print(f"Total hits: {strand_stats['total_sense'] + strand_stats['total_anti']:,}")
    print(f"Sense (+): {strand_stats['total_sense']:,} ({strand_stats['sense_pct']:.1f}%)")
    print(f"Antisense (-): {strand_stats['total_anti']:,} ({strand_stats['anti_pct']:.1f}%)")
    print(f"\nUTRs analyzed: {strand_stats['utrs_analyzed']:,}")
    print(f"  Strongly sense-biased (≥70%): {strand_stats['strongly_sense_utrs']:,}")
    print(f"  Strongly antisense-biased (≥70%): {strand_stats['strongly_anti_utrs']:,}")
    print(f"\nTEs analyzed (≥10 hits): {strand_stats['tes_analyzed']:,}")
    print(f"  Sense-biased (≥60%): {strand_stats['strongly_sense_tes']:,}")
    print(f"  Antisense-biased (≥60%): {strand_stats['strongly_anti_tes']:,}")

    print(f"\nDone! Results in {args.output}/")


if __name__ == '__main__':
    main()
