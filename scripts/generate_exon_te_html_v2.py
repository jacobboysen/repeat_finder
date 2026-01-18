#!/usr/bin/env python3
"""
Generate HTML visualizations of TE BLAST hits for top exon genes (v2).

Improved version with:
- Distinct colors for different annotation types
- Start/stop codon markers
- Proper transcript handling
- Clear visual separation of annotation layers
"""

import argparse
import sys
import re
from collections import defaultdict
from pathlib import Path
from html import escape

sys.path.insert(0, str(Path(__file__).parent))
from utils.te_domain_classifier import classify_te_domain, infer_te_class, get_domain_description

# DISTINCT domain colors - use very different hues
DOMAIN_COLORS = {
    # Coding domains - warm colors
    'gag': '#DC143C',         # Crimson red
    'pol': '#8B008B',         # Dark magenta
    'orf1': '#FF4500',        # Orange red
    'orf2_rt': '#FF6347',     # Tomato
    'orf2_en': '#FF8C00',     # Dark orange
    'transposase': '#FFD700', # Gold
    # Non-coding - cool colors
    '5_ltr': '#1E90FF',       # Dodger blue
    '3_ltr': '#00CED1',       # Dark turquoise
    '5_tir': '#32CD32',       # Lime green
    '3_tir': '#228B22',       # Forest green
    '5_utr': '#87CEEB',       # Sky blue
    '3_utr': '#98FB98',       # Pale green
    'unknown': '#A9A9A9',     # Dark gray
}

# Structural annotation colors - completely different palette
STRUCT_COLORS = {
    'start_codon': '#00FF00',   # Bright green
    'stop_codon': '#FF0000',    # Bright red
    'exon_boundary': '#000000', # Black
}


def load_exon_metadata(path):
    """Load exon metadata."""
    exons = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            exons[row['exon_id']] = row
    return exons


def load_blast_hits_for_genes(blast_path, target_genes, exon_meta):
    """Load BLAST hits for specific genes."""
    hits_by_gene = defaultdict(list)
    exon_to_gene = {eid: meta['fbgn'] for eid, meta in exon_meta.items()}

    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 16:
                continue

            exon_id = parts[0]
            if exon_id not in exon_to_gene:
                continue

            fbgn = exon_to_gene[exon_id]
            if fbgn not in target_genes:
                continue

            sstart = int(parts[8])
            send = int(parts[9])
            slen = int(parts[13])
            te_id = parts[1]

            domain_info = classify_te_domain(te_id, sstart, send, slen)

            hits_by_gene[fbgn].append({
                'exon_id': exon_id,
                'te_id': te_id,
                'pident': float(parts[2]),
                'length': int(parts[3]),
                'qstart': int(parts[6]),
                'qend': int(parts[7]),
                'sstart': sstart,
                'send': send,
                'evalue': float(parts[10]),
                'bitscore': float(parts[11]),
                'qseq': parts[14],
                'sseq': parts[15],
                'strand': '+' if sstart < send else '-',
                'te_class': domain_info['te_class'],
                'te_domain': domain_info['domain'],
            })

    return hits_by_gene


def load_exon_sequences(fasta_path, target_exons):
    """Load exon sequences for target exons."""
    seqs = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_id and current_id in target_exons:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id and current_id in target_exons:
            seqs[current_id] = ''.join(current_seq)

    return seqs


def find_codons(seq):
    """Find start (ATG) and stop codons in sequence."""
    start_codons = []
    stop_codons = []

    seq_upper = seq.upper()

    # Find ATG (start)
    for i in range(len(seq_upper) - 2):
        if seq_upper[i:i+3] == 'ATG':
            start_codons.append(i)

    # Find stop codons (TAA, TAG, TGA)
    for i in range(len(seq_upper) - 2):
        codon = seq_upper[i:i+3]
        if codon in ('TAA', 'TAG', 'TGA'):
            stop_codons.append((i, codon))

    return start_codons, stop_codons


def generate_gene_html(fbgn, symbol, hits, exon_meta, exon_seqs, output_dir):
    """Generate HTML for a single gene."""

    # Get all exons/transcripts for this gene
    gene_exons = {eid: meta for eid, meta in exon_meta.items() if meta['fbgn'] == fbgn}

    # Group by transcript
    transcripts = defaultdict(list)
    for eid, meta in gene_exons.items():
        transcripts[meta['fbtr']].append((eid, meta))

    # Sort exons within each transcript
    for fbtr in transcripts:
        transcripts[fbtr].sort(key=lambda x: int(x[1]['exon_num']))

    # Stats
    total_hits = len(hits)
    sense_hits = sum(1 for h in hits if h['strand'] == '+')

    domain_counts = defaultdict(int)
    te_family_counts = defaultdict(int)
    for h in hits:
        domain_counts[h['te_domain']] += 1
        te_family_counts[h['te_id']] += 1

    top_domains = sorted(domain_counts.items(), key=lambda x: -x[1])[:8]
    top_families = sorted(te_family_counts.items(), key=lambda x: -x[1])[:10]

    # Start HTML
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{symbol} ({fbgn}) - Exon TE Analysis v2</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 14px;
    line-height: 1.6;
    margin: 40px;
    background: #fafafa;
    color: #333;
    max-width: 1600px;
}}
h1, h2, h3 {{ font-weight: 600; }}
h1 {{ border-bottom: 3px solid #333; padding-bottom: 10px; }}
h2 {{ border-bottom: 1px solid #ccc; padding-bottom: 5px; margin-top: 30px; }}

.summary-box {{
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: white;
    padding: 20px 25px;
    border-radius: 12px;
    margin: 20px 0;
    display: flex;
    flex-wrap: wrap;
    gap: 40px;
}}
.summary-stat .label {{ font-size: 11px; opacity: 0.9; text-transform: uppercase; letter-spacing: 1px; }}
.summary-stat .value {{ font-size: 28px; font-weight: 700; }}

.legend-section {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    margin: 20px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}}
.legend-title {{
    font-weight: 700;
    font-size: 16px;
    margin-bottom: 15px;
    padding-bottom: 10px;
    border-bottom: 2px solid #eee;
}}
.legend-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
    gap: 12px;
}}
.legend-item {{
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 8px 12px;
    background: #f8f9fa;
    border-radius: 6px;
    font-size: 13px;
}}
.legend-color {{
    width: 24px;
    height: 24px;
    border-radius: 4px;
    border: 2px solid rgba(0,0,0,0.2);
    flex-shrink: 0;
}}

.transcript-section {{
    background: white;
    padding: 25px;
    border-radius: 12px;
    box-shadow: 0 2px 12px rgba(0,0,0,0.08);
    margin: 25px 0;
}}
.transcript-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding-bottom: 15px;
    border-bottom: 2px solid #eee;
    margin-bottom: 20px;
}}
.transcript-title {{
    font-size: 20px;
    font-weight: 700;
    color: #333;
}}
.transcript-meta {{
    color: #666;
    font-size: 13px;
}}

.sequence-container {{
    background: #1e1e1e;
    padding: 20px;
    border-radius: 8px;
    overflow-x: auto;
    margin: 15px 0;
}}
.seq-line {{
    white-space: pre;
    font-family: "JetBrains Mono", "Fira Code", "SF Mono", monospace;
    font-size: 13px;
    line-height: 1.6;
    color: #d4d4d4;
}}
.pos-num {{
    color: #6a9955;
    margin-right: 15px;
    display: inline-block;
    width: 60px;
    user-select: none;
}}

/* Domain hit styling - background colors */
.hit-gag {{ background-color: #DC143C; color: white; }}
.hit-pol {{ background-color: #8B008B; color: white; }}
.hit-orf1 {{ background-color: #FF4500; color: white; }}
.hit-orf2_rt {{ background-color: #FF6347; color: white; }}
.hit-orf2_en {{ background-color: #FF8C00; color: black; }}
.hit-transposase {{ background-color: #FFD700; color: black; }}
.hit-5_ltr {{ background-color: #1E90FF; color: white; }}
.hit-3_ltr {{ background-color: #00CED1; color: black; }}
.hit-5_tir {{ background-color: #32CD32; color: black; }}
.hit-3_tir {{ background-color: #228B22; color: white; }}
.hit-unknown {{ background-color: #A9A9A9; color: black; }}

/* Structural annotations - distinct styling */
.codon-start {{
    background-color: #00FF00 !important;
    color: black !important;
    font-weight: bold;
    border: 2px solid #006400;
    border-radius: 3px;
    padding: 0 2px;
}}
.codon-stop {{
    background-color: #FF0000 !important;
    color: white !important;
    font-weight: bold;
    border: 2px solid #8B0000;
    border-radius: 3px;
    padding: 0 2px;
}}
.exon-start {{
    border-left: 4px solid #000000;
    margin-left: 2px;
    padding-left: 3px;
}}
.exon-marker {{
    display: inline-block;
    background: #000000;
    color: #00FF00;
    font-size: 10px;
    font-weight: bold;
    padding: 2px 6px;
    border-radius: 3px;
    margin: 0 3px;
    vertical-align: middle;
}}

/* Annotation track above sequence */
.annotation-track {{
    background: #2d2d2d;
    padding: 10px 20px;
    border-radius: 8px 8px 0 0;
    margin-bottom: 0;
    font-family: monospace;
    font-size: 11px;
    color: #aaa;
    overflow-x: auto;
    white-space: pre;
}}
.track-label {{
    display: inline-block;
    width: 80px;
    color: #888;
}}

table {{ border-collapse: collapse; margin: 15px 0; width: 100%; }}
th, td {{ padding: 10px 14px; border: 1px solid #ddd; text-align: left; }}
th {{ background: #f0f0f0; font-weight: 600; }}
tr:nth-child(even) {{ background: #fafafa; }}
tr:hover {{ background: #f5f5f5; }}

.alignment-box {{
    background: #f8f9fa;
    padding: 15px;
    margin: 12px 0;
    border-radius: 8px;
    border-left: 5px solid;
}}
.alignment-box.strand-plus {{ border-color: #28a745; }}
.alignment-box.strand-minus {{ border-color: #ffc107; }}
.align-header {{
    font-weight: 700;
    margin-bottom: 8px;
    display: flex;
    align-items: center;
    gap: 10px;
}}
.align-meta {{ color: #666; font-size: 12px; margin-bottom: 10px; }}
.align-seq {{
    white-space: pre;
    font-family: "JetBrains Mono", monospace;
    font-size: 12px;
    line-height: 1.5;
    background: #2d2d2d;
    color: #d4d4d4;
    padding: 12px;
    border-radius: 4px;
}}
.match {{ color: #4ec9b0; }}
.mismatch {{ color: #f14c4c; }}
.gap {{ color: #666; }}

.domain-badge {{
    display: inline-block;
    padding: 4px 10px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}}
.strand-badge {{
    display: inline-block;
    padding: 3px 8px;
    border-radius: 50%;
    font-size: 14px;
    font-weight: bold;
    width: 24px;
    height: 24px;
    text-align: center;
    line-height: 18px;
}}
.strand-plus {{ background: #28a745; color: white; }}
.strand-minus {{ background: #ffc107; color: black; }}

.nav-links {{
    margin: 20px 0;
    padding: 12px 16px;
    background: white;
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}}
.nav-links a {{
    color: #0066cc;
    text-decoration: none;
    font-weight: 500;
}}
.nav-links a:hover {{ text-decoration: underline; }}

.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 15px;
    margin: 20px 0;
}}
.stat-card {{
    background: white;
    padding: 18px;
    border-radius: 10px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
}}
.stat-card h4 {{
    margin: 0 0 12px 0;
    font-size: 13px;
    color: #666;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}}
.stat-value {{ font-size: 24px; font-weight: 700; }}
</style>
</head>
<body>

<div class="nav-links">
<a href="index.html">← Back to Index</a>
</div>

<h1>{escape(symbol)} <span style="font-weight: 400; color: #666;">({fbgn})</span></h1>
<p style="color: #666; font-size: 15px;">Exon-level TE similarity analysis with domain and structural annotations</p>

<div class="summary-box">
    <div class="summary-stat">
        <div class="label">Transcripts</div>
        <div class="value">{len(transcripts)}</div>
    </div>
    <div class="summary-stat">
        <div class="label">Total Exons</div>
        <div class="value">{len(gene_exons)}</div>
    </div>
    <div class="summary-stat">
        <div class="label">TE Hits</div>
        <div class="value">{total_hits:,}</div>
    </div>
    <div class="summary-stat">
        <div class="label">Sense Strand</div>
        <div class="value">{100*sense_hits/max(total_hits,1):.1f}%</div>
    </div>
    <div class="summary-stat">
        <div class="label">TE Families</div>
        <div class="value">{len(te_family_counts)}</div>
    </div>
</div>
'''

    # LEGEND SECTION - Two separate legends
    html += '''
<div class="legend-section">
<div class="legend-title">🎨 TE Domain Colors (hit regions)</div>
<div class="legend-grid">
'''

    # Only show domains actually used
    domains_used = set(h['te_domain'] for h in hits)
    for domain in ['gag', 'pol', 'orf1', 'orf2_rt', 'orf2_en', 'transposase', '5_ltr', '3_ltr', '5_tir', '3_tir', 'unknown']:
        if domain in domains_used:
            color = DOMAIN_COLORS.get(domain, '#A9A9A9')
            desc = get_domain_description(domain)
            html += f'''<div class="legend-item">
    <div class="legend-color" style="background:{color}"></div>
    <div><strong>{domain}</strong><br><small style="color:#666">{desc}</small></div>
</div>
'''

    html += '''</div>
</div>

<div class="legend-section">
<div class="legend-title">🔬 Structural Annotations</div>
<div class="legend-grid">
    <div class="legend-item">
        <div class="legend-color" style="background:#00FF00; border-color:#006400"></div>
        <div><strong>ATG</strong><br><small style="color:#666">Start codon</small></div>
    </div>
    <div class="legend-item">
        <div class="legend-color" style="background:#FF0000; border-color:#8B0000"></div>
        <div><strong>TAA/TAG/TGA</strong><br><small style="color:#666">Stop codon</small></div>
    </div>
    <div class="legend-item">
        <div class="legend-color" style="background:#000000"></div>
        <div><strong>|E#|</strong><br><small style="color:#666">Exon boundary</small></div>
    </div>
    <div class="legend-item">
        <div class="legend-color" style="background:#28a745"></div>
        <div><strong>+ Sense</strong><br><small style="color:#666">Same as gene strand</small></div>
    </div>
    <div class="legend-item">
        <div class="legend-color" style="background:#ffc107"></div>
        <div><strong>- Antisense</strong><br><small style="color:#666">Opposite strand</small></div>
    </div>
</div>
</div>
'''

    # Domain distribution stats
    html += '<h2>TE Domain Distribution</h2>\n<div class="stats-grid">\n'

    for domain, count in top_domains:
        pct = 100 * count / total_hits if total_hits else 0
        color = DOMAIN_COLORS.get(domain, '#A9A9A9')
        desc = get_domain_description(domain)
        html += f'''<div class="stat-card">
    <h4><span class="domain-badge" style="background:{color}; color:{'white' if domain in ['gag','pol','5_ltr','3_tir'] else 'black'}">{domain}</span></h4>
    <div class="stat-value">{count:,} <span style="font-size:14px;color:#666">({pct:.1f}%)</span></div>
    <div style="color:#888; font-size:12px; margin-top:5px;">{desc}</div>
</div>
'''

    html += '</div>\n'

    # Top TE families
    html += '''<h2>Top TE Families</h2>
<table>
<tr><th>TE Family</th><th>Hits</th><th>Class</th><th>Percentage</th></tr>
'''
    for te_id, count in top_families:
        te_class = infer_te_class(te_id)
        pct = 100 * count / total_hits if total_hits else 0
        html += f'<tr><td><strong>{escape(te_id)}</strong></td><td>{count:,}</td><td>{te_class}</td><td>{pct:.1f}%</td></tr>\n'

    html += '</table>\n'

    # Group hits by exon
    hits_by_exon = defaultdict(list)
    for h in hits:
        hits_by_exon[h['exon_id']].append(h)

    # TRANSCRIPT SECTIONS
    html += '<h2>Transcript Details</h2>\n'

    for fbtr, exon_list in sorted(transcripts.items()):
        # Concatenate all exons for this transcript
        full_seq = ''
        exon_boundaries = []
        exon_info_list = []

        for exon_id, meta in exon_list:
            seq = exon_seqs.get(exon_id, '')
            if full_seq:  # Not first exon
                exon_boundaries.append((len(full_seq), meta['exon_num']))
            exon_info_list.append({
                'exon_id': exon_id,
                'exon_num': meta['exon_num'],
                'offset': len(full_seq),
                'length': len(seq),
                'position': meta['position'],
                'utr_overlap': meta['utr_overlap'],
            })
            full_seq += seq

        if not full_seq:
            continue

        # Find codons
        start_codons, stop_codons = find_codons(full_seq)

        # Build domain coverage map
        position_domains = {}
        for exon_id, meta in exon_list:
            exon_hits = hits_by_exon.get(exon_id, [])
            exon_offset = next((e['offset'] for e in exon_info_list if e['exon_id'] == exon_id), 0)

            for h in exon_hits:
                domain = h['te_domain']
                for p in range(h['qstart'] - 1, h['qend']):
                    global_pos = exon_offset + p
                    if global_pos < len(full_seq):
                        if global_pos not in position_domains:
                            position_domains[global_pos] = domain

        # Count hits for this transcript
        transcript_hits = sum(len(hits_by_exon.get(e['exon_id'], [])) for e in exon_info_list)

        html += f'''
<div class="transcript-section">
<div class="transcript-header">
    <div>
        <span class="transcript-title">{fbtr}</span>
        <span style="margin-left: 15px; color: #666;">({len(exon_list)} exon{"s" if len(exon_list) > 1 else ""})</span>
    </div>
    <div class="transcript-meta">
        {len(full_seq):,} bp | {transcript_hits:,} TE hits |
        {len(start_codons)} ATG | {len(stop_codons)} stop codons
    </div>
</div>
'''

        # Exon structure diagram
        if len(exon_list) > 1:
            html += '<div style="margin-bottom: 15px;">\n'
            html += '<strong style="font-size: 12px; color: #666;">Exon Structure:</strong> '
            for i, info in enumerate(exon_info_list):
                colors = {'first_exon': '#28a745', 'internal_exon': '#007bff', 'last_exon': '#ffc107', 'single_exon': '#dc3545'}
                color = colors.get(info['position'], '#6c757d')
                html += f'<span style="display:inline-block; background:{color}; color:white; padding:3px 10px; border-radius:4px; margin:2px; font-size:12px;">E{info["exon_num"]} ({info["length"]}bp)</span>'
                if i < len(exon_info_list) - 1:
                    html += '<span style="color:#999; margin:0 5px;">→</span>'
            html += '\n</div>\n'

        # Render sequence with all annotations
        html += '<div class="sequence-container">\n'

        display_len = min(len(full_seq), 4000)
        boundary_positions = {b[0]: b[1] for b in exon_boundaries}
        start_codon_set = set(start_codons)
        stop_codon_dict = {pos: codon for pos, codon in stop_codons}

        for line_start in range(0, display_len, 80):
            line_end = min(line_start + 80, display_len)
            html += f'<div class="seq-line"><span class="pos-num">{line_start+1:5d}</span>'

            i = line_start
            while i < line_end:
                char = full_seq[i]

                # Check for exon boundary
                if i in boundary_positions:
                    html += f'<span class="exon-marker">E{boundary_positions[i]}</span>'

                # Determine styling
                classes = []
                styles = []

                # Check for start codon (ATG)
                if i in start_codon_set and i + 2 < len(full_seq):
                    html += f'<span class="codon-start">{full_seq[i:i+3]}</span>'
                    i += 3
                    continue

                # Check for stop codon
                if i in stop_codon_dict and i + 2 < len(full_seq):
                    html += f'<span class="codon-stop">{full_seq[i:i+3]}</span>'
                    i += 3
                    continue

                # Check for domain hit
                domain = position_domains.get(i)
                if domain:
                    color = DOMAIN_COLORS.get(domain, '#A9A9A9')
                    text_color = 'white' if domain in ['gag', 'pol', '5_ltr', '3_ltr', '3_tir', 'orf1'] else 'black'
                    html += f'<span style="background-color:{color}; color:{text_color};">{char}</span>'
                else:
                    html += char

                i += 1

            html += '</div>\n'

        if len(full_seq) > 4000:
            html += f'<div class="seq-line" style="color:#888;"><em>... ({len(full_seq) - 4000:,} more bp not shown)</em></div>\n'

        html += '</div>\n'  # Close sequence-container

        # Top alignments for this transcript
        transcript_all_hits = []
        for info in exon_info_list:
            transcript_all_hits.extend(hits_by_exon.get(info['exon_id'], []))

        if transcript_all_hits:
            html += '<h4 style="margin-top: 20px;">Top Alignments</h4>\n'
            sorted_hits = sorted(transcript_all_hits, key=lambda x: -x['bitscore'])[:15]

            for h in sorted_hits:
                strand_class = 'strand-plus' if h['strand'] == '+' else 'strand-minus'
                domain = h['te_domain']
                domain_color = DOMAIN_COLORS.get(domain, '#A9A9A9')

                # Color alignment
                qseq = h['qseq'][:70]
                sseq = h['sseq'][:70]
                colored_match = ''
                for q, s in zip(qseq, sseq):
                    if q == s:
                        colored_match += '<span class="match">|</span>'
                    elif q == '-' or s == '-':
                        colored_match += '<span class="gap"> </span>'
                    else:
                        colored_match += '<span class="mismatch">.</span>'

                html += f'''<div class="alignment-box {strand_class}">
<div class="align-header">
    <span class="strand-badge {strand_class}">{h['strand']}</span>
    <strong>{escape(h['te_id'])}</strong>
    <span class="domain-badge" style="background:{domain_color}; color:{'white' if domain in ['gag','pol','5_ltr','3_ltr','3_tir','orf1'] else 'black'}">{domain}</span>
</div>
<div class="align-meta">
    Position {h['qstart']}-{h['qend']} in {h['exon_id']} |
    {h['pident']:.1f}% identity | {h['length']}bp |
    E-value: {h['evalue']:.2e} | Bitscore: {h['bitscore']:.1f}
</div>
<div class="align-seq">Query: {escape(qseq)}{"..." if len(h['qseq']) > 70 else ""}
       {colored_match}
   TE: {escape(sseq)}{"..." if len(h['sseq']) > 70 else ""}</div>
</div>
'''

        html += '</div>\n'  # Close transcript-section

    html += '</body>\n</html>'

    # Write file
    safe_symbol = symbol.replace('/', '_').replace(' ', '_')
    html_path = output_dir / f'{safe_symbol}_{fbgn}.html'
    with open(html_path, 'w') as f:
        f.write(html)

    return html_path.name


def generate_index_html(gene_info, output_dir):
    """Generate index page."""
    html = '''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Top 50 Genes - Exon TE Analysis v2</title>
<style>
body {
    font-family: -apple-system, sans-serif;
    margin: 40px;
    max-width: 1400px;
    background: #fafafa;
}
h1 { border-bottom: 3px solid #333; padding-bottom: 10px; }
table {
    border-collapse: collapse;
    margin: 20px 0;
    width: 100%;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    border-radius: 8px;
    overflow: hidden;
}
th, td { padding: 12px 16px; border-bottom: 1px solid #eee; text-align: left; }
th { background: #f5f5f5; font-weight: 600; position: sticky; top: 0; }
tr:hover { background: #f8f9fa; }
a { color: #0066cc; text-decoration: none; font-weight: 500; }
a:hover { text-decoration: underline; }
.domain-badge {
    display: inline-block;
    padding: 3px 8px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: 600;
}
</style>
</head>
<body>
<h1>Top 50 Genes by Exon TE Density (v2)</h1>
<p>Enhanced visualization with distinct domain colors, start/stop codons, and exon boundaries.</p>
<table>
<tr>
<th>Rank</th>
<th>Gene</th>
<th>Transcripts</th>
<th>Exons</th>
<th>TE Hits</th>
<th>Density</th>
<th>Sense %</th>
<th>Top Domain</th>
<th>View</th>
</tr>
'''

    for info in gene_info:
        domain = info['top_domain']
        color = DOMAIN_COLORS.get(domain, '#A9A9A9')
        text_color = 'white' if domain in ['gag','pol','5_ltr','3_ltr','3_tir','orf1'] else 'black'

        html += f'''<tr>
<td>{info['rank']}</td>
<td><strong>{escape(info['symbol'])}</strong><br><small style="color:#666">{info['fbgn']}</small></td>
<td>{info['transcripts']}</td>
<td>{info['exons']}</td>
<td>{info['hits']:,}</td>
<td>{info['density']:.1f}</td>
<td>{info['sense_pct']:.1f}%</td>
<td><span class="domain-badge" style="background:{color}; color:{text_color}">{domain}</span></td>
<td><a href="{info['html_file']}">View →</a></td>
</tr>
'''

    html += '</table>\n</body>\n</html>'

    with open(output_dir / 'index.html', 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--top', type=int, default=50)
    parser.add_argument('--output-dir', type=Path, default=Path('results/exon_analysis/html_v2'))
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    blast_path = Path('results/exon_analysis/genome_wide_exons.tsv')
    gene_stats_path = Path('results/exon_analysis/gene_exon_te_summary.tsv')
    exon_meta_path = Path('data/queries/genome_wide/exon_metadata.tsv')
    exon_fasta_path = Path('data/queries/genome_wide/exons_sense.fasta')

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading top {args.top} genes...")
    top_genes = {}
    with open(gene_stats_path) as f:
        f.readline()
        for i, line in enumerate(f):
            if i >= args.top:
                break
            parts = line.strip().split('\t')
            fbgn = parts[1]
            top_genes[fbgn] = {
                'rank': int(parts[0]),
                'symbol': parts[2],
                'exons': int(parts[3]),
                'hits': int(parts[5]),
                'density': float(parts[8]),
                'sense_pct': float(parts[11]) if len(parts) > 11 else 50,
            }

    print(f"  Found {len(top_genes)} genes")

    print("Loading exon metadata...")
    exon_meta = load_exon_metadata(exon_meta_path)

    target_exon_ids = {eid for eid, meta in exon_meta.items() if meta['fbgn'] in top_genes}
    print(f"  {len(target_exon_ids)} exons for target genes")

    print("Loading exon sequences...")
    exon_seqs = load_exon_sequences(exon_fasta_path, target_exon_ids)

    print("Loading BLAST hits...")
    hits_by_gene = load_blast_hits_for_genes(blast_path, top_genes, exon_meta)

    print(f"\nGenerating HTML pages...")
    gene_info = []

    for fbgn, info in sorted(top_genes.items(), key=lambda x: x[1]['rank']):
        symbol = info['symbol']
        hits = hits_by_gene.get(fbgn, [])

        # Count transcripts
        gene_exons = {eid: meta for eid, meta in exon_meta.items() if meta['fbgn'] == fbgn}
        transcripts = set(meta['fbtr'] for meta in gene_exons.values())

        domain_counts = defaultdict(int)
        for h in hits:
            domain_counts[h['te_domain']] += 1
        top_domain = max(domain_counts, key=domain_counts.get) if domain_counts else 'unknown'

        html_file = generate_gene_html(fbgn, symbol, hits, exon_meta, exon_seqs, args.output_dir)

        gene_info.append({
            'rank': info['rank'],
            'fbgn': fbgn,
            'symbol': symbol,
            'transcripts': len(transcripts),
            'exons': len(gene_exons),
            'hits': len(hits),
            'density': info['density'],
            'sense_pct': info['sense_pct'],
            'top_domain': top_domain,
            'html_file': html_file,
        })

        if args.verbose:
            print(f"  {info['rank']:3d}. {symbol}: {len(hits):,} hits, {len(transcripts)} transcripts -> {html_file}")

    print("Generating index page...")
    generate_index_html(gene_info, args.output_dir)

    print(f"\nDone! Open {args.output_dir / 'index.html'} to view results.")


if __name__ == '__main__':
    main()
