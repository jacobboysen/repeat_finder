#!/usr/bin/env python3
"""
Generate HTML visualizations of TE BLAST hits for top exon genes.

Creates one HTML page per gene showing:
- Full gene view with exon boundaries marked
- Sequence colored by TE domain (gag/pol/LTR/etc)
- Per-exon alignments with domain annotations
- Top TE families hitting each exon
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from html import escape

sys.path.insert(0, str(Path(__file__).parent))
from utils.te_domain_classifier import classify_te_domain, infer_te_class, get_domain_description

# Domain colors for visualization
DOMAIN_COLORS = {
    'gag': '#e74c3c',        # Red
    'pol': '#9b59b6',        # Purple
    '5_ltr': '#3498db',      # Blue
    '3_ltr': '#2980b9',      # Darker blue
    'orf1': '#1abc9c',       # Teal
    'orf2_rt': '#16a085',    # Darker teal
    'orf2_en': '#27ae60',    # Green
    'transposase': '#e67e22', # Orange
    '5_tir': '#f39c12',      # Yellow-orange
    '3_tir': '#d68910',      # Darker orange
    '5_utr': '#5dade2',      # Light blue
    '3_utr': '#85c1e9',      # Lighter blue
    'unknown': '#bdc3c7',    # Gray
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

    # Build exon_id -> fbgn mapping
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

            # Classify TE domain
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


def generate_gene_html(fbgn, symbol, hits, exon_meta, exon_seqs, output_dir):
    """Generate HTML for a single gene."""

    # Get exons for this gene
    gene_exons = {eid: meta for eid, meta in exon_meta.items()
                  if meta['fbgn'] == fbgn}

    # Sort exons by exon number
    sorted_exons = sorted(gene_exons.items(),
                          key=lambda x: int(x[1]['exon_num']))

    # Group hits by exon
    hits_by_exon = defaultdict(list)
    for h in hits:
        hits_by_exon[h['exon_id']].append(h)

    # Count stats
    total_hits = len(hits)
    sense_hits = sum(1 for h in hits if h['strand'] == '+')

    # TE domain distribution
    domain_counts = defaultdict(int)
    te_class_counts = defaultdict(int)
    te_family_counts = defaultdict(int)
    for h in hits:
        domain_counts[h['te_domain']] += 1
        te_class_counts[h['te_class']] += 1
        te_family_counts[h['te_id']] += 1

    top_domains = sorted(domain_counts.items(), key=lambda x: -x[1])[:5]
    top_families = sorted(te_family_counts.items(), key=lambda x: -x[1])[:10]

    # Start HTML
    html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{symbol} ({fbgn}) - Exon TE Analysis</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 14px;
    line-height: 1.6;
    margin: 40px;
    background: #fafafa;
    color: #333;
    max-width: 1400px;
}}
h1, h2, h3 {{ font-weight: 600; }}
h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
.summary-box {{
    background: #e8f4e8;
    padding: 15px 20px;
    border-radius: 8px;
    margin: 20px 0;
    display: flex;
    flex-wrap: wrap;
    gap: 30px;
}}
.summary-stat {{
    min-width: 150px;
}}
.summary-stat .label {{ font-size: 12px; color: #666; text-transform: uppercase; }}
.summary-stat .value {{ font-size: 24px; font-weight: 600; }}
.exon-container {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin: 20px 0;
}}
.exon-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    border-bottom: 1px solid #eee;
    padding-bottom: 10px;
    margin-bottom: 15px;
}}
.exon-title {{ font-size: 18px; font-weight: 600; }}
.exon-meta {{ color: #666; font-size: 13px; }}
.position-tag {{
    display: inline-block;
    padding: 3px 10px;
    border-radius: 12px;
    font-size: 12px;
    font-weight: 600;
    margin-left: 10px;
}}
.pos-first {{ background: #d4edda; color: #155724; }}
.pos-internal {{ background: #cce5ff; color: #004085; }}
.pos-last {{ background: #fff3cd; color: #856404; }}
.pos-single {{ background: #f8d7da; color: #721c24; }}
.utr-tag {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 3px;
    font-size: 11px;
    margin-left: 5px;
}}
.utr-5 {{ background: #28a745; color: white; }}
.utr-3 {{ background: #ffc107; color: #333; }}
.utr-both {{ background: #dc3545; color: white; }}
.utr-none {{ background: #6c757d; color: white; }}
.sequence-box {{
    background: #f8f9fa;
    padding: 15px;
    border-radius: 4px;
    overflow-x: auto;
    margin: 10px 0;
}}
.seq-line {{
    white-space: pre;
    font-family: "SF Mono", "Monaco", "Inconsolata", "Fira Mono", monospace;
    font-size: 12px;
    line-height: 1.4;
}}
.pos-num {{ color: #888; margin-right: 10px; display: inline-block; width: 50px; }}
.te-sense {{ background-color: #cce5ff; }}
.te-anti {{ background-color: #fff3cd; }}
.te-both {{ background-color: #d4edda; }}
.exon-boundary {{
    border-left: 3px solid #dc3545;
    padding-left: 2px;
    margin-left: 1px;
}}
.exon-boundary-marker {{
    display: inline-block;
    background: #dc3545;
    color: white;
    font-size: 9px;
    padding: 1px 4px;
    border-radius: 2px;
    vertical-align: super;
    margin: 0 2px;
}}
.full-gene-view {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin: 20px 0;
}}
.domain-legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin: 15px 0;
    padding: 10px;
    background: #f8f9fa;
    border-radius: 4px;
}}
.domain-legend-item {{
    display: flex;
    align-items: center;
    gap: 5px;
    font-size: 12px;
}}
.domain-color-box {{
    width: 14px;
    height: 14px;
    border-radius: 2px;
    border: 1px solid rgba(0,0,0,0.2);
}}
.view-toggle {{
    margin: 10px 0;
}}
.view-toggle button {{
    padding: 8px 16px;
    margin-right: 8px;
    border: 1px solid #ddd;
    background: white;
    border-radius: 4px;
    cursor: pointer;
}}
.view-toggle button.active {{
    background: #007bff;
    color: white;
    border-color: #007bff;
}}
.gene-schematic {{
    margin: 15px 0;
    padding: 15px;
    background: #f8f9fa;
    border-radius: 4px;
}}
.schematic-bar {{
    height: 30px;
    display: flex;
    margin: 10px 0;
    border-radius: 4px;
    overflow: hidden;
}}
.schematic-exon {{
    display: flex;
    align-items: center;
    justify-content: center;
    color: white;
    font-size: 11px;
    font-weight: 600;
    min-width: 30px;
    border-right: 2px solid white;
}}
.schematic-intron {{
    background: #eee;
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 10px;
    color: #666;
    min-width: 20px;
}}
.domain-tag {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 3px;
    font-size: 11px;
    font-weight: 600;
    margin-right: 4px;
}}
.domain-gag {{ background: #e74c3c; color: white; }}
.domain-pol {{ background: #9b59b6; color: white; }}
.domain-ltr {{ background: #3498db; color: white; }}
.domain-transposase {{ background: #e67e22; color: white; }}
.domain-orf {{ background: #1abc9c; color: white; }}
.domain-unknown {{ background: #95a5a6; color: white; }}
table {{ border-collapse: collapse; margin: 15px 0; width: 100%; }}
th, td {{ padding: 8px 12px; border: 1px solid #ddd; text-align: left; }}
th {{ background: #f5f5f5; font-weight: 600; }}
tr:nth-child(even) {{ background: #fafafa; }}
.alignment-box {{
    background: #f8f9fa;
    padding: 12px;
    margin: 10px 0;
    border-radius: 4px;
    border-left: 4px solid #007bff;
}}
.alignment-box.strand-minus {{ border-left-color: #ffc107; }}
.align-header {{ font-weight: 600; margin-bottom: 5px; }}
.align-meta {{ color: #666; font-size: 12px; margin-bottom: 8px; }}
.align-seq {{
    white-space: pre;
    font-family: "SF Mono", monospace;
    font-size: 11px;
    line-height: 1.5;
}}
.match {{ color: #28a745; }}
.mismatch {{ color: #dc3545; }}
.gap {{ color: #6c757d; }}
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
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px;
    margin: 15px 0;
}}
.stat-card {{
    background: white;
    padding: 15px;
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}}
.stat-card h4 {{ margin: 0 0 10px 0; font-size: 14px; color: #666; }}
.nav-links {{
    margin: 20px 0;
    padding: 10px;
    background: white;
    border-radius: 4px;
}}
.nav-links a {{ margin-right: 15px; color: #0066cc; }}
</style>
</head>
<body>

<div class="nav-links">
<a href="index.html">← Back to Index</a>
</div>

<h1>{escape(symbol)} ({fbgn})</h1>
<p>Exon-level TE similarity analysis</p>

<div class="summary-box">
    <div class="summary-stat">
        <div class="label">Exons</div>
        <div class="value">{len(sorted_exons)}</div>
    </div>
    <div class="summary-stat">
        <div class="label">Total TE Hits</div>
        <div class="value">{total_hits:,}</div>
    </div>
    <div class="summary-stat">
        <div class="label">Sense Orientation</div>
        <div class="value">{100*sense_hits/max(total_hits,1):.1f}%</div>
    </div>
    <div class="summary-stat">
        <div class="label">TE Families</div>
        <div class="value">{len(te_family_counts)}</div>
    </div>
</div>

<h2>TE Domain Distribution</h2>
<div class="stats-grid">
'''

    for domain, count in top_domains:
        pct = 100 * count / total_hits if total_hits else 0
        desc = get_domain_description(domain)
        domain_class = 'domain-' + domain.split('_')[0] if domain != 'unknown' else 'domain-unknown'
        html += f'''<div class="stat-card">
    <h4><span class="domain-tag {domain_class}">{domain}</span></h4>
    <div style="font-size: 20px; font-weight: 600;">{count:,} hits ({pct:.1f}%)</div>
    <div style="color: #666; font-size: 12px;">{desc}</div>
</div>
'''

    html += '</div>\n'

    # Top TE families table
    html += '''<h2>Top TE Families</h2>
<table>
<tr><th>TE Family</th><th>Hits</th><th>Class</th><th>%</th></tr>
'''
    for te_id, count in top_families:
        te_class = infer_te_class(te_id)
        pct = 100 * count / total_hits if total_hits else 0
        html += f'<tr><td>{escape(te_id)}</td><td>{count:,}</td><td>{te_class}</td><td>{pct:.1f}%</td></tr>\n'

    html += '</table>\n'

    # ===========================================
    # FULL GENE VIEW WITH EXON BOUNDARIES
    # ===========================================
    html += '''
<h2>Full Gene View</h2>
<p>Concatenated exon sequences with boundaries marked and hits colored by TE domain.</p>
'''

    # Gene schematic
    total_len = sum(int(meta['length']) for _, meta in sorted_exons)
    html += '<div class="gene-schematic">\n'
    html += '<div style="font-weight: 600; margin-bottom: 10px;">Gene Structure</div>\n'
    html += '<div class="schematic-bar">\n'

    for i, (exon_id, meta) in enumerate(sorted_exons):
        length = int(meta['length'])
        width_pct = max(3, 100 * length / total_len)
        exon_num = meta['exon_num']
        position = meta['position']

        # Color by position
        colors = {
            'first_exon': '#28a745',
            'internal_exon': '#007bff',
            'last_exon': '#ffc107',
            'single_exon': '#dc3545'
        }
        color = colors.get(position, '#6c757d')

        html += f'<div class="schematic-exon" style="width:{width_pct}%; background:{color};" title="Exon {exon_num}: {length}bp">E{exon_num}</div>\n'

    html += '</div>\n'

    # Legend for exon colors
    html += '''<div style="display: flex; gap: 15px; margin-top: 10px; font-size: 12px;">
<span><span style="display:inline-block;width:12px;height:12px;background:#28a745;border-radius:2px;"></span> First exon</span>
<span><span style="display:inline-block;width:12px;height:12px;background:#007bff;border-radius:2px;"></span> Internal</span>
<span><span style="display:inline-block;width:12px;height:12px;background:#ffc107;border-radius:2px;"></span> Last exon</span>
<span><span style="display:inline-block;width:12px;height:12px;background:#dc3545;border-radius:2px;"></span> Single exon</span>
</div>
</div>
'''

    # Domain color legend
    html += '<div class="domain-legend">\n'
    html += '<strong style="margin-right: 10px;">TE Domain Colors:</strong>\n'
    domains_used = set(h['te_domain'] for h in hits)
    for domain in sorted(domains_used):
        color = DOMAIN_COLORS.get(domain, '#bdc3c7')
        html += f'<div class="domain-legend-item"><div class="domain-color-box" style="background:{color}"></div> {domain}</div>\n'
    html += '</div>\n'

    # Build full gene sequence with domain-colored hits
    html += '<div class="full-gene-view">\n'
    html += '<div style="font-weight: 600; margin-bottom: 10px;">Sequence with TE Domain Hits</div>\n'

    # Concatenate sequences and track exon boundaries
    full_seq = ''
    exon_boundaries = []  # List of (position, exon_num)
    exon_offsets = {}  # exon_id -> offset in concatenated seq

    for exon_id, meta in sorted_exons:
        exon_offsets[exon_id] = len(full_seq)
        seq = exon_seqs.get(exon_id, '')
        if exon_boundaries or full_seq:  # Not the first exon
            exon_boundaries.append((len(full_seq), meta['exon_num']))
        full_seq += seq

    # Build domain coverage map (position -> domain with highest priority)
    domain_priority = ['gag', 'pol', 'orf2_rt', 'orf2_en', 'orf1', 'transposase', '5_ltr', '3_ltr', '5_tir', '3_tir', '5_utr', '3_utr', 'unknown']
    position_domains = {}  # position -> domain

    for h in hits:
        exon_offset = exon_offsets.get(h['exon_id'], 0)
        domain = h['te_domain']
        for p in range(h['qstart'] - 1, h['qend']):
            global_pos = exon_offset + p
            if global_pos < len(full_seq):
                current = position_domains.get(global_pos)
                if current is None:
                    position_domains[global_pos] = domain
                else:
                    # Keep higher priority domain
                    try:
                        if domain_priority.index(domain) < domain_priority.index(current):
                            position_domains[global_pos] = domain
                    except ValueError:
                        pass

    # Render sequence with colors and boundaries (limit to 3000bp for display)
    html += '<div class="sequence-box">\n'
    display_len = min(len(full_seq), 3000)
    boundary_set = {b[0] for b in exon_boundaries}

    for line_start in range(0, display_len, 80):
        line_end = min(line_start + 80, display_len)
        html += f'<div class="seq-line"><span class="pos-num">{line_start+1:5d}</span>'

        for i in range(line_start, line_end):
            # Check for exon boundary
            is_boundary = i in boundary_set
            boundary_marker = ''
            if is_boundary:
                for bpos, exon_num in exon_boundaries:
                    if bpos == i:
                        boundary_marker = f'<span class="exon-boundary-marker">E{exon_num}</span>'
                        break

            domain = position_domains.get(i)
            if domain:
                color = DOMAIN_COLORS.get(domain, '#bdc3c7')
                if is_boundary:
                    html += f'{boundary_marker}<span class="exon-boundary" style="background-color:{color};">{full_seq[i]}</span>'
                else:
                    html += f'<span style="background-color:{color};">{full_seq[i]}</span>'
            else:
                if is_boundary:
                    html += f'{boundary_marker}<span class="exon-boundary">{full_seq[i]}</span>'
                else:
                    html += full_seq[i]

        html += '</div>\n'

    if len(full_seq) > 3000:
        html += f'<div class="seq-line"><em>... ({len(full_seq) - 3000} more bp not shown)</em></div>\n'

    html += '</div>\n'  # Close sequence-box
    html += '</div>\n'  # Close full-gene-view

    # ===========================================
    # PER-EXON DETAILS
    # ===========================================

    # Legend
    html += '''
<h2>Per-Exon Details</h2>
<div class="legend">
<div class="legend-item"><div class="legend-color" style="background:#cce5ff"></div> TE sense (+)</div>
<div class="legend-item"><div class="legend-color" style="background:#fff3cd"></div> TE antisense (-)</div>
<div class="legend-item"><div class="legend-color" style="background:#d4edda"></div> Both strands</div>
<div class="legend-item"><div class="legend-color" style="background:#dc3545"></div> Exon boundary</div>
</div>
'''

    # Each exon
    for exon_id, meta in sorted_exons:
        exon_hits = hits_by_exon.get(exon_id, [])
        seq = exon_seqs.get(exon_id, '')

        position = meta['position']
        pos_class = {
            'first_exon': 'pos-first',
            'internal_exon': 'pos-internal',
            'last_exon': 'pos-last',
            'single_exon': 'pos-single'
        }.get(position, 'pos-internal')

        utr_overlap = meta['utr_overlap']
        utr_class = {
            '5utr': 'utr-5',
            '3utr': 'utr-3',
            'both': 'utr-both',
            'none': 'utr-none'
        }.get(utr_overlap, 'utr-none')

        exon_num = meta['exon_num']
        length = int(meta['length'])
        exon_start = meta['start']
        exon_end = meta['end']

        html += f'''
<div class="exon-container">
<div class="exon-header">
    <div>
        <span class="exon-title">Exon {exon_num}</span>
        <span class="position-tag {pos_class}">{position.replace('_', ' ')}</span>
        <span class="utr-tag {utr_class}">UTR: {utr_overlap}</span>
    </div>
    <div class="exon-meta">
        {length:,} bp | {len(exon_hits):,} TE hits | {meta['chrom']}:{exon_start}-{exon_end}
    </div>
</div>
'''

        if seq and len(exon_hits) > 0:
            # Build coverage map
            sense_pos = set()
            anti_pos = set()
            for h in exon_hits:
                for p in range(h['qstart'] - 1, h['qend']):
                    if h['strand'] == '+':
                        sense_pos.add(p)
                    else:
                        anti_pos.add(p)

            # Render sequence (max 1500 bp shown)
            html += '<div class="sequence-box">\n'
            display_len = min(len(seq), 1500)

            for line_start in range(0, display_len, 80):
                line_end = min(line_start + 80, display_len)
                html += f'<div class="seq-line"><span class="pos-num">{line_start+1:5d}</span>'

                for i in range(line_start, line_end):
                    in_sense = i in sense_pos
                    in_anti = i in anti_pos

                    if in_sense and in_anti:
                        html += f'<span class="te-both">{seq[i]}</span>'
                    elif in_sense:
                        html += f'<span class="te-sense">{seq[i]}</span>'
                    elif in_anti:
                        html += f'<span class="te-anti">{seq[i]}</span>'
                    else:
                        html += seq[i]

                html += '</div>\n'

            if len(seq) > 1500:
                html += f'<div class="seq-line"><em>... ({len(seq) - 1500} more bp not shown)</em></div>\n'

            html += '</div>\n'

            # Top alignments for this exon
            if exon_hits:
                html += '<h4>Top Alignments</h4>\n'
                sorted_hits = sorted(exon_hits, key=lambda x: -x['bitscore'])[:10]

                for h in sorted_hits:
                    strand_class = '' if h['strand'] == '+' else 'strand-minus'
                    strand_sym = '+' if h['strand'] == '+' else '-'

                    domain = h['te_domain']
                    domain_class = 'domain-' + domain.split('_')[0] if domain != 'unknown' else 'domain-unknown'

                    # Color the alignment
                    qseq = h['qseq'][:60]
                    sseq = h['sseq'][:60]
                    colored_match = ''
                    for q, s in zip(qseq, sseq):
                        if q == s:
                            colored_match += f'<span class="match">|</span>'
                        elif q == '-' or s == '-':
                            colored_match += f'<span class="gap"> </span>'
                        else:
                            colored_match += f'<span class="mismatch">.</span>'

                    html += f'''<div class="alignment-box {strand_class}">
<div class="align-header">{escape(h['te_id'])} ({strand_sym}) <span class="domain-tag {domain_class}">{domain}</span></div>
<div class="align-meta">Position {h['qstart']}-{h['qend']} | {h['pident']:.1f}% identity | {h['length']}bp | E={h['evalue']:.2e}</div>
<div class="align-seq">Query: {escape(qseq)}{"..." if len(h['qseq']) > 60 else ""}
       {colored_match}
TE:    {escape(sseq)}{"..." if len(h['sseq']) > 60 else ""}</div>
</div>
'''

        html += '</div>\n'  # Close exon-container

    html += '</body>\n</html>'

    # Write file
    safe_symbol = symbol.replace('/', '_').replace(' ', '_')
    html_path = output_dir / f'{safe_symbol}_{fbgn}.html'
    with open(html_path, 'w') as f:
        f.write(html)

    return html_path.name


def generate_index_html(gene_info, output_dir):
    """Generate index page linking to all gene pages."""

    html = '''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>Top 50 Genes - Exon TE Analysis</title>
<style>
body { font-family: -apple-system, sans-serif; margin: 40px; max-width: 1200px; }
h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
table { border-collapse: collapse; margin: 20px 0; width: 100%; }
th, td { padding: 10px 16px; border: 1px solid #ddd; text-align: left; }
th { background: #f5f5f5; position: sticky; top: 0; }
tr:nth-child(even) { background: #fafafa; }
tr:hover { background: #e8f4e8; }
a { color: #0066cc; text-decoration: none; }
a:hover { text-decoration: underline; }
.high-density { background: #ffe6e6; }
.domain-tag { display: inline-block; padding: 2px 6px; border-radius: 3px; font-size: 10px; color: white; margin: 1px; }
.domain-gag { background: #e74c3c; }
.domain-pol { background: #9b59b6; }
.domain-ltr { background: #3498db; }
.domain-unknown { background: #95a5a6; }
</style>
</head>
<body>
<h1>Top 50 Genes by Exon TE Density</h1>
<p>Click on a gene to view detailed exon-level TE alignments and sequence visualizations.</p>
<table>
<tr>
<th>Rank</th>
<th>Gene</th>
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
        domain_class = 'domain-' + domain.split('_')[0] if domain != 'unknown' else 'domain-unknown'

        html += f'''<tr>
<td>{info['rank']}</td>
<td><strong>{escape(info['symbol'])}</strong><br><small>{info['fbgn']}</small></td>
<td>{info['exons']}</td>
<td>{info['hits']:,}</td>
<td>{info['density']:.1f}</td>
<td>{info['sense_pct']:.1f}%</td>
<td><span class="domain-tag {domain_class}">{domain}</span></td>
<td><a href="{info['html_file']}">View →</a></td>
</tr>
'''

    html += '</table>\n</body>\n</html>'

    with open(output_dir / 'index.html', 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--top', type=int, default=50, help='Number of top genes')
    parser.add_argument('--output-dir', type=Path,
                        default=Path('results/exon_analysis/html'),
                        help='Output directory')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    # Paths
    blast_path = Path('results/exon_analysis/genome_wide_exons.tsv')
    gene_stats_path = Path('results/exon_analysis/gene_exon_te_summary.tsv')
    exon_meta_path = Path('data/queries/genome_wide/exon_metadata.tsv')
    exon_fasta_path = Path('data/queries/genome_wide/exons_sense.fasta')

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load top genes
    print(f"Loading top {args.top} genes...")
    top_genes = {}
    with open(gene_stats_path) as f:
        f.readline()  # header
        for i, line in enumerate(f):
            if i >= args.top:
                break
            parts = line.strip().split('\t')
            fbgn = parts[1]
            symbol = parts[2]
            top_genes[fbgn] = {
                'rank': int(parts[0]),
                'symbol': symbol,
                'exons': int(parts[3]),
                'hits': int(parts[5]),
                'density': float(parts[8]),
                'sense_pct': float(parts[11]) if len(parts) > 11 else 50,
            }

    print(f"  Found {len(top_genes)} genes")

    # Load exon metadata
    print("Loading exon metadata...")
    exon_meta = load_exon_metadata(exon_meta_path)

    # Get exon IDs for target genes
    target_exon_ids = {eid for eid, meta in exon_meta.items()
                       if meta['fbgn'] in top_genes}
    print(f"  {len(target_exon_ids)} exons for target genes")

    # Load exon sequences
    print("Loading exon sequences...")
    exon_seqs = load_exon_sequences(exon_fasta_path, target_exon_ids)
    print(f"  Loaded {len(exon_seqs)} sequences")

    # Load BLAST hits
    print("Loading BLAST hits for target genes...")
    hits_by_gene = load_blast_hits_for_genes(blast_path, top_genes, exon_meta)
    print(f"  Loaded hits for {len(hits_by_gene)} genes")

    # Generate HTML for each gene
    print(f"\nGenerating HTML pages...")
    gene_info = []

    for fbgn, info in sorted(top_genes.items(), key=lambda x: x[1]['rank']):
        symbol = info['symbol']
        hits = hits_by_gene.get(fbgn, [])

        # Get top domain
        domain_counts = defaultdict(int)
        for h in hits:
            domain_counts[h['te_domain']] += 1
        top_domain = max(domain_counts, key=domain_counts.get) if domain_counts else 'unknown'

        html_file = generate_gene_html(fbgn, symbol, hits, exon_meta, exon_seqs, args.output_dir)

        gene_info.append({
            'rank': info['rank'],
            'fbgn': fbgn,
            'symbol': symbol,
            'exons': info['exons'],
            'hits': len(hits),
            'density': info['density'],
            'sense_pct': info['sense_pct'],
            'top_domain': top_domain,
            'html_file': html_file,
        })

        if args.verbose:
            print(f"  {info['rank']:3d}. {symbol}: {len(hits):,} hits -> {html_file}")

    # Generate index
    print("Generating index page...")
    generate_index_html(gene_info, args.output_dir)

    print(f"\nDone! Open {args.output_dir / 'index.html'} to view results.")


if __name__ == '__main__':
    main()
