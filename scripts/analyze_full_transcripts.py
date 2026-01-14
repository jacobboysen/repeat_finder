#!/usr/bin/env python3
"""
Analyze TE content across full transcripts (5'UTR + CDS + 3'UTR) for top genes.
"""

import subprocess
from pathlib import Path
from collections import defaultdict
import tempfile
import os

# Top 10 genes by TE density
TOP_GENES = [
    ('FBgn0040959', 'Prt-15a'),
    ('FBgn0034403', 'CG18190'),
    ('FBgn0067905', 'Dso2'),
    ('FBgn0033948', 'CG12863'),
    ('FBgn0053093', 'CG33093'),
    ('FBgn0037514', 'CG10919'),
    ('FBgn0260995', 'dpr21'),
    ('FBgn0040534', 'Sf3b5'),
    ('FBgn0053458', 'CG33458'),
    ('FBgn0054003', 'NimB3'),
]

def parse_fasta_by_parent(fasta_path, is_cds=False):
    """Parse FASTA file and index by parent gene (FBgn)."""
    sequences = defaultdict(dict)
    current_id = None
    current_seq = []
    current_parent = None
    current_fbtr = None
    current_len = 0

    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current_fbtr and current_parent:
                    sequences[current_parent][current_fbtr] = {
                        'seq': ''.join(current_seq),
                        'length': current_len
                    }

                parts = line.strip().split()
                current_id = parts[0][1:]  # May be FBtr or protein name
                current_seq = []
                current_parent = None
                current_fbtr = None
                current_len = 0

                # Parse header attributes
                header = line.strip()
                for part in parts:
                    if part.startswith('parent='):
                        parent_val = part.split('=')[1].rstrip(';')
                        # Parent may have multiple values like "FBgn0040959,FBtr0079713"
                        for p in parent_val.split(','):
                            if p.startswith('FBgn'):
                                current_parent = p
                            elif p.startswith('FBtr'):
                                current_fbtr = p
                    if part.startswith('length='):
                        current_len = int(part.split('=')[1].rstrip(';'))

                # If ID is FBtr, use it
                if current_id.startswith('FBtr'):
                    current_fbtr = current_id
            else:
                current_seq.append(line.strip())

        if current_fbtr and current_parent:
            sequences[current_parent][current_fbtr] = {
                'seq': ''.join(current_seq),
                'length': current_len
            }

    return sequences


def build_full_transcripts(utr5_seqs, cds_seqs, utr3_seqs, genes):
    """Build full transcripts by concatenating 5'UTR + CDS + 3'UTR."""
    transcripts = {}

    for fbgn, symbol in genes:
        # Find transcripts that have all three regions
        # Or at least CDS + 3'UTR

        cds_trs = set(cds_seqs.get(fbgn, {}).keys())
        utr3_trs = set(utr3_seqs.get(fbgn, {}).keys())
        utr5_trs = set(utr5_seqs.get(fbgn, {}).keys())

        # Find transcripts with at least CDS and 3'UTR
        common_trs = cds_trs & utr3_trs

        if not common_trs:
            print(f"  Warning: No complete transcripts for {symbol} ({fbgn})")
            # Try to use just 3'UTR if available
            if utr3_trs:
                tr = list(utr3_trs)[0]
                transcripts[fbgn] = {
                    'symbol': symbol,
                    'transcript': tr,
                    'utr5_len': 0,
                    'cds_len': 0,
                    'utr3_len': utr3_seqs[fbgn][tr]['length'],
                    'utr5_seq': '',
                    'cds_seq': '',
                    'utr3_seq': utr3_seqs[fbgn][tr]['seq'],
                    'full_seq': utr3_seqs[fbgn][tr]['seq'],
                }
            continue

        # Pick the transcript with longest combined length
        best_tr = None
        best_len = 0

        for tr in common_trs:
            utr5_len = utr5_seqs.get(fbgn, {}).get(tr, {}).get('length', 0)
            cds_len = cds_seqs[fbgn][tr]['length']
            utr3_len = utr3_seqs[fbgn][tr]['length']
            total = utr5_len + cds_len + utr3_len

            if total > best_len:
                best_len = total
                best_tr = tr

        tr = best_tr
        utr5_seq = utr5_seqs.get(fbgn, {}).get(tr, {}).get('seq', '')
        cds_seq = cds_seqs[fbgn][tr]['seq']
        utr3_seq = utr3_seqs[fbgn][tr]['seq']

        transcripts[fbgn] = {
            'symbol': symbol,
            'transcript': tr,
            'utr5_len': len(utr5_seq),
            'cds_len': len(cds_seq),
            'utr3_len': len(utr3_seq),
            'utr5_seq': utr5_seq,
            'cds_seq': cds_seq,
            'utr3_seq': utr3_seq,
            'full_seq': utr5_seq + cds_seq + utr3_seq,
        }

        print(f"  {symbol}: 5'UTR={len(utr5_seq)}bp, CDS={len(cds_seq)}bp, 3'UTR={len(utr3_seq)}bp, Total={len(utr5_seq)+len(cds_seq)+len(utr3_seq)}bp")

    return transcripts


def run_blast(transcripts, te_db, blast_path):
    """Run BLAST for full transcripts against TE database."""
    results = {}

    for fbgn, data in transcripts.items():
        symbol = data['symbol']
        seq = data['full_seq']

        if not seq:
            continue

        # Write query to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">{fbgn}\n{seq}\n")
            query_file = f.name

        # Run BLAST
        cmd = [
            blast_path,
            '-query', query_file,
            '-db', te_db,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq',
            '-evalue', '10',
            '-word_size', '7',
            '-gapopen', '2',
            '-gapextend', '1',
            '-penalty', '-1',
            '-reward', '1',
            '-dust', 'yes',
            '-max_target_seqs', '1000',
            '-max_hsps', '100',
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            hits = []

            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 16:
                    hits.append({
                        'te_id': parts[1],
                        'pident': float(parts[2]),
                        'length': int(parts[3]),
                        'qstart': int(parts[6]),
                        'qend': int(parts[7]),
                        'sstart': int(parts[8]),
                        'send': int(parts[9]),
                        'evalue': float(parts[10]),
                        'bitscore': float(parts[11]),
                        'qseq': parts[14],
                        'sseq': parts[15],
                        'strand': '+' if int(parts[8]) < int(parts[9]) else '-'
                    })

            results[fbgn] = hits
            print(f"  {symbol}: {len(hits)} hits")

        except Exception as e:
            print(f"  {symbol}: BLAST error - {e}")
            results[fbgn] = []

        os.unlink(query_file)

    return results


def annotate_hits(transcripts, blast_results):
    """Annotate each hit with which transcript region it falls in."""
    annotated = {}

    for fbgn, hits in blast_results.items():
        data = transcripts[fbgn]
        utr5_end = data['utr5_len']
        cds_end = utr5_end + data['cds_len']
        utr3_end = cds_end + data['utr3_len']

        annotated_hits = []
        for hit in hits:
            qstart = hit['qstart']
            qend = hit['qend']

            # Determine which region(s) the hit spans
            regions = []
            region_bp = {'5UTR': 0, 'CDS': 0, '3UTR': 0}

            # 5'UTR portion
            if qstart <= utr5_end:
                utr5_start = qstart
                utr5_stop = min(qend, utr5_end)
                bp = utr5_stop - utr5_start + 1
                if bp > 0:
                    regions.append('5UTR')
                    region_bp['5UTR'] = bp

            # CDS portion
            if qend > utr5_end and qstart <= cds_end:
                cds_start = max(qstart, utr5_end + 1)
                cds_stop = min(qend, cds_end)
                bp = cds_stop - cds_start + 1
                if bp > 0:
                    regions.append('CDS')
                    region_bp['CDS'] = bp

            # 3'UTR portion
            if qend > cds_end:
                utr3_start = max(qstart, cds_end + 1)
                utr3_stop = min(qend, utr3_end)
                bp = utr3_stop - utr3_start + 1
                if bp > 0:
                    regions.append('3UTR')
                    region_bp['3UTR'] = bp

            hit['regions'] = regions
            hit['region_bp'] = region_bp
            hit['primary_region'] = max(region_bp, key=region_bp.get) if region_bp else 'Unknown'
            annotated_hits.append(hit)

        annotated[fbgn] = annotated_hits

    return annotated


def generate_html(transcripts, annotated_hits, output_dir):
    """Generate HTML visualizations."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for fbgn, data in transcripts.items():
        symbol = data['symbol']
        hits = annotated_hits.get(fbgn, [])

        utr5_len = data['utr5_len']
        cds_len = data['cds_len']
        utr3_len = data['utr3_len']
        total_len = utr5_len + cds_len + utr3_len

        if total_len == 0:
            continue

        # Build position coverage maps
        utr5_cov = [0] * utr5_len
        cds_cov = [0] * cds_len
        utr3_cov = [0] * utr3_len

        sense_pos = set()
        anti_pos = set()

        for hit in hits:
            for p in range(hit['qstart'] - 1, hit['qend']):
                if hit['strand'] == '+':
                    sense_pos.add(p)
                else:
                    anti_pos.add(p)

                # Count coverage
                if p < utr5_len:
                    utr5_cov[p] += 1
                elif p < utr5_len + cds_len:
                    cds_cov[p - utr5_len] += 1
                else:
                    idx = p - utr5_len - cds_len
                    if idx < utr3_len:
                        utr3_cov[idx] += 1

        # Calculate coverage stats
        utr5_covered = sum(1 for x in utr5_cov if x > 0) if utr5_cov else 0
        cds_covered = sum(1 for x in cds_cov if x > 0) if cds_cov else 0
        utr3_covered = sum(1 for x in utr3_cov if x > 0) if utr3_cov else 0

        utr5_pct = 100 * utr5_covered / utr5_len if utr5_len else 0
        cds_pct = 100 * cds_covered / cds_len if cds_len else 0
        utr3_pct = 100 * utr3_covered / utr3_len if utr3_len else 0

        # Count hits by region
        hits_5utr = sum(1 for h in hits if '5UTR' in h['regions'])
        hits_cds = sum(1 for h in hits if 'CDS' in h['regions'])
        hits_3utr = sum(1 for h in hits if '3UTR' in h['regions'])

        # Strand stats
        sense_count = sum(1 for h in hits if h['strand'] == '+')
        anti_count = len(hits) - sense_count

        # TE family stats
        te_counts = defaultdict(int)
        for h in hits:
            te_counts[h['te_id']] += 1
        top_tes = sorted(te_counts.items(), key=lambda x: x[1], reverse=True)[:10]

        # Generate HTML
        html = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>{symbol} ({fbgn}) Full Transcript TE Analysis</title>
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: 14px;
    line-height: 1.6;
    margin: 40px;
    background: #fafafa;
    color: #333;
}}
h1, h2, h3 {{ font-weight: 600; }}
h1 {{ border-bottom: 2px solid #333; padding-bottom: 10px; }}
.summary-box {{
    background: #e8f4e8;
    padding: 15px 20px;
    border-radius: 8px;
    margin: 20px 0;
}}
.region-box {{
    display: inline-block;
    padding: 15px 20px;
    border-radius: 8px;
    margin: 5px;
    text-align: center;
    min-width: 120px;
}}
.region-5utr {{ background: #d4edda; border: 2px solid #28a745; }}
.region-cds {{ background: #cce5ff; border: 2px solid #007bff; }}
.region-3utr {{ background: #fff3cd; border: 2px solid #ffc107; }}
.transcript-map {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    margin: 20px 0;
    overflow-x: auto;
}}
.region-bar {{
    height: 40px;
    display: inline-block;
    text-align: center;
    line-height: 40px;
    font-weight: 600;
    font-size: 12px;
    color: white;
}}
.bar-5utr {{ background: #28a745; }}
.bar-cds {{ background: #007bff; }}
.bar-3utr {{ background: #ffc107; color: #333; }}
.coverage-track {{
    height: 20px;
    background: #f0f0f0;
    margin-top: 5px;
    position: relative;
}}
.coverage-fill {{
    height: 100%;
    position: absolute;
    top: 0;
}}
.fill-5utr {{ background: rgba(40, 167, 69, 0.7); }}
.fill-cds {{ background: rgba(0, 123, 255, 0.7); }}
.fill-3utr {{ background: rgba(255, 193, 7, 0.7); }}
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
    font-size: 12px;
    line-height: 1.3;
}}
.pos-num {{ color: #888; margin-right: 10px; }}
.te-sense {{ background-color: #cce5ff; }}
.te-anti {{ background-color: #fff3cd; }}
.te-both {{ background-color: #d4edda; }}
.utr5-char {{ border-bottom: 3px solid #28a745; }}
.cds-char {{ border-bottom: 3px solid #007bff; }}
.utr3-char {{ border-bottom: 3px solid #ffc107; }}
table {{ border-collapse: collapse; margin: 20px 0; width: 100%; max-width: 900px; }}
th, td {{ padding: 8px 12px; border: 1px solid #ddd; text-align: left; }}
th {{ background: #f5f5f5; font-weight: 600; }}
tr:nth-child(even) {{ background: #fafafa; }}
.alignment-box {{
    background: white;
    padding: 15px;
    margin: 15px 0;
    border-radius: 4px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
}}
.alignment-box.in-5utr {{ border-left: 4px solid #28a745; }}
.alignment-box.in-cds {{ border-left: 4px solid #007bff; }}
.alignment-box.in-3utr {{ border-left: 4px solid #ffc107; }}
.alignment-box.multi-region {{ border-left: 4px solid #dc3545; }}
.align-header {{ font-weight: 600; margin-bottom: 8px; }}
.align-meta {{ color: #666; font-size: 12px; margin-bottom: 10px; }}
.align-seq {{
    white-space: pre;
    font-family: "SF Mono", monospace;
    font-size: 11px;
    line-height: 1.4;
}}
.region-tag {{
    display: inline-block;
    padding: 2px 8px;
    border-radius: 3px;
    font-size: 11px;
    font-weight: 600;
    color: white;
    margin-right: 4px;
}}
.tag-5utr {{ background: #28a745; }}
.tag-cds {{ background: #007bff; }}
.tag-3utr {{ background: #ffc107; color: #333; }}
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
</style>
</head>
<body>

<h1>{symbol} ({fbgn}) - Full Transcript TE Analysis</h1>

<div class="summary-box">
<strong>Transcript:</strong> {data['transcript']}<br>
<strong>Total length:</strong> {total_len:,} bp<br>
<strong>Total TE hits:</strong> {len(hits):,}<br>
<strong>Strand bias:</strong> {sense_count} sense ({100*sense_count/max(len(hits),1):.1f}%) / {anti_count} antisense ({100*anti_count/max(len(hits),1):.1f}%)
</div>

<h2>Transcript Structure & TE Coverage</h2>

<div class="region-box region-5utr">
    <strong>5' UTR</strong><br>
    {utr5_len:,} bp<br>
    {utr5_pct:.1f}% covered<br>
    {hits_5utr} hits
</div>
<div class="region-box region-cds">
    <strong>CDS</strong><br>
    {cds_len:,} bp<br>
    {cds_pct:.1f}% covered<br>
    {hits_cds} hits
</div>
<div class="region-box region-3utr">
    <strong>3' UTR</strong><br>
    {utr3_len:,} bp<br>
    {utr3_pct:.1f}% covered<br>
    {hits_3utr} hits
</div>

<div class="transcript-map">
<h3>Transcript Map (to scale)</h3>
'''

        # Draw transcript bar
        if total_len > 0:
            utr5_width = max(0.5, 100 * utr5_len / total_len) if utr5_len else 0
            cds_width = max(0.5, 100 * cds_len / total_len) if cds_len else 0
            utr3_width = max(0.5, 100 * utr3_len / total_len) if utr3_len else 0

            html += f'''<div style="width: 100%; max-width: 800px;">
<div class="region-bar bar-5utr" style="width: {utr5_width}%;">{utr5_len}bp</div>'''
            html += f'<div class="region-bar bar-cds" style="width: {cds_width}%;">{cds_len}bp</div>'
            html += f'<div class="region-bar bar-3utr" style="width: {utr3_width}%;">{utr3_len}bp</div>'

            # Coverage track
            html += '<div class="coverage-track">'
            if utr5_len:
                html += f'<div class="coverage-fill fill-5utr" style="left: 0; width: {utr5_width * utr5_pct / 100}%;"></div>'
            if cds_len:
                html += f'<div class="coverage-fill fill-cds" style="left: {utr5_width}%; width: {cds_width * cds_pct / 100}%;"></div>'
            if utr3_len:
                html += f'<div class="coverage-fill fill-3utr" style="left: {utr5_width + cds_width}%; width: {utr3_width * utr3_pct / 100}%;"></div>'
            html += '</div>'
            html += '<div style="font-size: 11px; color: #666; margin-top: 5px;">Coverage track shows TE-covered positions</div>'
            html += '</div>'

        html += '</div>'

        # TE family breakdown
        html += '''
<h2>Top TE Families</h2>
<table>
<tr><th>TE</th><th>Hits</th><th>Primary Region</th></tr>
'''
        for te_id, count in top_tes:
            # Find primary region for this TE
            te_regions = defaultdict(int)
            for h in hits:
                if h['te_id'] == te_id:
                    for r in h['regions']:
                        te_regions[r] += 1
            primary = max(te_regions, key=te_regions.get) if te_regions else 'Unknown'
            tag_class = {'5UTR': 'tag-5utr', 'CDS': 'tag-cds', '3UTR': 'tag-3utr'}.get(primary, '')
            html += f'<tr><td>{te_id}</td><td>{count}</td><td><span class="region-tag {tag_class}">{primary}</span></td></tr>\n'

        html += '</table>'

        # Sequence visualization
        html += '''
<h2>Sequence with TE Annotations</h2>
<div class="legend">
<div class="legend-item"><div class="legend-color" style="background:#cce5ff"></div> TE sense (+)</div>
<div class="legend-item"><div class="legend-color" style="background:#fff3cd"></div> TE antisense (-)</div>
<div class="legend-item"><div class="legend-color" style="background:#d4edda"></div> Both strands</div>
<div class="legend-item"><div class="legend-color" style="background:#28a745"></div> 5' UTR</div>
<div class="legend-item"><div class="legend-color" style="background:#007bff"></div> CDS</div>
<div class="legend-item"><div class="legend-color" style="background:#ffc107"></div> 3' UTR</div>
</div>
'''

        html += '<div class="sequence-container">\n'

        full_seq = data['full_seq']
        display_len = min(len(full_seq), 3000)

        for line_start in range(0, display_len, 80):
            line_end = min(line_start + 80, display_len)
            html += f'<div class="seq-line"><span class="pos-num">{line_start+1:5d}</span>'

            for i in range(line_start, line_end):
                in_sense = i in sense_pos
                in_anti = i in anti_pos

                # Determine region
                if i < utr5_len:
                    region_class = 'utr5-char'
                elif i < utr5_len + cds_len:
                    region_class = 'cds-char'
                else:
                    region_class = 'utr3-char'

                if in_sense and in_anti:
                    html += f'<span class="te-both {region_class}">{full_seq[i]}</span>'
                elif in_sense:
                    html += f'<span class="te-sense {region_class}">{full_seq[i]}</span>'
                elif in_anti:
                    html += f'<span class="te-anti {region_class}">{full_seq[i]}</span>'
                else:
                    html += f'<span class="{region_class}">{full_seq[i]}</span>'

            html += '</div>\n'

        if len(full_seq) > 3000:
            html += f'<div class="seq-line"><em>... ({len(full_seq) - 3000} more bp not shown)</em></div>\n'

        html += '</div>\n'

        # Top alignments
        html += '<h2>Top Alignments</h2>\n'

        # Sort by bitscore
        sorted_hits = sorted(hits, key=lambda x: x['bitscore'], reverse=True)[:20]

        for h in sorted_hits:
            regions = h['regions']
            if len(regions) > 1:
                box_class = 'multi-region'
            elif '5UTR' in regions:
                box_class = 'in-5utr'
            elif 'CDS' in regions:
                box_class = 'in-cds'
            else:
                box_class = 'in-3utr'

            region_tags = ''
            for r in regions:
                tag_class = {'5UTR': 'tag-5utr', 'CDS': 'tag-cds', '3UTR': 'tag-3utr'}.get(r, '')
                region_tags += f'<span class="region-tag {tag_class}">{r}</span>'

            strand_sym = '+' if h['strand'] == '+' else '-'

            html += f'''<div class="alignment-box {box_class}">
<div class="align-header">{h['te_id']} ({strand_sym}) {region_tags}</div>
<div class="align-meta">Position {h['qstart']}-{h['qend']} | {h['pident']:.1f}% identity | {h['length']}bp | E={h['evalue']:.2e}</div>
<div class="align-seq">Query: {h['qseq'][:80]}{"..." if len(h['qseq']) > 80 else ""}
TE:    {h['sseq'][:80]}{"..." if len(h['sseq']) > 80 else ""}</div>
</div>
'''

        html += '</body>\n</html>'

        # Write HTML
        html_path = output_dir / f'{symbol}_{fbgn}_full_transcript.html'
        with open(html_path, 'w') as f:
            f.write(html)

        print(f"  Created {html_path.name}")

    # Create index page
    with open(output_dir / 'index.html', 'w') as f:
        f.write('''<!DOCTYPE html>
<html>
<head>
<title>Full Transcript TE Analysis - Top 10 Genes</title>
<style>
body { font-family: -apple-system, sans-serif; margin: 40px; }
h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
table { border-collapse: collapse; margin: 20px 0; }
th, td { padding: 10px 16px; border: 1px solid #ddd; text-align: left; }
th { background: #f5f5f5; }
a { color: #0066cc; }
.high { background: #ffe6e6; }
.region-tag { display: inline-block; padding: 2px 6px; border-radius: 3px; font-size: 11px; color: white; margin: 1px; }
.tag-5utr { background: #28a745; }
.tag-cds { background: #007bff; }
.tag-3utr { background: #ffc107; color: #333; }
</style>
</head>
<body>
<h1>Full Transcript TE Analysis - Top 10 Genes by 3'UTR TE Density</h1>
<p>These genes have the highest TE content in their 3'UTRs. This analysis shows TE coverage across the entire transcript (5'UTR + CDS + 3'UTR) to determine if TE similarity extends beyond the UTR.</p>
<table>
<tr><th>Rank</th><th>Gene</th><th>5'UTR</th><th>CDS</th><th>3'UTR</th><th>Total Hits</th><th>View</th></tr>
''')
        for i, (fbgn, symbol) in enumerate(TOP_GENES, 1):
            if fbgn in transcripts:
                data = transcripts[fbgn]
                hits = annotated_hits.get(fbgn, [])

                utr5_pct = 100 * sum(1 for h in hits if '5UTR' in h['regions']) / max(len(hits), 1)
                cds_pct = 100 * sum(1 for h in hits if 'CDS' in h['regions']) / max(len(hits), 1)
                utr3_pct = 100 * sum(1 for h in hits if '3UTR' in h['regions']) / max(len(hits), 1)

                f.write(f'<tr class="high"><td>{i}</td><td><strong>{symbol}</strong><br><small>{fbgn}</small></td>')
                f.write(f'<td>{data["utr5_len"]}bp<br><span class="region-tag tag-5utr">{utr5_pct:.0f}%</span></td>')
                f.write(f'<td>{data["cds_len"]}bp<br><span class="region-tag tag-cds">{cds_pct:.0f}%</span></td>')
                f.write(f'<td>{data["utr3_len"]}bp<br><span class="region-tag tag-3utr">{utr3_pct:.0f}%</span></td>')
                f.write(f'<td>{len(hits)}</td>')
                f.write(f'<td><a href="{symbol}_{fbgn}_full_transcript.html">View</a></td></tr>\n')

        f.write('</table>\n</body>\n</html>')

    print(f"\nCreated index at {output_dir / 'index.html'}")


def main():
    base_dir = Path('/Users/jacobboysen/git_repos/repeat_finder')
    blast_path = '/Users/jacobboysen/miniconda3/envs/bioinformatics-program/bin/blastn'
    te_db = str(base_dir / 'data/blastdb/dmel_te_flybase')
    output_dir = base_dir / 'results/full_transcript_te'

    print("Loading sequence data...")
    utr5_seqs = parse_fasta_by_parent(base_dir / 'data/references/dmel_5utr.fasta')
    cds_seqs = parse_fasta_by_parent(base_dir / 'data/references/dmel_cds.fasta')
    utr3_seqs = parse_fasta_by_parent(base_dir / 'data/references/dmel_3utr.fasta')

    print(f"\nBuilding full transcripts for top 10 genes...")
    transcripts = build_full_transcripts(utr5_seqs, cds_seqs, utr3_seqs, TOP_GENES)

    print(f"\nRunning BLAST against TE database...")
    blast_results = run_blast(transcripts, te_db, blast_path)

    print(f"\nAnnotating hits by transcript region...")
    annotated = annotate_hits(transcripts, blast_results)

    print(f"\nGenerating HTML visualizations...")
    generate_html(transcripts, annotated, output_dir)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY: TE Coverage by Transcript Region")
    print("="*70)
    print(f"{'Gene':<12} {'5UTR':>10} {'CDS':>10} {'3UTR':>10} {'Spanning':>10}")
    print("-"*55)

    for fbgn, symbol in TOP_GENES:
        if fbgn not in annotated:
            continue
        hits = annotated[fbgn]

        only_5utr = sum(1 for h in hits if h['regions'] == ['5UTR'])
        only_cds = sum(1 for h in hits if h['regions'] == ['CDS'])
        only_3utr = sum(1 for h in hits if h['regions'] == ['3UTR'])
        spanning = sum(1 for h in hits if len(h['regions']) > 1)

        print(f"{symbol:<12} {only_5utr:>10} {only_cds:>10} {only_3utr:>10} {spanning:>10}")

    print("="*70)


if __name__ == '__main__':
    main()
