#!/usr/bin/env python3
"""
Generate comprehensive TE fossil analysis report.

Creates HTML and Markdown reports summarizing:
- Gene list and methodology
- Parameter optimization results
- TE signal analysis by gene
- Candidate TE fossils
- Follow-up recommendations
"""

import argparse
import base64
import json
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def load_gene_list(gene_list_file):
    """Load consolidated gene list."""
    if not gene_list_file.exists():
        return None
    return pd.read_csv(gene_list_file, sep='\t')


def load_density_data(density_file):
    """Load density data from pickle."""
    if not density_file.exists():
        return None
    with open(density_file, 'rb') as f:
        return pickle.load(f)


def load_cluster_data(cluster_file):
    """Load cluster summary."""
    if not cluster_file.exists():
        return None
    return pd.read_csv(cluster_file, sep='\t')


def load_sweep_summary(sweep_dir):
    """Load parameter sweep summary."""
    summary_file = sweep_dir / 'sweep_summary.tsv'
    if not summary_file.exists():
        return None
    return pd.read_csv(summary_file, sep='\t')


def image_to_base64(image_path):
    """Convert image to base64 for embedding in HTML."""
    if not image_path.exists():
        return None
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def generate_markdown_report(output_file, data):
    """Generate Markdown report."""
    lines = []

    # Header
    lines.append("# TE Fossil Analysis Report")
    lines.append("")
    lines.append(f"**Generated:** {data['timestamp']}")
    lines.append("")

    # Executive Summary
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("This report summarizes the search for transposable element (TE) \"fossils\" ")
    lines.append("in the 3' untranslated regions (3'UTRs) of Drosophila germ plasm genes.")
    lines.append("")

    if data.get('top_genes'):
        lines.append("### Key Findings")
        lines.append("")
        lines.append("| Gene | Max Signal | Total Hits | (+) Strand | (-) Strand | Top TE Family |")
        lines.append("|------|------------|------------|------------|------------|---------------|")
        for gene in data['top_genes'][:5]:
            plus_hits = gene.get('plus_strand_hits', '-')
            minus_hits = gene.get('minus_strand_hits', '-')
            lines.append(f"| {gene['gene']} | {gene['max_signal']:.2e} | {gene['total_hits']} | {plus_hits} | {minus_hits} | {gene.get('top_te', 'N/A')} |")
        lines.append("")

    if data.get('clusters_summary'):
        lines.append(f"**Total Clusters Detected:** {data['clusters_summary']['total_clusters']}")
        lines.append(f"**Genes with Clusters:** {data['clusters_summary']['genes_with_clusters']}")
        lines.append("")

    # Methods
    lines.append("## Methods")
    lines.append("")
    lines.append("### Gene Selection")
    lines.append("")
    if data.get('gene_list') is not None and len(data['gene_list']) > 0:
        lines.append(f"- **Total genes analyzed:** {len(data['gene_list'])}")
        lines.append("- **Source:** Canonical germ plasm genes from literature")
        lines.append("- **Tier 1 genes (primary analysis):** nos, osk, pgc, gcl, vas, aub, CycB")
        lines.append("")

    lines.append("### BLAST Parameters")
    lines.append("")
    lines.append("Sensitive parameters optimized for detecting ancient/degraded TE sequences:")
    lines.append("")
    lines.append("- **E-value threshold:** 10 (permissive)")
    lines.append("- **Word size:** 7 (short seeds for sensitive detection)")
    lines.append("- **Dust filter:** Disabled")
    lines.append("- **Soft masking:** Disabled")
    lines.append("")

    if data.get('best_params'):
        lines.append("### Optimal Parameter Combination")
        lines.append("")
        lines.append("```")
        for key, value in data['best_params'].items():
            lines.append(f"{key}: {value}")
        lines.append("```")
        lines.append("")

    # Results by Gene
    lines.append("## Results by Gene")
    lines.append("")

    if data.get('gene_results'):
        for gene_data in data['gene_results']:
            lines.append(f"### {gene_data['gene']}")
            lines.append("")
            lines.append(f"- **3'UTR Length:** {gene_data['length']} bp")
            lines.append(f"- **Total BLAST Hits:** {gene_data['total_hits']}")
            plus_hits = gene_data.get('plus_strand_hits', 0)
            minus_hits = gene_data.get('minus_strand_hits', 0)
            lines.append(f"  - Plus strand (TE sense): {plus_hits}")
            lines.append(f"  - Minus strand (TE antisense): {minus_hits}")
            lines.append(f"- **Maximum Signal:** {gene_data['max_signal']:.2e}")
            lines.append(f"- **Clusters Detected:** {gene_data.get('num_clusters', 0)}")
            lines.append("")

            if gene_data.get('clusters'):
                lines.append("**Cluster Details:**")
                lines.append("")
                lines.append("| Position | Length | Signal | Dominant TE |")
                lines.append("|----------|--------|--------|-------------|")
                for cluster in gene_data['clusters'][:5]:
                    lines.append(f"| {cluster['start']}-{cluster['end']} | {cluster['length']} bp | {cluster['max_signal']:.2e} | {cluster.get('dominant_te', 'N/A')} |")
                lines.append("")

    # TE Family Analysis
    lines.append("## TE Family Analysis")
    lines.append("")
    if data.get('te_families'):
        lines.append("### Top TE Families (by hit count)")
        lines.append("")
        lines.append("| TE Family | Total Hits | Genes with Hits |")
        lines.append("|-----------|------------|-----------------|")
        for te_data in data['te_families'][:10]:
            lines.append(f"| {te_data['family']} | {te_data['count']} | {te_data.get('genes', 'N/A')} |")
        lines.append("")

    # Candidate TE Fossils
    lines.append("## Candidate TE Fossils")
    lines.append("")
    lines.append("The following regions show concentrated TE signal and are candidates for ")
    lines.append("ancient TE insertions that may have been co-opted for regulatory function:")
    lines.append("")

    if data.get('top_clusters'):
        lines.append("| Gene | Position | Length | Signal Score | Dominant TE |")
        lines.append("|------|----------|--------|--------------|-------------|")
        for cluster in data['top_clusters'][:10]:
            lines.append(f"| {cluster['gene']} | {cluster['start']}-{cluster['end']} | {cluster['length']} bp | {cluster['max_signal']:.2e} | {cluster.get('dominant_te', 'N/A')} |")
        lines.append("")

    # Follow-up Recommendations
    lines.append("## Follow-up Analysis Recommendations")
    lines.append("")
    lines.append("1. **Secondary Structure Analysis:** Run RNAfold on candidate TE fossil regions")
    lines.append("2. **Conservation Analysis:** Check if TE regions are conserved across Drosophilids")
    lines.append("3. **Motif Overlap:** Cross-reference with known localization elements")
    lines.append("4. **Experimental Validation:** Reporter assays for candidate regions")
    lines.append("")

    # Technical Notes
    lines.append("## Technical Notes")
    lines.append("")
    lines.append("### Signal Density Calculation")
    lines.append("")
    lines.append("The combined signal score at each position is calculated as:")
    lines.append("")
    lines.append("```")
    lines.append("signal[i] = Σ (pident/100 × 1/e-value × 1/length)")
    lines.append("```")
    lines.append("")
    lines.append("for all BLAST HSPs overlapping position i, with Gaussian smoothing (σ=25bp).")
    lines.append("")
    lines.append("### Strand Terminology")
    lines.append("")
    lines.append("BLAST hits are classified by their orientation on the subject (TE) sequence:")
    lines.append("")
    lines.append("- **Plus strand (+)**: Query matches the sense strand of the TE (sstart < send)")
    lines.append("- **Minus strand (-)**: Query matches the antisense strand of the TE (sstart > send)")
    lines.append("")
    lines.append("For true TE fossil detection, plus strand hits from sense 3'UTR queries are most informative,")
    lines.append("as they indicate the 3'UTR contains sequence matching the functional orientation of the TE.")
    lines.append("")

    # Footer
    lines.append("---")
    lines.append("")
    lines.append("*Report generated by TE Fossil Mining Pipeline*")
    lines.append("")

    # Write file
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))


def generate_html_report(output_file, data, figures_dir=None):
    """Generate HTML report with embedded figures."""
    html = []

    # HTML header
    html.append("""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>TE Fossil Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
        h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        h3 { color: #7f8c8d; }
        table { border-collapse: collapse; width: 100%; margin: 15px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .summary-box { background-color: #ecf0f1; padding: 15px; border-radius: 5px; margin: 15px 0; }
        .highlight { background-color: #f1c40f; padding: 2px 5px; }
        img { max-width: 100%; height: auto; margin: 15px 0; }
        code { background-color: #f4f4f4; padding: 2px 5px; border-radius: 3px; }
        pre { background-color: #f4f4f4; padding: 15px; border-radius: 5px; overflow-x: auto; }
        .figure-container { text-align: center; margin: 20px 0; }
        .figure-caption { font-style: italic; color: #666; }
    </style>
</head>
<body>
""")

    # Title
    html.append(f"<h1>TE Fossil Analysis Report</h1>")
    html.append(f"<p><strong>Generated:</strong> {data['timestamp']}</p>")

    # Executive Summary
    html.append("<h2>Executive Summary</h2>")
    html.append("<div class='summary-box'>")
    html.append("<p>This report summarizes the search for transposable element (TE) \"fossils\" ")
    html.append("in the 3' untranslated regions (3'UTRs) of Drosophila germ plasm genes.</p>")

    if data.get('top_genes'):
        html.append("<h3>Key Findings</h3>")
        html.append("<table>")
        html.append("<tr><th>Gene</th><th>Max Signal</th><th>Total Hits</th><th>(+) Strand</th><th>(-) Strand</th><th>Top TE Family</th></tr>")
        for gene in data['top_genes'][:5]:
            plus_hits = gene.get('plus_strand_hits', '-')
            minus_hits = gene.get('minus_strand_hits', '-')
            html.append(f"<tr><td>{gene['gene']}</td><td>{gene['max_signal']:.2e}</td>")
            html.append(f"<td>{gene['total_hits']}</td><td>{plus_hits}</td><td>{minus_hits}</td><td>{gene.get('top_te', 'N/A')}</td></tr>")
        html.append("</table>")

    if data.get('clusters_summary'):
        html.append(f"<p><strong>Total Clusters Detected:</strong> {data['clusters_summary']['total_clusters']}</p>")
        html.append(f"<p><strong>Genes with Clusters:</strong> {data['clusters_summary']['genes_with_clusters']}</p>")

    html.append("</div>")

    # Embed figures if available
    if figures_dir:
        gene_comparison = figures_dir / 'gene_comparison' / 'gene_ranking.png'
        if gene_comparison.exists():
            b64 = image_to_base64(gene_comparison)
            if b64:
                html.append("<div class='figure-container'>")
                html.append(f"<img src='data:image/png;base64,{b64}' alt='Gene Ranking'>")
                html.append("<p class='figure-caption'>Figure 1: TE signal ranking by gene</p>")
                html.append("</div>")

    # Methods section
    html.append("<h2>Methods</h2>")
    html.append("<h3>Gene Selection</h3>")
    if data.get('gene_list') is not None and len(data['gene_list']) > 0:
        html.append(f"<p><strong>Total genes analyzed:</strong> {len(data['gene_list'])}</p>")
        html.append("<p><strong>Source:</strong> Canonical germ plasm genes from literature</p>")
        html.append("<p><strong>Tier 1 genes:</strong> nos, osk, pgc, gcl, vas, aub, CycB</p>")

    html.append("<h3>BLAST Parameters</h3>")
    html.append("<ul>")
    html.append("<li>E-value threshold: 10 (permissive)</li>")
    html.append("<li>Word size: 7 (short seeds)</li>")
    html.append("<li>Dust filter: Disabled</li>")
    html.append("<li>Soft masking: Disabled</li>")
    html.append("</ul>")

    # Results by gene (abbreviated for HTML)
    html.append("<h2>Results by Gene</h2>")
    if data.get('gene_results'):
        html.append("<table>")
        html.append("<tr><th>Gene</th><th>3'UTR Length</th><th>Hits</th><th>(+)</th><th>(-)</th><th>Max Signal</th><th>Clusters</th></tr>")
        for gene_data in data['gene_results']:
            plus_hits = gene_data.get('plus_strand_hits', '-')
            minus_hits = gene_data.get('minus_strand_hits', '-')
            html.append(f"<tr><td>{gene_data['gene']}</td><td>{gene_data['length']} bp</td>")
            html.append(f"<td>{gene_data['total_hits']}</td><td>{plus_hits}</td><td>{minus_hits}</td>")
            html.append(f"<td>{gene_data['max_signal']:.2e}</td>")
            html.append(f"<td>{gene_data.get('num_clusters', 0)}</td></tr>")
        html.append("</table>")

    # Candidate TE Fossils
    html.append("<h2>Candidate TE Fossils</h2>")
    if data.get('top_clusters'):
        html.append("<table>")
        html.append("<tr><th>Gene</th><th>Position</th><th>Length</th><th>Signal</th><th>Dominant TE</th></tr>")
        for cluster in data['top_clusters'][:10]:
            html.append(f"<tr><td>{cluster['gene']}</td><td>{cluster['start']}-{cluster['end']}</td>")
            html.append(f"<td>{cluster['length']} bp</td><td>{cluster['max_signal']:.2e}</td>")
            html.append(f"<td>{cluster.get('dominant_te', 'N/A')}</td></tr>")
        html.append("</table>")

    # Follow-up recommendations
    html.append("<h2>Follow-up Recommendations</h2>")
    html.append("<ol>")
    html.append("<li><strong>Secondary Structure:</strong> Run RNAfold on candidate regions</li>")
    html.append("<li><strong>Conservation:</strong> Check cross-species conservation</li>")
    html.append("<li><strong>Motif Overlap:</strong> Compare with known localization elements</li>")
    html.append("<li><strong>Experimental:</strong> Reporter assays for candidates</li>")
    html.append("</ol>")

    # Footer
    html.append("<hr>")
    html.append("<p><em>Report generated by TE Fossil Mining Pipeline</em></p>")
    html.append("</body></html>")

    # Write file
    with open(output_file, 'w') as f:
        f.write('\n'.join(html))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--density-dir',
        type=Path,
        help='Directory containing density_data.pkl'
    )
    parser.add_argument(
        '--clusters-file',
        type=Path,
        default=Path('data/te_clusters/clusters_summary.tsv'),
        help='Cluster summary file'
    )
    parser.add_argument(
        '--gene-list',
        type=Path,
        default=Path('data/gene_lists/germ_plasm_genes_consolidated.tsv'),
        help='Consolidated gene list'
    )
    parser.add_argument(
        '--figures-dir',
        type=Path,
        default=Path('figures'),
        help='Directory containing figure files'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('reports'),
        help='Output directory'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find density data if not specified
    if args.density_dir is None:
        sweep_dir = Path('results/parameter_sweep')
        if sweep_dir.exists():
            subdirs = sorted([d for d in sweep_dir.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                combo_dirs = sorted([d for d in subdirs[0].iterdir()
                                    if d.is_dir() and d.name.startswith('combo_')])
                if combo_dirs:
                    args.density_dir = combo_dirs[0] / 'density'

    print("TE Fossil Report Generation")
    print("=" * 60)

    # Collect data
    data = {
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }

    # Load gene list
    if args.gene_list.exists():
        data['gene_list'] = load_gene_list(args.gene_list)
        print(f"Loaded gene list: {len(data['gene_list'])} genes")

    # Load density data
    if args.density_dir and (args.density_dir / 'density_data.pkl').exists():
        density_dict = load_density_data(args.density_dir / 'density_data.pkl')
        if density_dict:
            print(f"Loaded density data: {len(density_dict)} sequences")

            # Extract gene-level results
            gene_results = []
            for qseqid, d in density_dict.items():
                gene = qseqid.split('_')[0]
                strand_counts = d.get('strand_counts', {'plus': 0, 'minus': 0})
                gene_results.append({
                    'gene': gene,
                    'qseqid': qseqid,
                    'length': d['length'],
                    'total_hits': d['total_hits'],
                    'plus_strand_hits': strand_counts.get('plus', 0),
                    'minus_strand_hits': strand_counts.get('minus', 0),
                    'max_signal': float(np.max(d['combined'])) if len(d['combined']) > 0 else 0
                })

            data['gene_results'] = gene_results

            # Top genes
            data['top_genes'] = sorted(gene_results, key=lambda x: -x['max_signal'])

    # Load cluster data
    if args.clusters_file.exists():
        clusters_df = load_cluster_data(args.clusters_file)
        if clusters_df is not None and len(clusters_df) > 0:
            print(f"Loaded cluster data: {len(clusters_df)} clusters")

            data['clusters_summary'] = {
                'total_clusters': len(clusters_df),
                'genes_with_clusters': clusters_df['qseqid'].nunique()
            }

            # Top clusters
            clusters_df['gene'] = clusters_df['qseqid'].apply(lambda x: x.split('_')[0])
            top_clusters = clusters_df.nlargest(20, 'max_signal').to_dict('records')
            data['top_clusters'] = top_clusters

            # Add cluster info to gene results
            if 'gene_results' in data:
                for gene_data in data['gene_results']:
                    gene_clusters = clusters_df[clusters_df['qseqid'] == gene_data['qseqid']]
                    gene_data['num_clusters'] = len(gene_clusters)
                    gene_data['clusters'] = gene_clusters.to_dict('records')

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate timestamp for filenames
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Generate reports
    print("\nGenerating reports...")

    md_file = args.output_dir / f'te_fossil_analysis_{timestamp}.md'
    generate_markdown_report(md_file, data)
    print(f"  Saved: {md_file.name}")

    html_file = args.output_dir / f'te_fossil_analysis_{timestamp}.html'
    generate_html_report(html_file, data, args.figures_dir)
    print(f"  Saved: {html_file.name}")

    print("\n" + "=" * 60)
    print("Report generation complete!")
    print(f"\nView HTML report at: {html_file}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
