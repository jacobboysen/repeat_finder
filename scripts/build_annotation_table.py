#!/usr/bin/env python3
"""
Build unified gene annotation table from external data sources.

Merges:
- FlyBase RNA-Seq expression data
- FlyBase GO annotations
- FlyFISH localization patterns
- FlyBase gene groups

Output: data/annotations/gene_annotations.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_references_dir
from utils.annotation_loaders import (
    load_fbgn_to_symbol,
    build_symbol_to_fbgn_map,
    build_comprehensive_symbol_map,
    load_rnaseq_expression,
    load_rnaseq_matrix,
    load_go_annotations,
    load_gene_groups,
    load_flyfish_localization,
    categorize_flyfish_localization,
    filter_go_by_aspect,
    summarize_annotations,
)


def get_all_genome_genes(references_dir: Path) -> set:
    """
    Get set of all genes in the genome from reference files.

    Uses fbgn_to_symbol.tsv as the authoritative gene list.
    """
    fbgn_to_symbol_path = references_dir / "fbgn_to_symbol.tsv"
    if fbgn_to_symbol_path.exists():
        mapping = load_fbgn_to_symbol(fbgn_to_symbol_path)
        return set(mapping.keys())

    # Fallback: parse gene_info.tsv
    gene_info_path = references_dir / "gene_info.tsv"
    if gene_info_path.exists():
        genes = set()
        with open(gene_info_path) as f:
            for line in f:
                if line.startswith('FBgn'):
                    genes.add(line.split('\t')[0])
        return genes

    return set()


def select_key_tissues(expression: dict, tissues: list = None) -> list:
    """
    Select key tissues for the annotation table.

    Returns list of (tissue_key, display_name) tuples.
    Uses substring matching to find tissues in the available list.
    """
    # Priority tissue patterns - biologically meaningful for TE/germ plasm analysis
    # Format: (patterns_to_match, display_name)
    priority_patterns = [
        # Germline - most important
        (['ovary'], 'ovary_rpkm'),
        (['testis'], 'testis_rpkm'),

        # Embryonic stages
        (['embryo_0_2', 'em0_2'], 'embryo_0_2h_rpkm'),
        (['embryo_2_4', 'em2_4'], 'embryo_2_4h_rpkm'),
        (['embryo_4_6', 'em4_6'], 'embryo_4_6h_rpkm'),
        (['embryo_6_8', 'em6_8'], 'embryo_6_8h_rpkm'),
        (['embryo_8_10', 'em8_10'], 'embryo_8_10h_rpkm'),
        (['embryo_10_12', 'em10_12'], 'embryo_10_12h_rpkm'),
        (['embryo_12_14', 'em12_14'], 'embryo_12_14h_rpkm'),
        (['embryo_14_16', 'em14_16'], 'embryo_14_16h_rpkm'),
        (['embryo_16_18', 'em16_18'], 'embryo_16_18h_rpkm'),
        (['embryo_18_20', 'em18_20'], 'embryo_18_20h_rpkm'),
        (['embryo_20_22', 'em20_22'], 'embryo_20_22h_rpkm'),
        (['embryo_22_24', 'em22_24'], 'embryo_22_24h_rpkm'),

        # Adult tissues
        (['_head', 'adult_head', 'female_head', 'male_head'], 'head_rpkm'),
        (['carcass'], 'carcass_rpkm'),
        (['digestive', 'dig_sys', 'midgut'], 'digestive_rpkm'),
        (['salivary', 'saliv'], 'salivary_rpkm'),
        (['accessory_gland', 'acc_gland'], 'accessory_gland_rpkm'),
        (['fat_body', 'fatbody', '_fat'], 'fat_body_rpkm'),

        # Larval stages
        (['l1_', '_l1'], 'l1_rpkm'),
        (['l2_', '_l2'], 'l2_rpkm'),
        (['l3_', '_l3', 'larva3'], 'l3_rpkm'),
        (['wandering', 'wand'], 'wandering_l3_rpkm'),

        # Neural
        (['cns', 'brain'], 'cns_rpkm'),
        (['neuron'], 'neuron_rpkm'),

        # Other stages
        (['pupae', 'pupa'], 'pupae_rpkm'),
        (['prepupae', 'wpp'], 'prepupae_rpkm'),
    ]

    # Get available tissues from data
    available = []
    if tissues:
        available = list(tissues)
    else:
        for gene_expr in expression.values():
            available.extend(gene_expr.keys())
        available = list(set(available))

    # Select tissues using pattern matching
    selected = []
    seen_names = set()
    used_tissues = set()

    for patterns, display_name in priority_patterns:
        if display_name in seen_names:
            continue

        # Find first tissue matching any pattern
        for tissue in available:
            if tissue in used_tissues:
                continue
            tissue_lower = tissue.lower()
            if any(pattern in tissue_lower for pattern in patterns):
                selected.append((tissue, display_name))
                seen_names.add(display_name)
                used_tissues.add(tissue)
                break

    # If we got very few, add any others
    if len(selected) < 10:
        for tissue in sorted(available):
            if tissue not in used_tissues:
                display_name = f"{tissue}_rpkm"
                if display_name not in seen_names:
                    selected.append((tissue, display_name))
                    seen_names.add(display_name)
                    used_tissues.add(tissue)
                    if len(selected) >= 25:
                        break

    return selected


def build_annotation_table(
    expression: dict,
    go_annotations: dict,
    flyfish: dict,
    gene_groups: dict,
    all_genes: set,
    fbgn_to_symbol: dict,
    key_tissues: list,
) -> pd.DataFrame:
    """
    Build unified annotation DataFrame.

    Args:
        expression: FBgn -> {tissue: RPKM}
        go_annotations: FBgn -> [GO annotations]
        flyfish: FBgn -> [localization patterns]
        gene_groups: group_name -> [FBgn list]
        all_genes: Set of all FBgn IDs
        fbgn_to_symbol: FBgn -> symbol mapping
        key_tissues: List of (tissue_key, display_name) tuples

    Returns:
        DataFrame with one row per gene
    """
    # Invert gene_groups to get per-gene membership
    gene_to_groups = {}
    for group_name, members in gene_groups.items():
        for fbgn in members:
            if fbgn not in gene_to_groups:
                gene_to_groups[fbgn] = []
            gene_to_groups[fbgn].append(group_name)

    # Build rows
    rows = []
    for fbgn in sorted(all_genes):
        row = {
            'fbgn': fbgn,
            'symbol': fbgn_to_symbol.get(fbgn, ''),
        }

        # Expression values
        gene_expr = expression.get(fbgn, {})
        for tissue_key, display_name in key_tissues:
            row[display_name] = gene_expr.get(tissue_key, 0.0)

        # GO annotations - collapsed by aspect
        gene_go = go_annotations.get(fbgn, [])
        go_mf = set()
        go_bp = set()
        go_cc = set()
        for annot in gene_go:
            aspect = annot.get('aspect', '').upper()
            go_id = annot.get('go_id', '')
            if aspect == 'F':
                go_mf.add(go_id)
            elif aspect == 'P':
                go_bp.add(go_id)
            elif aspect == 'C':
                go_cc.add(go_id)

        row['go_mf'] = '|'.join(sorted(go_mf)) if go_mf else ''
        row['go_bp'] = '|'.join(sorted(go_bp)) if go_bp else ''
        row['go_cc'] = '|'.join(sorted(go_cc)) if go_cc else ''

        # FlyFISH localization
        gene_flyfish = flyfish.get(fbgn, [])
        row['flyfish_patterns'] = '|'.join(gene_flyfish) if gene_flyfish else ''
        categories = categorize_flyfish_localization(gene_flyfish)
        row['flyfish_categories'] = '|'.join(categories) if categories else ''

        # Gene groups
        groups = gene_to_groups.get(fbgn, [])
        row['gene_groups'] = '|'.join(sorted(groups)) if groups else ''

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--annotations-dir',
        type=Path,
        default=None,
        help='Directory containing downloaded annotations (default: data/annotations/raw)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output path (default: data/annotations/gene_annotations.tsv)'
    )
    parser.add_argument(
        '--all-tissues',
        action='store_true',
        help='Include all available tissues (not just key ones)'
    )

    args = parser.parse_args()

    project_root = get_project_root()
    references_dir = get_references_dir()

    # Set paths
    if args.annotations_dir:
        annotations_dir = args.annotations_dir
    else:
        annotations_dir = project_root / "data" / "annotations" / "raw"

    if args.output:
        output_path = args.output
    else:
        output_path = project_root / "data" / "annotations" / "gene_annotations.tsv"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check for required files
    required_files = {
        'expression': annotations_dir / "gene_rpkm_report.tsv",
        'expression_matrix': annotations_dir / "gene_rpkm_matrix.tsv",
        'go': annotations_dir / "gene_association.gaf",
        'groups': annotations_dir / "gene_group_data.tsv",
        'flyfish': annotations_dir / "flyfish_localization.csv",
    }

    print("Loading annotation data...")
    print(f"  Source directory: {annotations_dir}")
    print()

    # Load FBgn to symbol mapping
    fbgn_symbol_path = references_dir / "fbgn_to_symbol.tsv"
    if fbgn_symbol_path.exists():
        fbgn_to_symbol = load_fbgn_to_symbol(fbgn_symbol_path)
        symbol_to_fbgn = build_symbol_to_fbgn_map(fbgn_symbol_path)
        print(f"  Loaded {len(fbgn_to_symbol):,} gene ID mappings")
    else:
        fbgn_to_symbol = {}
        symbol_to_fbgn = {}
        print("  Warning: fbgn_to_symbol.tsv not found")

    # Get all genome genes
    all_genes = get_all_genome_genes(references_dir)
    print(f"  Found {len(all_genes):,} genes in genome")

    # Load expression data
    expression = {}
    tissues = []
    if required_files['expression_matrix'].exists():
        print(f"  Loading expression matrix...")
        expression, tissues = load_rnaseq_matrix(required_files['expression_matrix'])
        print(f"    {len(expression):,} genes, {len(tissues)} tissues")
    elif required_files['expression'].exists():
        print(f"  Loading expression report...")
        expression = load_rnaseq_expression(required_files['expression'])
        print(f"    {len(expression):,} genes")
    else:
        print("  Warning: No expression data found")

    # Load GO annotations
    if required_files['go'].exists():
        print(f"  Loading GO annotations...")
        go_annotations = load_go_annotations(required_files['go'])
        total_annots = sum(len(v) for v in go_annotations.values())
        print(f"    {len(go_annotations):,} genes, {total_annots:,} annotations")
    else:
        go_annotations = {}
        print("  Warning: GO annotations not found")

    # Load gene groups
    if required_files['groups'].exists():
        print(f"  Loading gene groups...")
        gene_groups = load_gene_groups(required_files['groups'])
        total_members = sum(len(v) for v in gene_groups.values())
        print(f"    {len(gene_groups):,} groups, {total_members:,} memberships")
    else:
        gene_groups = {}
        print("  Warning: Gene groups not found")

    # Load FlyFISH with comprehensive symbol mapping
    if required_files['flyfish'].exists():
        print(f"  Loading FlyFISH data...")
        # Try to use comprehensive mapping from FlyBase annotation ID file
        annotation_id_path = annotations_dir / "fbgn_annotation_ID.tsv"
        if annotation_id_path.exists():
            print("    Using comprehensive FBgn mapping...")
            comprehensive_map = build_comprehensive_symbol_map(annotation_id_path)
            flyfish = load_flyfish_localization(required_files['flyfish'], comprehensive_map)
        else:
            flyfish = load_flyfish_localization(required_files['flyfish'], symbol_to_fbgn)
        print(f"    {len(flyfish):,} genes with localization data")
    else:
        flyfish = {}
        print("  Warning: FlyFISH data not found")

    # Update all_genes with any genes found in annotations
    all_genes.update(expression.keys())
    all_genes.update(go_annotations.keys())
    # Only add FlyFISH genes that have proper FBgn IDs
    all_genes.update(k for k in flyfish.keys() if k.startswith('FBgn'))
    for members in gene_groups.values():
        all_genes.update(members)

    print()
    print(f"Total genes to annotate: {len(all_genes):,}")

    # Select key tissues
    if args.all_tissues:
        key_tissues = [(t, f"{t}_rpkm") for t in sorted(tissues if tissues else
                       set(t for g in expression.values() for t in g.keys()))]
    else:
        key_tissues = select_key_tissues(expression, tissues)

    print(f"Including {len(key_tissues)} tissue columns")

    # Build table
    print()
    print("Building annotation table...")
    df = build_annotation_table(
        expression=expression,
        go_annotations=go_annotations,
        flyfish=flyfish,
        gene_groups=gene_groups,
        all_genes=all_genes,
        fbgn_to_symbol=fbgn_to_symbol,
        key_tissues=key_tissues,
    )

    # Save
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved to: {output_path}")
    print(f"  {len(df):,} genes")
    print(f"  {len(df.columns)} columns")

    # Summary statistics
    print()
    print("Summary:")
    print(f"  Genes with expression: {(df[[c for c in df.columns if '_rpkm' in c]].sum(axis=1) > 0).sum():,}")
    print(f"  Genes with GO MF: {(df['go_mf'] != '').sum():,}")
    print(f"  Genes with GO BP: {(df['go_bp'] != '').sum():,}")
    print(f"  Genes with GO CC: {(df['go_cc'] != '').sum():,}")
    print(f"  Genes with FlyFISH: {(df['flyfish_patterns'] != '').sum():,}")
    print(f"  Genes in groups: {(df['gene_groups'] != '').sum():,}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
