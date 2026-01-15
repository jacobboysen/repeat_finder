#!/usr/bin/env python3
"""
Build functional gene sets from unified annotation table.

Creates gene sets based on:
- Tissue expression (high expression, tissue-enriched)
- GO term categories (molecular function, biological process, cellular component)
- FlyFISH localization patterns

Output: data/gene_lists/functional/
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root


# Expression-based gene set definitions
EXPRESSION_SETS = {
    # High expression thresholds (absolute RPKM)
    'high_threshold': 10.0,
    'medium_threshold': 1.0,

    # Enrichment ratio (vs mean across tissues)
    'enrichment_ratio': 2.0,

    # Ubiquitous expression
    'ubiquitous_min_tissues': 5,
    'ubiquitous_min_rpkm': 5.0,
}

# Key tissues to generate sets for (patterns to match in column names)
TISSUE_PATTERNS = {
    'ovary': ['ovary', 'ovar'],
    'testis': ['testis', 'teste'],
    'embryo_early': ['embryo_0', 'embryo_2', 'em0_2', 'em2_4'],
    'embryo_mid': ['embryo_4', 'embryo_6', 'embryo_8', 'em4_6', 'em6_8', 'em8_10'],
    'embryo_late': ['embryo_10', 'embryo_12', 'embryo_14', 'embryo_16', 'embryo_18', 'embryo_20'],
    'head': ['head'],
    'carcass': ['carcass'],
    'larva': ['larva', 'l3_', '_l3', 'l1_', 'l2_'],
    'fat_body': ['fat', 'fatbody'],
    'digestive': ['digestive', 'dig_', 'midgut', 'hindgut', 'crop'],
    'cns': ['cns', 'brain', 'neuron', 'neuroblast'],
    'salivary': ['saliv', 'salivary'],
    'accessory_gland': ['accessory', 'acc_gland'],
}

# GO term categories (curated for biological relevance)
GO_CATEGORIES = {
    # Molecular Function
    'rna_binding': ['GO:0003723', 'GO:0000166'],  # RNA binding, nucleotide binding
    'dna_binding': ['GO:0003677', 'GO:0003700', 'GO:0043565'],  # DNA binding, TF activity
    'kinase': ['GO:0016301', 'GO:0004672', 'GO:0004674'],  # kinase activity
    'phosphatase': ['GO:0016791', 'GO:0004721'],  # phosphatase activity
    'helicase': ['GO:0004386', 'GO:0008026'],  # helicase activity
    'protease': ['GO:0008233', 'GO:0004175'],  # peptidase activity

    # Biological Process
    'transcription': ['GO:0006355', 'GO:0006351', 'GO:0045944'],  # transcription regulation
    'translation': ['GO:0006412', 'GO:0006414'],  # translation
    'rna_processing': ['GO:0006396', 'GO:0006397', 'GO:0008380'],  # RNA processing, splicing
    'cell_cycle': ['GO:0007049', 'GO:0000278', 'GO:0051301'],  # cell cycle, mitosis
    'dna_repair': ['GO:0006281', 'GO:0006974'],  # DNA repair
    'transport': ['GO:0006810', 'GO:0015031'],  # transport
    'signaling': ['GO:0007165', 'GO:0023052'],  # signal transduction
    'development': ['GO:0007275', 'GO:0009790'],  # development
    'apoptosis': ['GO:0006915', 'GO:0008219'],  # apoptosis
    'chromatin': ['GO:0006325', 'GO:0016568'],  # chromatin organization

    # Cellular Component
    'nucleus': ['GO:0005634', 'GO:0005654'],  # nucleus, nucleoplasm
    'cytoplasm': ['GO:0005737', 'GO:0005829'],  # cytoplasm, cytosol
    'mitochondria': ['GO:0005739', 'GO:0005743'],  # mitochondrion
    'er': ['GO:0005783', 'GO:0005789'],  # ER
    'golgi': ['GO:0005794'],  # Golgi
    'membrane': ['GO:0016020', 'GO:0005886'],  # membrane, plasma membrane
    'ribosome': ['GO:0005840', 'GO:0022626'],  # ribosome
    'centrosome': ['GO:0005813', 'GO:0000793'],  # centrosome
}

# FlyFISH localization categories (based on actual annotation terms)
FLYFISH_CATEGORIES = {
    # Expression timing
    'maternal': ['maternal'],
    'zygotic': ['zygotic'],
    'ubiquitous': ['ubiquitous'],
    # Subcellular localization
    'localized': ['localized'],
    'cytoplasmic_foci': ['cytoplasmic foci', 'foci'],
    'perinuclear': ['perinuclear'],
    'apical': ['apical enrichment', 'apical localization', 'apical exclusion'],
    'basal': ['basal enrichment', 'basal localization'],
    'membrane': ['membrane associated'],
    # Germ plasm / pole cells
    'pole_cell': ['pole cell localization', 'pole cell enrichment', 'pole plasm'],
    'yolk_plasm': ['yolk plasm localization', 'yolk plasm enrichment', 'yolk cortex'],
    'posterior': ['posterior localization'],
    # Tissue-specific (from embryo stages)
    'blastoderm_nuclei': ['blastoderm nuclei'],
    'ectoderm': ['ectoderm'],
    'mesoderm': ['mesoderm'],
    'nervous_system': ['nervous system', 'cns', 'brain', 'ventral nerve cord'],
}


def get_expression_columns(df: pd.DataFrame) -> list:
    """Get all expression (RPKM) columns from dataframe."""
    return [col for col in df.columns if col.endswith('_rpkm')]


def find_matching_columns(df: pd.DataFrame, patterns: list) -> list:
    """Find all expression columns matching any of the patterns."""
    expr_cols = get_expression_columns(df)
    matches = []
    for col in expr_cols:
        col_lower = col.lower()
        if any(pattern in col_lower for pattern in patterns):
            matches.append(col)
    return matches


def find_matching_column(df: pd.DataFrame, tissue: str) -> str:
    """Find first expression column matching tissue name."""
    candidates = find_matching_columns(df, [tissue.lower()])
    if candidates:
        return candidates[0]
    return None


def build_expression_sets(df: pd.DataFrame) -> dict:
    """
    Build gene sets based on expression thresholds.

    Returns:
        Dictionary mapping set_name -> list of FBgn IDs
    """
    gene_sets = {}
    expr_cols = get_expression_columns(df)

    if not expr_cols:
        print("  Warning: No expression columns found")
        return gene_sets

    print(f"  Found {len(expr_cols)} expression columns")

    # Calculate mean expression per gene (across all tissues)
    df['mean_expr'] = df[expr_cols].mean(axis=1)

    # Per-tissue sets using pattern matching
    for tissue_name, patterns in TISSUE_PATTERNS.items():
        matching_cols = find_matching_columns(df, patterns)
        if not matching_cols:
            continue

        print(f"    {tissue_name}: {len(matching_cols)} columns matched")

        # Take mean across all matching columns for this tissue
        tissue_expr = df[matching_cols].mean(axis=1)

        # High expression (absolute threshold)
        high_mask = tissue_expr >= EXPRESSION_SETS['high_threshold']
        set_name = f"expr_{tissue_name}_high"
        gene_sets[set_name] = df.loc[high_mask, 'fbgn'].tolist()

        # Tissue-enriched (relative to mean across all tissues)
        enriched_mask = (tissue_expr >= EXPRESSION_SETS['medium_threshold']) & \
                        (tissue_expr >= EXPRESSION_SETS['enrichment_ratio'] * df['mean_expr'])
        set_name = f"expr_{tissue_name}_enriched"
        gene_sets[set_name] = df.loc[enriched_mask, 'fbgn'].tolist()

    # Ubiquitous expression
    expr_matrix = df[expr_cols].values
    n_expressed = (expr_matrix >= EXPRESSION_SETS['ubiquitous_min_rpkm']).sum(axis=1)
    ubiq_mask = n_expressed >= EXPRESSION_SETS['ubiquitous_min_tissues']
    gene_sets['expr_ubiquitous'] = df.loc[ubiq_mask, 'fbgn'].tolist()

    # Low/no expression
    max_expr = df[expr_cols].max(axis=1)
    low_mask = max_expr < EXPRESSION_SETS['medium_threshold']
    gene_sets['expr_low'] = df.loc[low_mask, 'fbgn'].tolist()

    # Germline-specific (high in ovary/testis, low in somatic tissues)
    ovary_cols = find_matching_columns(df, ['ovary'])
    testis_cols = find_matching_columns(df, ['testis'])
    head_cols = find_matching_columns(df, ['head'])
    carcass_cols = find_matching_columns(df, ['carcass'])

    somatic_cols = head_cols + carcass_cols
    if ovary_cols and somatic_cols:
        ovary_expr = df[ovary_cols].mean(axis=1)
        somatic_expr = df[somatic_cols].mean(axis=1)
        germline_f_mask = (ovary_expr >= EXPRESSION_SETS['high_threshold']) & \
                          (ovary_expr >= 5 * (somatic_expr + 0.1))
        gene_sets['expr_ovary_specific'] = df.loc[germline_f_mask, 'fbgn'].tolist()

    if testis_cols and somatic_cols:
        testis_expr = df[testis_cols].mean(axis=1)
        somatic_expr = df[somatic_cols].mean(axis=1)
        germline_m_mask = (testis_expr >= EXPRESSION_SETS['high_threshold']) & \
                          (testis_expr >= 5 * (somatic_expr + 0.1))
        gene_sets['expr_testis_specific'] = df.loc[germline_m_mask, 'fbgn'].tolist()

    # Maternal (high in early embryo, low in later stages)
    early_cols = find_matching_columns(df, ['embryo_0', 'embryo_2', 'em0_2', 'em2_4'])
    late_cols = find_matching_columns(df, ['larva', 'l3_', 'adult'])
    if early_cols and late_cols:
        early_mean = df[early_cols].mean(axis=1)
        late_mean = df[late_cols].mean(axis=1)
        maternal_mask = (early_mean >= EXPRESSION_SETS['high_threshold']) & \
                        (early_mean >= 3 * (late_mean + 0.1))
        gene_sets['expr_maternal'] = df.loc[maternal_mask, 'fbgn'].tolist()

    return gene_sets


def build_go_sets(df: pd.DataFrame) -> dict:
    """
    Build gene sets based on GO term categories.

    Returns:
        Dictionary mapping set_name -> list of FBgn IDs
    """
    gene_sets = {}

    # Check which GO columns exist
    go_cols = {'mf': 'go_mf', 'bp': 'go_bp', 'cc': 'go_cc'}
    available_cols = {k: v for k, v in go_cols.items() if v in df.columns}

    if not available_cols:
        print("  Warning: No GO columns found")
        return gene_sets

    for category, go_ids in GO_CATEGORIES.items():
        matching_genes = set()

        for col in available_cols.values():
            for idx, row in df.iterrows():
                go_str = str(row[col]) if pd.notna(row[col]) else ''
                if any(go_id in go_str for go_id in go_ids):
                    matching_genes.add(row['fbgn'])

        set_name = f"go_{category}"
        gene_sets[set_name] = list(matching_genes)

    return gene_sets


def build_flyfish_sets(df: pd.DataFrame) -> dict:
    """
    Build gene sets based on FlyFISH localization.

    Returns:
        Dictionary mapping set_name -> list of FBgn IDs
    """
    gene_sets = {}

    # Check for FlyFISH columns - prefer flyfish_patterns for raw annotation terms
    flyfish_col = None
    for col in ['flyfish_patterns', 'flyfish_categories']:
        if col in df.columns:
            flyfish_col = col
            break

    if flyfish_col is None:
        print("  Warning: No FlyFISH columns found")
        return gene_sets

    for category, patterns in FLYFISH_CATEGORIES.items():
        matching_genes = []

        for idx, row in df.iterrows():
            loc_str = str(row[flyfish_col]).lower() if pd.notna(row[flyfish_col]) else ''
            if any(pattern in loc_str for pattern in patterns):
                matching_genes.append(row['fbgn'])

        set_name = f"flyfish_{category}"
        gene_sets[set_name] = matching_genes

    # Any FlyFISH annotation
    has_flyfish = df[flyfish_col].notna() & (df[flyfish_col] != '')
    gene_sets['flyfish_any'] = df.loc[has_flyfish, 'fbgn'].tolist()

    return gene_sets


def build_gene_group_sets(df: pd.DataFrame) -> dict:
    """
    Build gene sets from FlyBase gene groups.

    Returns:
        Dictionary mapping set_name -> list of FBgn IDs
    """
    gene_sets = {}

    if 'gene_groups' not in df.columns:
        return gene_sets

    # Collect all unique groups
    all_groups = defaultdict(list)
    for idx, row in df.iterrows():
        groups_str = str(row['gene_groups']) if pd.notna(row['gene_groups']) else ''
        if groups_str:
            for group in groups_str.split('|'):
                group = group.strip()
                if group:
                    all_groups[group].append(row['fbgn'])

    # Filter to groups with meaningful size
    for group, genes in all_groups.items():
        if 10 <= len(genes) <= 2000:  # Reasonable group size
            # Clean group name for filename
            clean_name = group.lower().replace(' ', '_').replace('-', '_')
            clean_name = ''.join(c for c in clean_name if c.isalnum() or c == '_')
            set_name = f"group_{clean_name}"
            gene_sets[set_name] = genes

    return gene_sets


def save_gene_set(genes: list, output_path: Path, set_name: str):
    """Save gene set to file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        for fbgn in sorted(set(genes)):  # Deduplicate and sort
            f.write(f"{fbgn}\n")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--annotations',
        type=Path,
        default=None,
        help='Path to gene_annotations.tsv (default: data/annotations/gene_annotations.tsv)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        help='Output directory for gene sets (default: data/gene_lists/functional)'
    )
    parser.add_argument(
        '--min-genes',
        type=int,
        default=10,
        help='Minimum genes per set to include (default: 10)'
    )
    parser.add_argument(
        '--max-genes',
        type=int,
        default=5000,
        help='Maximum genes per set to include (default: 5000)'
    )
    parser.add_argument(
        '--include-groups',
        action='store_true',
        help='Include FlyBase gene group sets (may create many files)'
    )

    args = parser.parse_args()

    project_root = get_project_root()

    # Set paths
    if args.annotations:
        annotations_path = args.annotations
    else:
        annotations_path = project_root / "data" / "annotations" / "gene_annotations.tsv"

    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = project_root / "data" / "gene_lists" / "functional"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load annotations
    if not annotations_path.exists():
        print(f"Error: Annotation file not found: {annotations_path}")
        print("Run build_annotation_table.py first.")
        return 1

    print(f"Loading annotations from: {annotations_path}")
    df = pd.read_csv(annotations_path, sep='\t')
    print(f"  {len(df):,} genes loaded")
    print()

    # Build gene sets
    all_sets = {}

    print("Building expression-based gene sets...")
    expr_sets = build_expression_sets(df)
    all_sets.update(expr_sets)
    print(f"  Created {len(expr_sets)} sets")

    print("Building GO-based gene sets...")
    go_sets = build_go_sets(df)
    all_sets.update(go_sets)
    print(f"  Created {len(go_sets)} sets")

    print("Building FlyFISH-based gene sets...")
    flyfish_sets = build_flyfish_sets(df)
    all_sets.update(flyfish_sets)
    print(f"  Created {len(flyfish_sets)} sets")

    if args.include_groups:
        print("Building gene group sets...")
        group_sets = build_gene_group_sets(df)
        all_sets.update(group_sets)
        print(f"  Created {len(group_sets)} sets")

    # Filter by size and save
    print()
    print(f"Saving gene sets to: {output_dir}")

    saved_count = 0
    skipped_small = 0
    skipped_large = 0

    for set_name, genes in sorted(all_sets.items()):
        n_genes = len(set(genes))

        if n_genes < args.min_genes:
            skipped_small += 1
            continue
        if n_genes > args.max_genes:
            skipped_large += 1
            continue

        output_path = output_dir / f"{set_name}_fbgn_ids.txt"
        save_gene_set(genes, output_path, set_name)
        saved_count += 1

    print()
    print(f"Saved {saved_count} gene sets")
    if skipped_small:
        print(f"Skipped {skipped_small} sets with < {args.min_genes} genes")
    if skipped_large:
        print(f"Skipped {skipped_large} sets with > {args.max_genes} genes")

    # Create summary file
    summary_path = output_dir / "gene_sets_summary.tsv"
    with open(summary_path, 'w') as f:
        f.write("set_name\tn_genes\tcategory\n")
        for set_name, genes in sorted(all_sets.items()):
            n_genes = len(set(genes))
            if args.min_genes <= n_genes <= args.max_genes:
                if set_name.startswith('expr_'):
                    category = 'expression'
                elif set_name.startswith('go_'):
                    category = 'go'
                elif set_name.startswith('flyfish_'):
                    category = 'flyfish'
                elif set_name.startswith('group_'):
                    category = 'gene_group'
                else:
                    category = 'other'
                f.write(f"{set_name}\t{n_genes}\t{category}\n")

    print(f"Summary saved to: {summary_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
