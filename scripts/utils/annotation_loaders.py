"""
Annotation loading utilities for external data sources.

Parses:
- FlyBase RNA-Seq expression data (RPKM matrix)
- FlyBase GO annotations (GAF format)
- FlyBase gene groups and pathways
- FlyFISH RNA localization patterns
"""

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


def load_fbgn_to_symbol(path: Union[str, Path]) -> Dict[str, str]:
    """
    Load FBgn to symbol mapping from reference file.

    Args:
        path: Path to fbgn_to_symbol.tsv

    Returns:
        Dictionary mapping FBgn ID to gene symbol
    """
    path = Path(path)
    mapping = {}

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                fbgn = parts[0].strip()
                symbol = parts[1].strip()
                if fbgn.startswith('FBgn'):
                    mapping[fbgn] = symbol

    return mapping


def build_symbol_to_fbgn_map(path: Union[str, Path]) -> Dict[str, str]:
    """
    Build symbol to FBgn mapping (reverse of fbgn_to_symbol).

    Args:
        path: Path to fbgn_to_symbol.tsv

    Returns:
        Dictionary mapping gene symbol to FBgn ID
    """
    fbgn_to_symbol = load_fbgn_to_symbol(path)
    return {symbol: fbgn for fbgn, symbol in fbgn_to_symbol.items()}


def build_comprehensive_symbol_map(annotation_id_path: Union[str, Path]) -> Dict[str, str]:
    """
    Build comprehensive symbol/CG-number to FBgn mapping from FlyBase annotation ID file.

    This file contains mappings for:
    - Gene symbols -> FBgn
    - CG numbers (annotation IDs) -> FBgn
    - Secondary annotation IDs -> FBgn

    Args:
        annotation_id_path: Path to fbgn_annotation_ID.tsv

    Returns:
        Dictionary mapping various gene identifiers to FBgn ID
    """
    path = Path(annotation_id_path)
    mapping = {}

    with open(path, encoding='utf-8') as f:
        header_found = False
        for line in f:
            # Skip comment lines
            if line.startswith('##'):
                continue

            parts = line.strip().split('\t')

            # Skip header line
            if not header_found:
                header_found = True
                continue

            if len(parts) < 5:
                continue

            symbol = parts[0].strip()
            organism = parts[1].strip() if len(parts) > 1 else ''
            fbgn = parts[2].strip() if len(parts) > 2 else ''
            annotation_id = parts[4].strip() if len(parts) > 4 else ''
            secondary_ids = parts[5].strip() if len(parts) > 5 else ''

            # Only use Dmel genes
            if organism != 'Dmel':
                continue

            if not fbgn or not fbgn.startswith('FBgn'):
                continue

            # Map symbol (case-insensitive)
            if symbol:
                mapping[symbol] = fbgn
                mapping[symbol.lower()] = fbgn
                mapping[symbol.upper()] = fbgn

            # Map primary annotation ID (CG number)
            if annotation_id:
                mapping[annotation_id] = fbgn

            # Map secondary annotation IDs
            if secondary_ids:
                for sec_id in secondary_ids.split(','):
                    sec_id = sec_id.strip()
                    if sec_id:
                        mapping[sec_id] = fbgn

    return mapping


def load_rnaseq_expression(path: Union[str, Path]) -> Dict[str, Dict[str, float]]:
    """
    Load FlyBase RNA-Seq RPKM data.

    Handles the gene_rpkm_report format with columns:
    Release_ID, FBgn#, GeneSymbol, Parent_library_FBlc#, Parent_library_name,
    RNASource_FBlc#, RNASource_name, RPKM_value, ...

    Args:
        path: Path to gene_rpkm_report.tsv

    Returns:
        Dictionary mapping FBgn -> {tissue_name: RPKM_value}
    """
    path = Path(path)
    expression = defaultdict(dict)

    with open(path, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')

        header = None
        for row in reader:
            if not row or row[0].startswith('#'):
                continue

            # Find header row
            if header is None:
                if any('FBgn' in cell or 'RPKM' in cell.upper() for cell in row):
                    header = row
                    # Find column indices
                    fbgn_col = None
                    source_col = None
                    rpkm_col = None

                    for i, col in enumerate(header):
                        col_lower = col.lower()
                        if 'fbgn' in col_lower and fbgn_col is None:
                            fbgn_col = i
                        elif 'rnasource_name' in col_lower or 'source_name' in col_lower:
                            source_col = i
                        elif 'rpkm' in col_lower and rpkm_col is None:
                            rpkm_col = i

                    if fbgn_col is None or rpkm_col is None:
                        # Try alternative format
                        for i, col in enumerate(header):
                            if col.startswith('FBgn'):
                                fbgn_col = i
                            elif 'value' in col.lower() or 'rpkm' in col.lower():
                                rpkm_col = i
                continue

            if len(row) <= max(fbgn_col or 0, rpkm_col or 0):
                continue

            fbgn = row[fbgn_col] if fbgn_col is not None else None
            if not fbgn or not fbgn.startswith('FBgn'):
                continue

            # Get tissue/source name
            if source_col is not None and source_col < len(row):
                tissue = row[source_col]
            else:
                tissue = 'unknown'

            # Get RPKM value
            try:
                rpkm = float(row[rpkm_col]) if rpkm_col is not None else 0.0
            except (ValueError, IndexError):
                rpkm = 0.0

            # Clean tissue name for use as key
            tissue_key = _clean_tissue_name(tissue)
            expression[fbgn][tissue_key] = rpkm

    return dict(expression)


def load_rnaseq_matrix(path: Union[str, Path]) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    """
    Load FlyBase RNA-Seq RPKM matrix format.

    Matrix format: genes as rows, samples as columns.
    First columns are: gene_primary_id, gene_symbol, gene_fullname, gene_type
    Then tissue/sample columns follow.

    Args:
        path: Path to gene_rpkm_matrix.tsv

    Returns:
        Tuple of (expression_dict, tissue_list)
        - expression_dict: FBgn -> {tissue: RPKM}
        - tissue_list: List of tissue/sample names in order
    """
    path = Path(path)
    expression = {}
    tissues = []

    # Metadata columns to skip (not tissue data)
    metadata_cols = {'gene_primary_id', 'gene_symbol', 'gene_fullname', 'gene_type',
                     '#gene_primary_id'}

    with open(path, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')

        header = None
        fbgn_col = 0
        tissue_col_indices = []

        for row in reader:
            if not row or row[0].startswith('##'):
                continue

            # First non-comment row is header
            if header is None:
                header = row
                # Find FBgn column and identify tissue columns
                for i, col in enumerate(header):
                    col_clean = col.lstrip('#')
                    if col_clean in metadata_cols:
                        if 'primary_id' in col_clean or 'fbgn' in col_clean.lower():
                            fbgn_col = i
                        continue
                    # This is a tissue/sample column
                    tissue_col_indices.append(i)
                    tissues.append(_clean_tissue_name(col))
                continue

            if len(row) <= fbgn_col:
                continue

            fbgn = row[fbgn_col]
            if not fbgn.startswith('FBgn'):
                continue

            expression[fbgn] = {}
            for tissue_idx, col_idx in enumerate(tissue_col_indices):
                if col_idx < len(row) and tissue_idx < len(tissues):
                    try:
                        expression[fbgn][tissues[tissue_idx]] = float(row[col_idx])
                    except ValueError:
                        expression[fbgn][tissues[tissue_idx]] = 0.0

    return expression, tissues


def _clean_tissue_name(name: str) -> str:
    """
    Clean tissue/sample name for use as dictionary key.

    Handles FlyBase column names like:
    - mE_mRNA_A_MateF_4d_ovary_(FBlc0000208) -> adult_female_4d_ovary
    - RNA-Seq_Profile_FlyAtlas2_Adult_Female_Ovary_(FBlc0003627) -> adult_female_ovary
    - BCM_1_E2-4hr_(FBlc0000061) -> embryo_2_4hr
    """
    # Remove FlyBase library IDs
    name = re.sub(r'\(FBlc\d+\)', '', name).strip('_')

    # Handle FlyAtlas2 format
    if 'FlyAtlas2' in name:
        name = re.sub(r'RNA-Seq_Profile_FlyAtlas2_', '', name)
        name = name.replace('_', ' ').lower()
        name = re.sub(r'\s+', '_', name)
        return name

    # Handle modENCODE format (mE_mRNA_...)
    if name.startswith('mE_mRNA_'):
        name = name.replace('mE_mRNA_', '')
        # Parse tissue/stage info
        # A_MateF_4d_ovary -> adult_mated_female_4d_ovary
        name = name.replace('_A_', '_adult_')
        name = name.replace('_L3_', '_larva3_')
        name = name.replace('_L1_', '_larva1_')
        name = name.replace('_L2_', '_larva2_')
        name = name.replace('MateF', 'mated_female')
        name = name.replace('MateM', 'mated_male')
        name = name.replace('VirF', 'virgin_female')
        name = name.replace('AdF', 'adult_female')
        name = name.replace('AdM', 'adult_male')
        name = name.replace('_em', '_embryo_')
        name = name.replace('_dig_sys', '_digestive')
        name = name.replace('_acc_gland', '_accessory_gland')
        name = name.replace('_saliv', '_salivary')
        name = name.replace('_imag_disc', '_imaginal_disc')
        name = name.replace('_Wand_', '_wandering_')
        name = name.replace('_CNS', '_cns')
        name = name.replace('_WPP', '_white_prepupae')

    # Handle BCM format
    if name.startswith('BCM_'):
        name = name.replace('BCM_1_', '')
        name = name.replace('E2-4hr', 'embryo_2_4hr')
        name = name.replace('E14-16hr', 'embryo_14_16hr')
        name = name.replace('E2-16hr', 'embryo_2_16hr')
        name = name.replace('L3i', 'larva3')
        name = name.replace('P3d', 'pupae_3d')
        name = name.replace('FA3d', 'female_adult_3d')
        name = name.replace('MA3d', 'male_adult_3d')

    # Remove common prefixes
    name = re.sub(r'^BCM_HGSC_', '', name)
    name = re.sub(r'^modENCODE_', '', name)
    name = re.sub(r'^Knoblich_mRNA_', '', name)

    # Normalize separators
    name = name.lower()
    name = re.sub(r'[,;:\s]+', '_', name)
    name = re.sub(r'-+', '_', name)
    name = re.sub(r'_+', '_', name)
    name = name.strip('_')

    return name


def load_go_annotations(path: Union[str, Path]) -> Dict[str, List[Dict]]:
    """
    Load GO annotations from GAF (Gene Association Format) file.

    GAF 2.2 format columns:
    1. DB (e.g., FB)
    2. DB Object ID (e.g., FBgn0000003)
    3. DB Object Symbol
    4. Qualifier (NOT, enables, etc.)
    5. GO ID (e.g., GO:0003674)
    6. DB:Reference
    7. Evidence Code (IDA, IEA, etc.)
    8. With/From
    9. Aspect (F=molecular function, P=biological process, C=cellular component)
    10. DB Object Name
    11. DB Object Synonym
    12. DB Object Type
    13. Taxon
    14. Date
    15. Assigned By
    16+ Extensions, Annotations

    Args:
        path: Path to gene_association.gaf file

    Returns:
        Dictionary mapping FBgn -> list of GO annotation dicts
    """
    path = Path(path)
    annotations = defaultdict(list)

    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('!'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            db = parts[0]
            fbgn = parts[1]
            symbol = parts[2]
            qualifier = parts[3]
            go_id = parts[4]
            evidence = parts[6] if len(parts) > 6 else ''
            aspect = parts[8] if len(parts) > 8 else ''

            # Skip negative annotations (NOT)
            if 'NOT' in qualifier.upper():
                continue

            if not fbgn.startswith('FBgn'):
                continue

            annotations[fbgn].append({
                'go_id': go_id,
                'aspect': aspect,
                'evidence': evidence,
                'symbol': symbol,
                'qualifier': qualifier,
            })

    return dict(annotations)


def get_go_aspect_name(aspect: str) -> str:
    """Convert GO aspect code to full name."""
    aspects = {
        'F': 'molecular_function',
        'P': 'biological_process',
        'C': 'cellular_component',
    }
    return aspects.get(aspect.upper(), 'unknown')


def filter_go_by_aspect(
    annotations: Dict[str, List[Dict]],
    aspect: str
) -> Dict[str, Set[str]]:
    """
    Filter GO annotations to single aspect.

    Args:
        annotations: Full GO annotations dict
        aspect: 'F', 'P', or 'C'

    Returns:
        Dictionary mapping FBgn -> set of GO IDs for that aspect
    """
    filtered = defaultdict(set)

    for fbgn, annots in annotations.items():
        for annot in annots:
            if annot['aspect'].upper() == aspect.upper():
                filtered[fbgn].add(annot['go_id'])

    return dict(filtered)


def load_gene_groups(path: Union[str, Path]) -> Dict[str, List[str]]:
    """
    Load FlyBase gene group data.

    Format: FB_group_id, FB_group_symbol, FB_group_name, Parent_FB_group_id,
            Parent_FB_group_symbol, Group_member_FB_gene_id, Group_member_FB_gene_symbol

    Args:
        path: Path to gene_group_data.tsv

    Returns:
        Dictionary mapping group_symbol -> list of FBgn IDs
    """
    path = Path(path)
    groups = defaultdict(list)

    with open(path, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            if not row or row[0].startswith('#'):
                continue

            # Skip header
            if 'FB_group' in row[0]:
                continue

            if len(row) >= 7:
                group_symbol = row[1]
                member_fbgn = row[5]

                if member_fbgn.startswith('FBgn'):
                    groups[group_symbol].append(member_fbgn)

    return dict(groups)


def load_flyfish_localization(
    path: Union[str, Path],
    symbol_to_fbgn: Optional[Dict[str, str]] = None
) -> Dict[str, List[str]]:
    """
    Load FlyFISH subcellular localization data.

    Handles the FlyFISH annotation.csv format:
    - Probe, gene, FBgn ID, stage/tissue, term

    Note: The "FBgn ID" column contains clone IDs (FBcl), not gene IDs.
    We use the "gene" column (symbol) and map to FBgn.

    Args:
        path: Path to flyfish_localization.csv
        symbol_to_fbgn: Optional symbol->FBgn mapping for ID conversion

    Returns:
        Dictionary mapping FBgn (or symbol if no mapping) -> list of localization patterns
    """
    path = Path(path)
    localizations = defaultdict(set)

    with open(path, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)

        header = None
        gene_col = None
        term_col = None

        for row in reader:
            if not row:
                continue

            # Find header
            if header is None:
                header = [col.lower() for col in row]

                # Look for gene/symbol column
                for i, col in enumerate(header):
                    if col == 'gene' or col == 'symbol':
                        gene_col = i
                        break

                # Look for term/pattern column
                for i, col in enumerate(header):
                    if col == 'term' or col == 'pattern' or col == 'localization':
                        term_col = i
                        break

                continue

            if gene_col is None or term_col is None:
                continue

            if gene_col >= len(row) or term_col >= len(row):
                continue

            symbol = row[gene_col].strip()
            term = row[term_col].strip()

            if not symbol or not term:
                continue

            # Skip non-informative annotations
            if term.lower() in ('', 'na', 'n/a', 'none'):
                continue
            if term.startswith('Annotated by'):
                continue

            # Clean symbol - remove parenthetical annotations like "CG31764 (RE18101)"
            clean_symbol = re.sub(r'\s*\([^)]+\)', '', symbol).strip()

            # Convert to FBgn if mapping available
            gene_id = None
            if symbol_to_fbgn:
                # Try exact match first
                if symbol in symbol_to_fbgn:
                    gene_id = symbol_to_fbgn[symbol]
                # Try cleaned symbol
                elif clean_symbol in symbol_to_fbgn:
                    gene_id = symbol_to_fbgn[clean_symbol]
                # Try uppercase
                elif clean_symbol.upper() in symbol_to_fbgn:
                    gene_id = symbol_to_fbgn[clean_symbol.upper()]

            if gene_id is None:
                gene_id = clean_symbol

            localizations[gene_id].add(term.lower())

    # Convert sets to sorted lists
    return {gene: sorted(patterns) for gene, patterns in localizations.items()}


# Common localization patterns to look for
FLYFISH_LOCALIZATION_CATEGORIES = {
    'apical': ['apical', 'apico-lateral', 'apical enrichment'],
    'basal': ['basal', 'basal enrichment'],
    'posterior': ['posterior', 'posterior pole', 'pole plasm'],
    'anterior': ['anterior', 'anterior pole'],
    'nuclear': ['nuclear', 'nucleus', 'perinuclear'],
    'cytoplasmic': ['cytoplasmic', 'cytoplasm', 'uniform cytoplasmic'],
    'membrane': ['membrane', 'cortical', 'plasma membrane'],
    'granular': ['granular', 'punctate', 'foci'],
}


def categorize_flyfish_localization(patterns: List[str]) -> List[str]:
    """
    Map FlyFISH patterns to standardized categories.

    Args:
        patterns: List of raw localization patterns

    Returns:
        List of standardized category names
    """
    categories = set()

    for pattern in patterns:
        pattern_lower = pattern.lower()
        for category, keywords in FLYFISH_LOCALIZATION_CATEGORIES.items():
            if any(kw in pattern_lower for kw in keywords):
                categories.add(category)

    return sorted(categories)


def summarize_annotations(
    expression: Dict[str, Dict[str, float]],
    go_annotations: Dict[str, List[Dict]],
    flyfish: Dict[str, List[str]],
    gene_groups: Dict[str, List[str]]
) -> Dict[str, int]:
    """
    Generate summary statistics for loaded annotations.

    Returns:
        Dictionary of counts and statistics
    """
    all_genes = set()
    all_genes.update(expression.keys())
    all_genes.update(go_annotations.keys())
    all_genes.update(flyfish.keys())
    for members in gene_groups.values():
        all_genes.update(members)

    # Count GO terms by aspect
    go_counts = {'F': 0, 'P': 0, 'C': 0}
    for fbgn, annots in go_annotations.items():
        for annot in annots:
            aspect = annot.get('aspect', '').upper()
            if aspect in go_counts:
                go_counts[aspect] += 1

    return {
        'total_genes': len(all_genes),
        'genes_with_expression': len(expression),
        'genes_with_go': len(go_annotations),
        'genes_with_flyfish': len(flyfish),
        'go_mf_annotations': go_counts['F'],
        'go_bp_annotations': go_counts['P'],
        'go_cc_annotations': go_counts['C'],
        'gene_groups': len(gene_groups),
        'unique_tissues': len(set(
            tissue for tissues in expression.values() for tissue in tissues
        )),
    }
