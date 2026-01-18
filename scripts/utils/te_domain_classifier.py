"""
TE Domain Classifier - Annotate BLAST hits by TE internal domain type.

Uses position-based heuristics to classify which region of a TE element
(gag, pol, env, LTR, transposase, etc.) a BLAST hit maps to.

Domain boundaries are approximate and based on typical TE structures:
- LTR retrotransposons: 5'LTR - gag - pol - (env) - 3'LTR
- LINE elements: 5'UTR - ORF1 - ORF2 (RT + endonuclease)
- DNA transposons: TIR - transposase - TIR

Note: These are heuristic approximations. For precise domain mapping,
use Dfam HMMs or curated per-family annotations.
"""

from typing import Dict, Optional, Tuple

# Domain boundaries as fractions of total element length
# Based on typical Drosophila TE structures

LTR_RETROTRANSPOSON_DOMAINS = {
    # Standard LTR structure: 5'LTR - gag - pol - 3'LTR
    # (env is often missing or degenerate in Drosophila elements)
    '5_ltr': (0.0, 0.12),      # First ~12% is 5' LTR
    'gag': (0.12, 0.35),       # Gag capsid protein ~12-35%
    'pol': (0.35, 0.88),       # Pol (RT, RNaseH, INT) ~35-88%
    '3_ltr': (0.88, 1.0),      # Last ~12% is 3' LTR
}

# More detailed pol subdomains (optional refinement)
POL_SUBDOMAINS = {
    'protease': (0.35, 0.42),     # Protease
    'rt': (0.42, 0.62),           # Reverse transcriptase
    'rnaseH': (0.62, 0.72),       # RNase H
    'integrase': (0.72, 0.88),    # Integrase
}

LINE_DOMAINS = {
    # Non-LTR retrotransposons
    '5_utr': (0.0, 0.08),         # 5' UTR region
    'orf1': (0.08, 0.35),         # ORF1 (RNA binding, chaperone)
    'orf2_rt': (0.35, 0.70),      # ORF2 RT domain
    'orf2_en': (0.70, 0.92),      # ORF2 endonuclease
    '3_utr': (0.92, 1.0),         # 3' UTR/poly-A
}

DNA_TRANSPOSON_DOMAINS = {
    # DNA transposons with terminal inverted repeats
    '5_tir': (0.0, 0.10),         # 5' Terminal Inverted Repeat
    'transposase': (0.10, 0.90),  # Central transposase coding region
    '3_tir': (0.90, 1.0),         # 3' Terminal Inverted Repeat
}

HELITRON_DOMAINS = {
    # Rolling-circle transposons
    '5_end': (0.0, 0.05),         # 5' TC motif
    'rep_helicase': (0.05, 0.60), # Rep/helicase domain
    'internal': (0.60, 0.95),     # Variable internal region
    '3_hairpin': (0.95, 1.0),     # 3' hairpin structure
}

# Map TE class to domain structure
CLASS_TO_DOMAINS = {
    'LTR': LTR_RETROTRANSPOSON_DOMAINS,
    'LINE': LINE_DOMAINS,
    'DNA': DNA_TRANSPOSON_DOMAINS,
    'Helitron': HELITRON_DOMAINS,
}

# LTR family patterns for more specific classification
LTR_FAMILIES = {
    'gypsy', 'copia', 'bel', 'pao', 'mdg', 'roo', '412', '297',
    'blood', 'accord', 'tirant', 'springer', 'opus', 'diver',
    'quasimodo', 'idefix', 'invader', 'gtwin', 'tabor', 'stalker',
    'max', 'transpac', 'micropia', 'blastopia', '17.6', '1731'
}

LINE_FAMILIES = {
    'jockey', 'doc', 'i-element', 'f-element', 'g-element',
    'x-element', 'het-a', 'tart', 'tahre', 'r1', 'r2', 'cr1',
    'rt1', 'baggins', 'juan', 'ivk', 'bs', 'fw'
}

DNA_FAMILIES = {
    'p-element', 'hobo', 'pogo', 'bari', 's-element',
    'transib', 'tc1', 'mariner', 'piggybac', 'mite',
    'protop', 'dnarep', 'looper', '1360'
}

HELITRON_FAMILIES = {'helitron', 'dine'}


def infer_te_class(te_name: str) -> str:
    """
    Infer TE class from name/family.

    Args:
        te_name: TE family name or ID

    Returns:
        TE class: 'LTR', 'LINE', 'DNA', 'Helitron', or 'Unknown'
    """
    name_lower = te_name.lower()

    # Check explicit class in name (e.g., "Gypsy-2_DM#LTR/Gypsy")
    if '#' in name_lower:
        class_part = name_lower.split('#')[1].split('/')[0]
        if 'ltr' in class_part:
            return 'LTR'
        elif 'line' in class_part:
            return 'LINE'
        elif 'dna' in class_part:
            return 'DNA'

    # Check family patterns
    for family in LTR_FAMILIES:
        if family in name_lower:
            return 'LTR'

    for family in LINE_FAMILIES:
        if family in name_lower:
            return 'LINE'

    for family in DNA_FAMILIES:
        if family in name_lower:
            return 'DNA'

    for family in HELITRON_FAMILIES:
        if family in name_lower:
            return 'Helitron'

    return 'Unknown'


def get_relative_position(sstart: int, send: int, te_length: int) -> Tuple[float, float]:
    """
    Calculate relative position within TE as fraction (0-1).

    Args:
        sstart: Subject (TE) alignment start position
        send: Subject (TE) alignment end position
        te_length: Total length of TE element

    Returns:
        Tuple of (start_fraction, end_fraction)
    """
    # Normalize to 0-1 range
    # Handle both orientations (sstart < send or sstart > send)
    pos1 = min(sstart, send) / te_length
    pos2 = max(sstart, send) / te_length

    # Clamp to valid range
    pos1 = max(0.0, min(1.0, pos1))
    pos2 = max(0.0, min(1.0, pos2))

    return (pos1, pos2)


def classify_te_domain(
    te_id: str,
    sstart: int,
    send: int,
    te_length: int,
    te_class: Optional[str] = None
) -> Dict[str, str]:
    """
    Classify which TE domain a BLAST hit maps to.

    Args:
        te_id: TE identifier (used to infer class if not provided)
        sstart: Subject alignment start position
        send: Subject alignment end position
        te_length: Total length of TE element
        te_class: Optional pre-computed TE class

    Returns:
        Dict with keys:
            - 'te_class': LTR, LINE, DNA, Helitron, Unknown
            - 'domain': Primary domain hit (gag, pol, ltr, transposase, etc.)
            - 'domain_detail': More specific subdomain if applicable
            - 'position': 'start', 'middle', 'end' based on position
    """
    # Infer class if not provided
    if te_class is None:
        te_class = infer_te_class(te_id)

    # Get relative position
    rel_start, rel_end = get_relative_position(sstart, send, te_length)
    midpoint = (rel_start + rel_end) / 2

    # Get domain map for this class
    domains = CLASS_TO_DOMAINS.get(te_class, {})

    # Find which domain(s) the hit overlaps
    primary_domain = 'unknown'
    domain_detail = None

    for domain_name, (dom_start, dom_end) in domains.items():
        # Check if hit midpoint falls in this domain
        if dom_start <= midpoint <= dom_end:
            primary_domain = domain_name
            break

    # For LTR elements, check pol subdomains
    if te_class == 'LTR' and primary_domain == 'pol':
        for subdomain_name, (sub_start, sub_end) in POL_SUBDOMAINS.items():
            if sub_start <= midpoint <= sub_end:
                domain_detail = subdomain_name
                break

    # General position classification
    if midpoint < 0.2:
        position = 'start'
    elif midpoint > 0.8:
        position = 'end'
    else:
        position = 'middle'

    return {
        'te_class': te_class,
        'domain': primary_domain,
        'domain_detail': domain_detail,
        'position': position,
        'rel_start': round(rel_start, 3),
        'rel_end': round(rel_end, 3),
    }


def get_domain_description(domain: str) -> str:
    """Get human-readable description of a TE domain."""
    descriptions = {
        '5_ltr': "5' Long Terminal Repeat (regulatory)",
        '3_ltr': "3' Long Terminal Repeat (regulatory)",
        'gag': 'Gag capsid protein (structural)',
        'pol': 'Pol polyprotein (RT, integrase)',
        'protease': 'Protease domain',
        'rt': 'Reverse transcriptase domain',
        'rnaseH': 'RNase H domain',
        'integrase': 'Integrase domain',
        '5_utr': "5' UTR (non-coding)",
        '3_utr': "3' UTR (non-coding)",
        'orf1': 'ORF1 (RNA binding/chaperone)',
        'orf2_rt': 'ORF2 reverse transcriptase',
        'orf2_en': 'ORF2 endonuclease',
        '5_tir': "5' Terminal Inverted Repeat",
        '3_tir': "3' Terminal Inverted Repeat",
        'transposase': 'Transposase (DNA cutting/joining)',
        '5_end': "5' end motif",
        'rep_helicase': 'Rep/helicase domain',
        'internal': 'Internal variable region',
        '3_hairpin': "3' hairpin structure",
        'unknown': 'Unknown region',
    }
    return descriptions.get(domain, domain)


def is_coding_domain(domain: str) -> bool:
    """Check if domain is a coding/ORF region."""
    coding_domains = {
        'gag', 'pol', 'protease', 'rt', 'rnaseH', 'integrase',
        'orf1', 'orf2_rt', 'orf2_en', 'transposase', 'rep_helicase'
    }
    return domain in coding_domains


def is_regulatory_domain(domain: str) -> bool:
    """Check if domain is a regulatory/non-coding region."""
    regulatory_domains = {
        '5_ltr', '3_ltr', '5_utr', '3_utr', '5_tir', '3_tir',
        '5_end', '3_hairpin', 'internal'
    }
    return domain in regulatory_domains
