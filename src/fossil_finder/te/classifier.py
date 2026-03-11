"""TE domain position classifier.

Uses position-based heuristics to classify which region of a TE element
(gag, pol, env, LTR, transposase, etc.) a BLAST hit maps to.

Domain boundaries are approximate and based on typical TE structures:
- LTR retrotransposons: 5'LTR - gag - pol - (env) - 3'LTR
- LINE elements: 5'UTR - ORF1 - ORF2 (RT + endonuclease)
- DNA transposons: TIR - transposase - TIR
- Helitrons: 5'end - Rep/helicase - internal - 3'hairpin

Note: These are heuristic approximations. For precise domain mapping,
use Dfam HMMs or curated per-family annotations.
"""

from fossil_finder.te.taxonomy import infer_te_class

# Domain boundaries as fractions of total element length.
# Based on typical TE structures (not organism-specific).

LTR_DOMAINS = {
    "5_ltr": (0.0, 0.12),
    "gag": (0.12, 0.35),
    "pol": (0.35, 0.88),
    "3_ltr": (0.88, 1.0),
}

POL_SUBDOMAINS = {
    "protease": (0.35, 0.42),
    "rt": (0.42, 0.62),
    "rnaseH": (0.62, 0.72),
    "integrase": (0.72, 0.88),
}

LINE_DOMAINS = {
    "5_utr": (0.0, 0.08),
    "orf1": (0.08, 0.35),
    "orf2_rt": (0.35, 0.70),
    "orf2_en": (0.70, 0.92),
    "3_utr": (0.92, 1.0),
}

DNA_DOMAINS = {
    "5_tir": (0.0, 0.10),
    "transposase": (0.10, 0.90),
    "3_tir": (0.90, 1.0),
}

HELITRON_DOMAINS = {
    "5_end": (0.0, 0.05),
    "rep_helicase": (0.05, 0.60),
    "internal": (0.60, 0.95),
    "3_hairpin": (0.95, 1.0),
}

CLASS_TO_DOMAINS = {
    "LTR": LTR_DOMAINS,
    "LINE": LINE_DOMAINS,
    "DNA": DNA_DOMAINS,
    "Helitron": HELITRON_DOMAINS,
}

_CODING_DOMAINS = {
    "gag", "pol", "protease", "rt", "rnaseH", "integrase",
    "orf1", "orf2_rt", "orf2_en", "transposase", "rep_helicase",
}

_REGULATORY_DOMAINS = {
    "5_ltr", "3_ltr", "5_utr", "3_utr", "5_tir", "3_tir",
    "5_end", "3_hairpin", "internal",
}

_DOMAIN_DESCRIPTIONS = {
    "5_ltr": "5' Long Terminal Repeat (regulatory)",
    "3_ltr": "3' Long Terminal Repeat (regulatory)",
    "gag": "Gag capsid protein (structural)",
    "pol": "Pol polyprotein (RT, integrase)",
    "protease": "Protease domain",
    "rt": "Reverse transcriptase domain",
    "rnaseH": "RNase H domain",
    "integrase": "Integrase domain",
    "5_utr": "5' UTR (non-coding)",
    "3_utr": "3' UTR (non-coding)",
    "orf1": "ORF1 (RNA binding/chaperone)",
    "orf2_rt": "ORF2 reverse transcriptase",
    "orf2_en": "ORF2 endonuclease",
    "5_tir": "5' Terminal Inverted Repeat",
    "3_tir": "3' Terminal Inverted Repeat",
    "transposase": "Transposase (DNA cutting/joining)",
    "5_end": "5' end motif",
    "rep_helicase": "Rep/helicase domain",
    "internal": "Internal variable region",
    "3_hairpin": "3' hairpin structure",
    "unknown": "Unknown region",
}


def get_relative_position(sstart: int, send: int, te_length: int) -> tuple[float, float]:
    """Calculate relative position within TE as fraction (0-1).

    Handles both orientations (sstart < send and sstart > send).
    Returns (0.0, 0.0) if te_length <= 0.
    """
    if te_length <= 0:
        return (0.0, 0.0)
    pos1 = min(sstart, send) / te_length
    pos2 = max(sstart, send) / te_length
    return (
        max(0.0, min(1.0, pos1)),
        max(0.0, min(1.0, pos2)),
    )


def classify_te_domain(
    te_id: str,
    sstart: int,
    send: int,
    te_length: int,
    te_class: str | None = None,
) -> dict:
    """Classify which TE domain a BLAST hit maps to.

    Args:
        te_id: TE identifier (used to infer class if not provided).
        sstart: Subject alignment start position.
        send: Subject alignment end position.
        te_length: Total length of TE element.
        te_class: Optional pre-computed TE class.

    Returns:
        Dict with keys: te_class, domain, domain_detail, position,
        rel_start, rel_end.
    """
    if te_class is None:
        te_class = infer_te_class(te_id)

    rel_start, rel_end = get_relative_position(sstart, send, te_length)
    midpoint = (rel_start + rel_end) / 2

    domains = CLASS_TO_DOMAINS.get(te_class, {})
    primary_domain = "unknown"
    domain_detail = None

    for domain_name, (dom_start, dom_end) in domains.items():
        if dom_start <= midpoint <= dom_end:
            primary_domain = domain_name
            break

    if te_class == "LTR" and primary_domain == "pol":
        for sub_name, (sub_start, sub_end) in POL_SUBDOMAINS.items():
            if sub_start <= midpoint <= sub_end:
                domain_detail = sub_name
                break

    if midpoint < 0.2:
        position = "start"
    elif midpoint > 0.8:
        position = "end"
    else:
        position = "middle"

    return {
        "te_class": te_class,
        "domain": primary_domain,
        "domain_detail": domain_detail,
        "position": position,
        "rel_start": round(rel_start, 3),
        "rel_end": round(rel_end, 3),
    }


def get_domain_description(domain: str) -> str:
    """Get human-readable description of a TE domain."""
    return _DOMAIN_DESCRIPTIONS.get(domain, domain)


def is_coding_domain(domain: str) -> bool:
    """Check if domain is a coding/ORF region."""
    return domain in _CODING_DOMAINS


def is_regulatory_domain(domain: str) -> bool:
    """Check if domain is a regulatory/non-coding region."""
    return domain in _REGULATORY_DOMAINS
