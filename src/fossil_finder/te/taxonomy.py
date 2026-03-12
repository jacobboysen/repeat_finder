"""TE class inference from names and identifiers.

Supports three classification strategies in priority order:
1. Dfam notation: ID#class/family (most reliable, used by RepeatMasker/Dfam)
2. Configurable family name matching (defaults cover common cross-species families)
3. Keyword fallback (SINE, LTR in name)

Additionally supports consensus FASTA-based classification (RepBase format) and
FlyBase instance FASTA parsing for higher-accuracy classification when reference
files are available.

Family lists can be overridden per-genome via config for organisms with
non-standard TE naming conventions.
"""

from __future__ import annotations

import re
from pathlib import Path

# Manual classifications for Drosophila-specific families not present in
# consensus FASTA files (verified against FlyBase/RepBase annotations).
FLYBASE_OVERRIDES: dict[str, str] = {
    "h": "DNA",                # HB-related hAT superfamily
    "het-tag": "LINE",         # telomeric, HeT-A-related
    "antonia": "LTR",          # LTR/Gypsy
    "ninja-dsim-like": "LTR",  # LTR/Gypsy (ninja family)
    "y": "LINE",               # telomeric Y-element
}

# Default family lists -- covers common families found across many species.
# These are biased toward Drosophila (where the pipeline was developed) but
# include widely-distributed families found in vertebrates and plants too.
DEFAULT_FAMILY_LISTS: dict[str, set[str]] = {
    "LTR": {
        "gypsy", "copia", "bel", "pao", "mdg", "roo", "412", "297",
        "blood", "accord", "tirant", "springer", "opus", "diver",
        "quasimodo", "idefix", "invader", "gtwin", "tabor", "stalker",
        "max", "transpac", "micropia", "blastopia", "17.6", "1731",
    },
    "LINE": {
        "jockey", "doc", "i-element", "f-element", "g-element",
        "x-element", "het-a", "tart", "tahre", "r1", "r2", "cr1",
        "rt1", "baggins", "juan", "ivk", "bs", "fw",
    },
    "DNA": {
        "p-element", "hobo", "pogo", "bari", "s-element",
        "transib", "tc1", "mariner", "piggybac", "mite",
        "protop", "dnarep", "looper", "1360",
    },
    "Helitron": {"helitron", "dine"},
}


def parse_consensus_fasta(path: str | Path) -> dict[str, str]:
    """Parse consensus FASTA headers in RepBase format.

    Headers follow the pattern: ``>family_name#class/subclass``
    e.g. ``>gypsy#LTR/Gypsy``, ``>hobo#DNA/hAT``

    Args:
        path: Path to consensus FASTA file.

    Returns:
        Mapping of ``{lowercase_name: top_level_class}``
        e.g. ``{"gypsy": "LTR", "hobo": "DNA"}``.
        "RC" class is normalized to "Helitron".

    Raises:
        FileNotFoundError: If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Consensus FASTA not found: {path}")

    class_map: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if "#" not in header:
                continue
            name, classification = header.split("#", 1)
            top_class = classification.split("/")[0]
            if top_class == "RC":
                top_class = "Helitron"
            class_map[name.lower()] = top_class
    return class_map


def parse_instance_fasta(path: str | Path) -> dict[str, str]:
    """Parse FlyBase-style TE instance FASTA headers.

    Headers follow the pattern:
    ``>FBti0019256 ... name=invader2{}555; ...``

    Args:
        path: Path to TE instance FASTA file.

    Returns:
        Mapping of ``{te_id: family_name}``
        e.g. ``{"FBti0019256": "invader2{}555"}``.

    Raises:
        FileNotFoundError: If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"TE instance FASTA not found: {path}")

    te_map: dict[str, str] = {}
    name_re = re.compile(r"name=([^;]+)")
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            te_id = line.split()[0][1:]
            match = name_re.search(line)
            if match:
                te_map[te_id] = match.group(1).strip()
    return te_map


def strip_instance_suffix(name: str) -> str:
    """Strip FlyBase instance suffix from a TE family name.

    FlyBase encodes instances as ``family{}instance``, e.g.
    ``"invader2{}555"`` or ``"1360{}Eph[1520]"``.  This function
    removes everything from ``{}`` onward.

    Examples:
        >>> strip_instance_suffix("invader2{}555")
        'invader2'
        >>> strip_instance_suffix("gypsy")
        'gypsy'
    """
    idx = name.find("{}")
    return name[:idx] if idx >= 0 else name


def infer_te_class(
    te_name: str,
    family_lists: dict[str, set[str]] | None = None,
    consensus_map: dict[str, str] | None = None,
) -> str:
    """Infer TE class from a name or identifier.

    Args:
        te_name: TE family name, ID, or Dfam-annotated name.
        family_lists: Optional custom family name lists. If provided,
            replaces DEFAULT_FAMILY_LISTS entirely (use for genomes
            with non-standard TE naming).
        consensus_map: Optional mapping from lowercase family name to
            TE class, as returned by :func:`parse_consensus_fasta`.
            When provided, consensus-based lookup is attempted first
            (after stripping instance suffixes).

    Returns:
        One of: 'LTR', 'LINE', 'DNA', 'SINE', 'Helitron', 'Unknown'.
    """
    name_lower = te_name.lower()

    # Strategy 0 (new): Consensus FASTA-based classification
    if consensus_map is not None:
        base_lower = strip_instance_suffix(name_lower)
        # Exact match in consensus map
        if base_lower in consensus_map:
            return consensus_map[base_lower]
        # Try with -element suffix (FlyBase uses "Doc2" but consensus
        # may have "doc2-element")
        if base_lower + "-element" in consensus_map:
            return consensus_map[base_lower + "-element"]
        # Try FLYBASE_OVERRIDES for Drosophila-specific families
        if base_lower in FLYBASE_OVERRIDES:
            return FLYBASE_OVERRIDES[base_lower]

    # Strategy 1: Dfam notation -- ID#class/family
    if "#" in name_lower:
        class_part = name_lower.split("#")[1].split("/")[0]
        if "ltr" in class_part:
            return "LTR"
        if "line" in class_part:
            return "LINE"
        if "dna" in class_part:
            return "DNA"
        if "sine" in class_part:
            return "SINE"
        if "rc" in class_part or "helitron" in class_part:
            return "Helitron"

    # Strategy 2: Family name matching (longest match wins to avoid
    # false positives from short names like "r2" matching "mariner2")
    families = family_lists if family_lists is not None else DEFAULT_FAMILY_LISTS
    best_match: tuple[str, int] | None = None  # (te_class, match_length)
    for te_class, names in families.items():
        for fam in names:
            if fam in name_lower and (best_match is None or len(fam) > best_match[1]):
                best_match = (te_class, len(fam))
    if best_match is not None:
        return best_match[0]

    # Strategy 3: Keyword fallback (only when using defaults)
    if family_lists is None:
        if "sine" in name_lower or name_lower.startswith("ine-"):
            return "SINE"
        if "ltr" in name_lower:
            return "LTR"

    return "Unknown"
