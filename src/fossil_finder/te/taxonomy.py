"""TE class inference from names and identifiers.

Supports three classification strategies in priority order:
1. Dfam notation: ID#class/family (most reliable, used by RepeatMasker/Dfam)
2. Configurable family name matching (defaults cover common cross-species families)
3. Keyword fallback (SINE, LTR in name)

Family lists can be overridden per-genome via config for organisms with
non-standard TE naming conventions.
"""

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


def infer_te_class(
    te_name: str,
    family_lists: dict[str, set[str]] | None = None,
) -> str:
    """Infer TE class from a name or identifier.

    Args:
        te_name: TE family name, ID, or Dfam-annotated name.
        family_lists: Optional custom family name lists. If provided,
            replaces DEFAULT_FAMILY_LISTS entirely (use for genomes
            with non-standard TE naming).

    Returns:
        One of: 'LTR', 'LINE', 'DNA', 'SINE', 'Helitron', 'Unknown'.
    """
    name_lower = te_name.lower()

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
