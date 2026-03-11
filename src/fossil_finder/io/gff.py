"""Standards-compliant GFF3 parser for fossil_finder.

Parses any valid GFF3 file (FlyBase, Ensembl, NCBI, custom).
No organism-specific assumptions.

GFF3 spec: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

from pathlib import Path
from urllib.parse import unquote


def _parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 attribute column (col 9) into key-value dict."""
    attrs = {}
    if attr_string == ".":
        return attrs
    for pair in attr_string.split(";"):
        pair = pair.strip()
        if "=" not in pair:
            continue
        key, value = pair.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def parse_gff3(path: str | Path) -> list[dict]:
    """Parse GFF3 file into list of feature dicts.

    Each feature dict has: seqid, source, type, start (int), end (int),
    score, strand, phase, attributes (dict).

    Coordinates are 1-based inclusive (GFF3 standard).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")

    features = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            feature = {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": parts[6],
                "phase": None if parts[7] == "." else int(parts[7]),
                "attributes": _parse_attributes(parts[8]),
            }
            features.append(feature)

    return features


def iter_gff3(path: str | Path):
    """Iterate over GFF3 features one at a time (memory-efficient)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            yield {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": parts[6],
                "phase": None if parts[7] == "." else int(parts[7]),
                "attributes": _parse_attributes(parts[8]),
            }


def get_features_by_type(features: list[dict], feature_type: str) -> list[dict]:
    """Filter features to those matching a specific type (e.g. 'gene', 'exon')."""
    return [f for f in features if f["type"] == feature_type]


def get_children(features: list[dict], parent_id: str) -> list[dict]:
    """Get all features whose Parent attribute matches the given ID."""
    return [
        f for f in features
        if parent_id in f["attributes"].get("Parent", "").split(",")
    ]


def get_gene_to_transcripts(features: list[dict]) -> dict[str, list[str]]:
    """Build mapping from gene IDs to their transcript (mRNA) IDs."""
    genes = {
        f["attributes"]["ID"]: []
        for f in features
        if f["type"] == "gene" and "ID" in f["attributes"]
    }

    for f in features:
        if f["type"] == "mRNA" and "Parent" in f["attributes"]:
            parent = f["attributes"]["Parent"]
            if parent in genes and "ID" in f["attributes"]:
                genes[parent].append(f["attributes"]["ID"])

    return genes
