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


def _parse_feature(parts: list[str]) -> dict:
    """Parse a split GFF3 line (9 tab-separated fields) into a feature dict."""
    return {
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


def parse_gff3(
    path: str | Path,
    feature_types: set[str] | None = None,
    infer_missing_utrs: bool = True,
) -> list[dict]:
    """Parse GFF3 file into list of feature dicts.

    Each feature dict has: seqid, source, type, start (int), end (int),
    score, strand, phase, attributes (dict).

    Coordinates are 1-based inclusive (GFF3 standard).

    Args:
        path: Path to GFF3 file.
        feature_types: If provided, only load features whose type is in
            this set. Dramatically reduces memory for large GFF3 files
            (e.g., FlyBase dmel has 33M lines but only ~29k three_prime_UTR).
        infer_missing_utrs: If True and the requested feature_types include
            UTR types but the GFF3 has none, infer them from mRNA/CDS
            coordinates. Common for NCBI Gnomon annotations.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")

    utr_types = {"three_prime_UTR", "five_prime_UTR"}
    wants_utrs = feature_types is not None and bool(feature_types & utr_types)

    # If UTRs requested, also load mRNA + CDS for potential inference
    load_types = feature_types
    if wants_utrs and infer_missing_utrs and load_types is not None:
        load_types = load_types | {"mRNA", "CDS"}

    features = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            if load_types is not None and parts[2] not in load_types:
                continue

            features.append(_parse_feature(parts))

    # Infer UTRs if none were found in the GFF3
    if wants_utrs and infer_missing_utrs:
        existing_utrs = [f for f in features if f["type"] in utr_types]
        if not existing_utrs:
            inferred = infer_utrs(features)
            features.extend(inferred)

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
            yield _parse_feature(parts)


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


def infer_utrs(features: list[dict]) -> list[dict]:
    """Infer 3'UTR and 5'UTR features from mRNA and CDS coordinates.

    For GFF3 files that lack explicit UTR annotations (e.g., NCBI Gnomon),
    UTRs are computed as the mRNA regions beyond CDS boundaries.

    For plus-strand transcripts:
      5'UTR = mRNA start to first CDS start - 1
      3'UTR = last CDS end + 1 to mRNA end

    For minus-strand transcripts:
      5'UTR = last CDS end + 1 to mRNA end  (5' in biological sense)
      3'UTR = mRNA start to first CDS start - 1  (3' in biological sense)

    Returns synthetic feature dicts with type 'three_prime_UTR' or
    'five_prime_UTR', using the same dict structure as parse_gff3.
    """
    # Index mRNAs and their CDS children
    mrnas = {}
    cds_by_mrna: dict[str, list[dict]] = {}

    for f in features:
        if f["type"] == "mRNA" and "ID" in f["attributes"]:
            mrna_id = f["attributes"]["ID"]
            mrnas[mrna_id] = f
            cds_by_mrna[mrna_id] = []
        elif f["type"] == "CDS" and "Parent" in f["attributes"]:
            parent = f["attributes"]["Parent"]
            if parent in cds_by_mrna:
                cds_by_mrna[parent].append(f)

    utrs = []
    for mrna_id, mrna in mrnas.items():
        cds_list = cds_by_mrna.get(mrna_id, [])
        if not cds_list:
            continue  # non-coding transcript

        cds_starts = [c["start"] for c in cds_list]
        cds_ends = [c["end"] for c in cds_list]
        cds_min = min(cds_starts)
        cds_max = max(cds_ends)

        strand = mrna["strand"]
        gene_parent = mrna["attributes"].get("Parent", "")

        # Region before CDS (in genome coordinates)
        if mrna["start"] < cds_min:
            utr_type = "five_prime_UTR" if strand == "+" else "three_prime_UTR"
            utrs.append({
                "seqid": mrna["seqid"],
                "source": mrna["source"],
                "type": utr_type,
                "start": mrna["start"],
                "end": cds_min - 1,
                "score": None,
                "strand": strand,
                "phase": None,
                "attributes": {
                    "ID": f"{mrna_id}:{utr_type}:pre_cds",
                    "Parent": mrna_id,
                    "gene_parent": gene_parent,
                },
            })

        # Region after CDS (in genome coordinates)
        if cds_max < mrna["end"]:
            utr_type = "three_prime_UTR" if strand == "+" else "five_prime_UTR"
            utrs.append({
                "seqid": mrna["seqid"],
                "source": mrna["source"],
                "type": utr_type,
                "start": cds_max + 1,
                "end": mrna["end"],
                "score": None,
                "strand": strand,
                "phase": None,
                "attributes": {
                    "ID": f"{mrna_id}:{utr_type}:post_cds",
                    "Parent": mrna_id,
                    "gene_parent": gene_parent,
                },
            })

    return utrs
