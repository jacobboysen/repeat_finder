"""Shared helpers for comparative analysis scripts.

Handles region-to-gene mapping and ortholog table loading for both
dmel (FlyBase) and dpse (NCBI RefSeq) annotations.
"""

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "comparative"

BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq", "strand",
]


def load_blast(path: Path) -> pd.DataFrame:
    """Load BLAST TSV with standard 17-column format."""
    df = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLUMNS)
    # Normalize TE coordinates so sstart < send always
    swap = df["sstart"] > df["send"]
    df.loc[swap, ["sstart", "send"]] = df.loc[swap, ["send", "sstart"]].values
    return df


def build_region_to_gene_dmel(
    regions_tsv: Path,
    gff_path: Path,
) -> dict[str, str]:
    """Map dmel region_id -> FBgn gene ID via FBtr transcript parents."""
    regions = pd.read_csv(regions_tsv, sep="\t")

    # Build FBtr -> FBgn from GFF
    fbtr_to_fbgn: dict[str, str] = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = {}
            for pair in parts[8].strip().split(";"):
                if "=" in pair:
                    k, v = pair.split("=", 1)
                    attrs[k] = v
            mrna_id = attrs.get("ID", "")
            gene_id = attrs.get("Parent", "")
            if mrna_id and gene_id:
                fbtr_to_fbgn[mrna_id] = gene_id

    mapping = {}
    for _, row in regions.iterrows():
        # parent_id may be comma-separated (multiple transcripts)
        for parent in str(row["parent_id"]).split(","):
            parent = parent.strip()
            gene = fbtr_to_fbgn.get(parent, parent)
            mapping[row["region_id"]] = gene
            break  # use first match
    return mapping


def build_region_to_gene_dpse(
    regions_tsv: Path,
    gff_path: Path,
) -> dict[str, str]:
    """Map dpse region_id -> NCBI GeneID via transcript parents."""
    regions = pd.read_csv(regions_tsv, sep="\t")

    # Build rna-XM_* -> gene-LOC* from GFF
    rna_to_gene: dict[str, str] = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = {}
            for pair in parts[8].strip().split(";"):
                if "=" in pair:
                    k, v = pair.split("=", 1)
                    attrs[k] = v
            mrna_id = attrs.get("ID", "")
            gene_id = attrs.get("Parent", "")
            if mrna_id and gene_id:
                rna_to_gene[mrna_id] = gene_id

    mapping = {}
    for _, row in regions.iterrows():
        parent = str(row["parent_id"]).split(",")[0].strip()
        gene = rna_to_gene.get(parent, parent)
        mapping[row["region_id"]] = gene
    return mapping


def load_orthologs() -> pd.DataFrame:
    """Load dmel<->dpse ortholog table."""
    path = DATA_DIR / "references" / "orthologs" / "dmel_dpse_orthologs.tsv"
    return pd.read_csv(path, sep="\t")


def build_ortholog_map(orthologs: pd.DataFrame) -> dict[str, str]:
    """Build dmel FBgn -> dpse gene-LOC* mapping.

    The ortholog table has dmel_fbgn and dpse_entrez columns.
    dpse GFF uses gene-LOC{entrez} as gene IDs.
    """
    return {
        row["dmel_fbgn"]: f"gene-LOC{row['dpse_entrez']}"
        for _, row in orthologs.iterrows()
        if pd.notna(row["dmel_fbgn"]) and row["dmel_fbgn"] != ""
    }
