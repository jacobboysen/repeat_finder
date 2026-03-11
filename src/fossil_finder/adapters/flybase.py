"""FlyBase genome adapter.

Handles FlyBase-specific file formats:
- Gene IDs with FBgn prefix
- FlyBase symbol maps (TSV: FBgn -> symbol)
- FlyBase RNA-Seq RPKM reports and matrix format
- FlyBase GAF (GO annotations)
- FlyBase gene group data
- FlyFISH localization data
- FlyBase TE FASTA with name=...; headers
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

from fossil_finder.adapters.base import GenomeAdapter
from fossil_finder.config.schema import GenomeConfig


class FlyBaseAdapter(GenomeAdapter):
    """Adapter for FlyBase data sources (Drosophila)."""

    def __init__(self, config: GenomeConfig):
        super().__init__(config)
        self._gene_id_prefix = config.source.gene_id_prefix or "FBgn"

    # --- Gene IDs ---

    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load FBgn IDs from a gene list file.

        Handles simple one-ID-per-line files and TSV files where the
        FBgn ID may be in any column.
        """
        path = Path(path)
        prefix = self._gene_id_prefix
        genes = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "\t" in line:
                    for part in line.split("\t"):
                        if part.startswith(prefix):
                            genes.append(part)
                            break
                elif line.startswith(prefix):
                    genes.append(line)

        return genes

    # --- Symbol mapping ---

    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """Load FBgn -> symbol mapping from configured gene_symbol_map."""
        if not self.config.source.gene_symbol_map:
            return {}

        path = Path(self.config.source.gene_symbol_map)
        if not path.exists():
            return {}

        prefix = self._gene_id_prefix
        mapping = {}

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    fbgn = parts[0].strip()
                    symbol = parts[1].strip()
                    if fbgn.startswith(prefix):
                        mapping[fbgn] = symbol

        return mapping

    # --- Expression ---

    def load_expression(self) -> dict[str, dict[str, float]] | None:
        """Load FlyBase RNA-Seq expression data (matrix format).

        Supports the gene_rpkm_matrix format: genes as rows, tissues as columns.
        First columns are metadata (gene_primary_id, gene_symbol, etc.),
        followed by tissue/sample RPKM columns.
        """
        func = self.config.source.functional
        if not func or not func.expression:
            return None

        path = Path(func.expression)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        expression = {}
        metadata_cols = {
            "gene_primary_id", "gene_symbol", "gene_fullname", "gene_type",
            "#gene_primary_id",
        }

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            header = None
            fbgn_col = 0
            tissue_names: list[str] = []
            tissue_indices: list[int] = []

            for row in reader:
                if not row or row[0].startswith("##"):
                    continue

                if header is None:
                    header = row
                    for i, col in enumerate(header):
                        col_clean = col.lstrip("#")
                        if col_clean in metadata_cols:
                            if "primary_id" in col_clean or "fbgn" in col_clean.lower():
                                fbgn_col = i
                            continue
                        tissue_indices.append(i)
                        tissue_names.append(self._clean_tissue_name(col))
                    continue

                if len(row) <= fbgn_col:
                    continue

                gene_id = row[fbgn_col]
                if not gene_id.startswith(prefix):
                    continue

                expression[gene_id] = {}
                for t_idx, col_idx in enumerate(tissue_indices):
                    if col_idx < len(row) and t_idx < len(tissue_names):
                        try:
                            expression[gene_id][tissue_names[t_idx]] = float(row[col_idx])
                        except ValueError:
                            expression[gene_id][tissue_names[t_idx]] = 0.0

        return expression

    @staticmethod
    def _clean_tissue_name(name: str) -> str:
        """Clean FlyBase tissue/sample names for use as dictionary keys.

        Handles FlyBase column naming conventions:
        - FlyAtlas2: RNA-Seq_Profile_FlyAtlas2_Adult_Female_Ovary_(FBlc...)
        - modENCODE: mE_mRNA_A_MateF_4d_ovary_(FBlc...)
        - BCM: BCM_1_E2-4hr_(FBlc...)
        """
        # Remove FlyBase library IDs
        name = re.sub(r"\(FBlc\d+\)", "", name).strip("_")

        # Handle FlyAtlas2 format
        if "FlyAtlas2" in name:
            name = re.sub(r"RNA-Seq_Profile_FlyAtlas2_", "", name)
            return re.sub(r"\s+", "_", name.replace("_", " ").lower())

        # Handle modENCODE format
        if name.startswith("mE_mRNA_"):
            name = name.replace("mE_mRNA_", "")
            for old, new in [
                ("_A_", "_adult_"), ("_L3_", "_larva3_"),
                ("_L1_", "_larva1_"), ("_L2_", "_larva2_"),
                ("MateF", "mated_female"), ("MateM", "mated_male"),
                ("VirF", "virgin_female"), ("AdF", "adult_female"),
                ("AdM", "adult_male"),
            ]:
                name = name.replace(old, new)

        # Handle BCM format
        if name.startswith("BCM_"):
            name = name.replace("BCM_1_", "")
            for old, new in [
                ("E2-4hr", "embryo_2_4hr"), ("E14-16hr", "embryo_14_16hr"),
                ("E2-16hr", "embryo_2_16hr"),
            ]:
                name = name.replace(old, new)

        # Remove common prefixes
        for prefix in [r"^BCM_HGSC_", r"^modENCODE_", r"^Knoblich_mRNA_"]:
            name = re.sub(prefix, "", name)

        # Normalize
        name = name.lower()
        name = re.sub(r"[,;:\s]+", "_", name)
        name = re.sub(r"-+", "_", name)
        name = re.sub(r"_+", "_", name)
        return name.strip("_")

    # --- GO annotations ---

    def load_go_annotations(self) -> dict[str, list[dict]] | None:
        """Load GO annotations from FlyBase GAF file.

        GAF 2.2 format. Filters to gene IDs matching the configured
        prefix and skips negative (NOT) annotations.
        """
        func = self.config.source.functional
        if not func or not func.go_annotations:
            return None

        path = Path(func.go_annotations)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        annotations: dict[str, list[dict]] = defaultdict(list)

        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("!"):
                    continue

                parts = line.split("\t")
                if len(parts) < 9:
                    continue

                gene_id = parts[1]
                qualifier = parts[3]

                if "NOT" in qualifier.upper():
                    continue
                if not gene_id.startswith(prefix):
                    continue

                annotations[gene_id].append({
                    "go_id": parts[4],
                    "aspect": parts[8] if len(parts) > 8 else "",
                    "evidence": parts[6] if len(parts) > 6 else "",
                    "symbol": parts[2],
                    "qualifier": qualifier,
                })

        return dict(annotations)

    # --- Gene groups ---

    def load_gene_groups(self) -> dict[str, list[str]] | None:
        """Load FlyBase gene group data.

        Format: FB_group_id, FB_group_symbol, ..., Group_member_FB_gene_id, ...
        """
        func = self.config.source.functional
        if not func or not func.gene_groups:
            return None

        path = Path(func.gene_groups)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        groups: dict[str, list[str]] = defaultdict(list)

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#") or "FB_group" in row[0]:
                    continue
                if len(row) >= 7:
                    group_symbol = row[1]
                    member_id = row[5]
                    if member_id.startswith(prefix):
                        groups[group_symbol].append(member_id)

        return dict(groups)

    # --- Localization ---

    def load_localization(self) -> dict[str, list[str]] | None:
        """Load FlyFISH subcellular localization data.

        Uses the symbol column (not FBcl clone IDs) and optionally
        maps to FBgn via the configured symbol map.
        """
        func = self.config.source.functional
        if not func or not func.localization:
            return None

        path = Path(func.localization)
        if not path.exists():
            return None

        localizations: dict[str, set[str]] = defaultdict(set)

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f)
            header = None
            gene_col = None
            term_col = None

            for row in reader:
                if not row:
                    continue
                # Skip comment lines
                if row[0].startswith("#"):
                    continue

                if header is None:
                    header = [col.lower() for col in row]
                    for i, col in enumerate(header):
                        if col in ("gene", "symbol"):
                            gene_col = i
                        if col in ("term", "pattern", "localization"):
                            term_col = i
                    continue

                if gene_col is None or term_col is None:
                    continue
                if gene_col >= len(row) or term_col >= len(row):
                    continue

                symbol = row[gene_col].strip()
                term = row[term_col].strip()

                if not symbol or not term:
                    continue
                if term.lower() in ("", "na", "n/a", "none"):
                    continue

                localizations[symbol].add(term.lower())

        return {gene: sorted(pats) for gene, pats in localizations.items()}

    # --- TE metadata ---

    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE metadata from FlyBase-format FASTA.

        FlyBase TE headers use: >FBteXXXXXXX name=...; class=...; ...
        """
        path = Path(path)
        from fossil_finder.io.fasta import parse_fasta

        sequences = parse_fasta(path)
        te_info = {}

        with open(path) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                te_id = line.strip().split()[0][1:]
                meta = self.parse_fasta_metadata(line)
                name = meta.get("name", te_id)
                te_info[te_id] = {
                    "name": name,
                    "te_class": meta.get("class", "Unknown"),
                    "length": len(sequences.get(te_id, "")),
                }

        return te_info

    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract key=value metadata from FlyBase FASTA headers.

        Handles headers like:
            >FBte0000001 name=gypsy12; class=LTR; family=Gypsy; length=200
            >FBtr0070000 type=three_prime_UTR; loc=2L:100..200; parent=FBgn0031081;
        """
        meta = {}
        for match in re.finditer(r"(\w+)=([^;]+)", header):
            key = match.group(1)
            value = match.group(2).strip()
            meta[key] = value
        return meta
