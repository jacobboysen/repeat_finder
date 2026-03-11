"""Custom genome adapter for standard GFF3 + FASTA inputs.

Handles genomes with no organism-specific database. Extracts what it
can from standard file formats and FASTA header conventions.
"""

import re
from pathlib import Path

from fossil_finder.adapters.base import GenomeAdapter
from fossil_finder.config.schema import GenomeConfig


class CustomAdapter(GenomeAdapter):
    """Adapter for genomes with standard GFF3 + FASTA (no database)."""

    def __init__(self, config: GenomeConfig):
        super().__init__(config)
        self._gene_id_prefix = config.source.gene_id_prefix or ""

    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load gene IDs, filtering by configured prefix (if any)."""
        path = Path(path)
        prefix = self._gene_id_prefix
        genes = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "\t" in line:
                    gene_id = line.split("\t")[0]
                else:
                    gene_id = line

                if not prefix or gene_id.startswith(prefix):
                    genes.append(gene_id)

        return genes

    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """No symbol map available for custom genomes."""
        return {}

    def load_expression(self) -> None:
        """No expression data for custom genomes."""
        return None

    def load_go_annotations(self) -> None:
        """No GO annotations for custom genomes."""
        return None

    def load_gene_groups(self) -> None:
        """No gene groups for custom genomes."""
        return None

    def load_localization(self) -> None:
        """No localization data for custom genomes."""
        return None

    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE metadata from FASTA with key=value or Dfam-style headers."""
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
                te_info[te_id] = {
                    "name": meta.get("name", te_id),
                    "te_class": meta.get("class", "Unknown"),
                    "length": len(sequences.get(te_id, "")),
                }

        return te_info

    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract metadata from FASTA headers.

        Supports two conventions:
        1. Key=value pairs: >ID name=foo; class=LTR; family=gypsy
        2. Dfam notation: >Gypsy-2_DM#LTR/Gypsy
        """
        meta = {}

        # Try key=value pairs first
        for match in re.finditer(r"(\w+)=([^;]+)", header):
            meta[match.group(1)] = match.group(2).strip()

        # Try Dfam notation: ID#class/family
        if not meta:
            first_token = header.strip().split()[0].lstrip(">")
            if "#" in first_token:
                _, class_family = first_token.split("#", 1)
                parts = class_family.split("/", 1)
                meta["class"] = parts[0]
                if len(parts) > 1:
                    meta["family"] = parts[1]

        return meta
