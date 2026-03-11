"""Abstract base class for genome-specific data adapters.

Each supported data source (FlyBase, Ensembl, NCBI, custom) provides a
concrete adapter that knows how to parse that source's file formats and
ID conventions. The rest of the pipeline uses only the adapter interface,
never organism-specific parsing directly.
"""

from abc import ABC, abstractmethod
from pathlib import Path

from fossil_finder.config.schema import GenomeConfig


class GenomeAdapter(ABC):
    """Interface for organism-specific data loading.

    Adapters are constructed with a validated GenomeConfig and provide
    uniform access to gene IDs, symbols, expression, GO annotations,
    gene groups, localization, and TE metadata.
    """

    def __init__(self, config: GenomeConfig | None):
        self.config = config

    @abstractmethod
    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load gene identifiers from a gene list file.

        Args:
            path: Path to gene list (one ID per line, or TSV with IDs).

        Returns:
            List of gene identifiers matching this organism's ID format.
        """

    @abstractmethod
    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """Load gene ID -> symbol mapping from the configured source.

        Returns:
            Dict mapping gene IDs (e.g. FBgn0000003) to symbols (e.g. Act5C).
            Empty dict if no symbol map is configured.
        """

    @abstractmethod
    def load_expression(self) -> dict[str, dict[str, float]] | None:
        """Load expression data (e.g. RPKM/TPM matrix).

        Returns:
            Dict mapping gene_id -> {tissue_name: value}, or None if
            expression data is not configured.
        """

    @abstractmethod
    def load_go_annotations(self) -> dict[str, list[dict]] | None:
        """Load Gene Ontology annotations.

        Returns:
            Dict mapping gene_id -> list of GO annotation dicts,
            or None if not configured. Each annotation dict has keys:
            go_id, aspect, evidence, symbol, qualifier.
        """

    @abstractmethod
    def load_gene_groups(self) -> dict[str, list[str]] | None:
        """Load gene groups / pathway memberships.

        Returns:
            Dict mapping group_name -> list of gene IDs,
            or None if not configured.
        """

    @abstractmethod
    def load_localization(self) -> dict[str, list[str]] | None:
        """Load subcellular localization data.

        Returns:
            Dict mapping gene_id -> list of localization patterns,
            or None if not configured.
        """

    @abstractmethod
    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE database metadata from a FASTA file.

        Returns:
            Dict mapping te_id -> {name, te_class, length}.
        """

    @abstractmethod
    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract structured metadata from a FASTA header line.

        Different databases encode metadata differently in headers
        (FlyBase uses key=value; pairs, Ensembl uses space-separated
        key:value, NCBI uses pipe-delimited fields).

        Args:
            header: Full FASTA header line (including '>').

        Returns:
            Dict of extracted key-value metadata.
        """
