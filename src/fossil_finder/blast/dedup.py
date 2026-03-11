"""Coordinate-based BLAST hit deduplication.

Removes duplicate hits that arise when multiple transcripts of the same
gene produce identical BLAST alignments against the same TE region.
Replaces the legacy scripts/deduplicate_*_te_hits.py.
"""

import pandas as pd


DEFAULT_KEY_COLUMNS = [
    "gene_id", "chrom", "genomic_start", "genomic_end",
    "te_id", "te_start", "te_end",
]


class HitDeduplicator:
    """Remove duplicate BLAST hits by coordinate key.

    Args:
        key_columns: Columns forming the dedup key. Defaults to the
            legacy 7-tuple (gene_id, chrom, genomic_start, genomic_end,
            te_id, te_start, te_end).
        normalize_te_coords: If True, normalize te_start/te_end so
            start < end before computing the key (handles forward/reverse
            hits to the same TE region).
    """

    def __init__(
        self,
        key_columns: list[str] | None = None,
        normalize_te_coords: bool = True,
    ):
        self.key_columns = key_columns or DEFAULT_KEY_COLUMNS
        self.normalize_te_coords = normalize_te_coords
        self.stats: dict = {
            "total_input": 0,
            "unique_hits": 0,
            "duplicates_removed": 0,
            "duplication_rate": 0.0,
        }

    def deduplicate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove duplicate hits from annotated BLAST results.

        Args:
            df: DataFrame with columns matching key_columns.

        Returns:
            Deduplicated DataFrame (first occurrence kept).
        """
        if df.empty:
            self.stats = {
                "total_input": 0,
                "unique_hits": 0,
                "duplicates_removed": 0,
                "duplication_rate": 0.0,
            }
            return df.copy()

        work = df.copy()

        if self.normalize_te_coords and "te_start" in work.columns and "te_end" in work.columns:
            te_min = work[["te_start", "te_end"]].min(axis=1)
            te_max = work[["te_start", "te_end"]].max(axis=1)
            work["te_start"] = te_min
            work["te_end"] = te_max

        result = work.drop_duplicates(subset=self.key_columns, keep="first")

        total = len(df)
        unique = len(result)
        self.stats = {
            "total_input": total,
            "unique_hits": unique,
            "duplicates_removed": total - unique,
            "duplication_rate": (total - unique) / total if total > 0 else 0.0,
        }

        return result.reset_index(drop=True)

    def per_gene_stats(
        self, df: pd.DataFrame, gene_col: str = "gene_id",
    ) -> dict[str, dict]:
        """Compute per-gene duplication statistics.

        Args:
            df: The original (pre-dedup) DataFrame.
            gene_col: Column containing gene identifiers.

        Returns:
            Dict mapping gene_id -> {raw_hits, unique_hits, duplicates_removed}.
        """
        if df.empty or gene_col not in df.columns:
            return {}

        result = {}
        for gene_id, group in df.groupby(gene_col):
            # Use a temporary deduplicator to avoid corrupting self.stats
            temp = HitDeduplicator(
                key_columns=self.key_columns,
                normalize_te_coords=self.normalize_te_coords,
            )
            deduped = temp.deduplicate(group)
            raw = len(group)
            unique = len(deduped)
            result[gene_id] = {
                "raw_hits": raw,
                "unique_hits": unique,
                "duplicates_removed": raw - unique,
            }

        return result
