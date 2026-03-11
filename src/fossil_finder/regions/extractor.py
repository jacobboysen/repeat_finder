"""Config-driven genomic region extractor.

Replaces 6 Dmel-specific extraction scripts with a single generic class
that works with any genome. Reads GFF3 + genome FASTA, extracts regions
by feature type or as windows around anchor features, with optional
deduplication and metadata output.
"""

from pathlib import Path

from fossil_finder.config.schema import GenomeConfig
from fossil_finder.io.fasta import parse_fasta, write_fasta
from fossil_finder.io.gff import get_features_by_type, parse_gff3
from fossil_finder.regions.sequence import (
    SequenceDeduplicator,
    extract_subsequence,
)


class RegionExtractor:
    """Extract genomic regions from GFF3 annotation + genome FASTA.

    Two extraction modes:
    1. Feature-based: extract exact coordinates of GFF3 features
       (three_prime_UTR, exon, CDS, regulatory_region, etc.)
    2. Window-based: extract a window around anchor features
       (e.g., 2kb upstream + 500bp downstream of gene/TSS starts)
    """

    def __init__(self, config: GenomeConfig, base_dir: Path | None = None):
        """Initialize with genome config.

        Args:
            config: Validated GenomeConfig.
            base_dir: Optional base directory for resolving relative paths
                in the config. If None, paths are used as-is.
        """
        self.config = config
        self._base_dir = base_dir

        # Load genome
        genome_path = self._resolve(config.source.genome_fasta)
        self.genome = parse_fasta(genome_path)

        # Filter to configured chromosomes if specified
        if config.genome.chromosomes:
            self.genome = {
                k: v for k, v in self.genome.items()
                if k in config.genome.chromosomes
            }

        # Load annotation
        gff_path = self._resolve(config.source.annotation_gff)
        self.features = parse_gff3(gff_path)

    def _resolve(self, path_str: str) -> Path:
        """Resolve a path, optionally relative to base_dir."""
        p = Path(path_str)
        if not p.is_absolute() and self._base_dir:
            p = self._base_dir / p
        return p

    def extract_features(
        self,
        feature_type: str,
        gene_ids: list[str] | None = None,
        deduplicate: bool = False,
    ) -> list[dict]:
        """Extract regions matching a GFF3 feature type.

        Args:
            feature_type: GFF3 feature type (e.g., 'three_prime_UTR', 'exon').
            gene_ids: Optional list of gene IDs to filter by. Only features
                whose ancestor gene matches will be included.
            deduplicate: If True, collapse identical sequences.

        Returns:
            List of region dicts with keys: region_id, sequence, chrom,
            start, end, strand, parent_id, feature_type.
        """
        matched = get_features_by_type(self.features, feature_type)

        gene_filter = set(gene_ids) if gene_ids else None
        mrna_to_gene = self._build_mrna_to_gene_map() if gene_filter else {}

        dedup = SequenceDeduplicator() if deduplicate else None
        regions = []

        for feat in matched:
            chrom = feat["seqid"]
            if chrom not in self.genome:
                continue

            parent_id = feat["attributes"].get("Parent", "")
            region_id = feat["attributes"].get(
                "ID", f"{chrom}:{feat['start']}-{feat['end']}"
            )

            # Filter by gene if requested
            if gene_filter:
                gene_id = mrna_to_gene.get(parent_id, parent_id)
                if gene_id not in gene_filter:
                    continue

            seq = extract_subsequence(
                self.genome, chrom, feat["start"], feat["end"], feat["strand"]
            )

            region = {
                "region_id": region_id,
                "sequence": seq,
                "chrom": chrom,
                "start": feat["start"],
                "end": feat["end"],
                "strand": feat["strand"],
                "parent_id": parent_id,
                "feature_type": feature_type,
            }

            if dedup:
                if dedup.add(region_id, seq, metadata=region):
                    regions.append(region)
            else:
                regions.append(region)

        return regions

    def extract_windows(
        self,
        anchor_type: str,
        upstream: int,
        downstream: int,
        gene_ids: list[str] | None = None,
        deduplicate: bool = False,
    ) -> list[dict]:
        """Extract windows around anchor features.

        For plus strand: window is [anchor_start - upstream, anchor_start + downstream].
        For minus strand: window is [anchor_start - downstream, anchor_start + upstream].
        The extracted sequence is always reverse-complemented for minus strand.

        Args:
            anchor_type: GFF3 feature type to use as anchor (e.g., 'gene', 'TSS').
            upstream: Base pairs upstream of anchor start.
            downstream: Base pairs downstream of anchor start.
            gene_ids: Optional gene ID filter.
            deduplicate: Collapse identical sequences.

        Returns:
            List of region dicts.
        """
        anchors = get_features_by_type(self.features, anchor_type)
        gene_filter = set(gene_ids) if gene_ids else None
        dedup = SequenceDeduplicator() if deduplicate else None
        regions = []

        for feat in anchors:
            chrom = feat["seqid"]
            if chrom not in self.genome:
                continue

            anchor_id = feat["attributes"].get("ID", "")
            if gene_filter and anchor_id not in gene_filter:
                continue

            strand = feat["strand"]

            if strand == "-":
                # Biological start (TSS) is at feat["end"] for minus-strand
                anchor_pos = feat["end"]
                # "Upstream" is to the right in genome coords
                win_start = anchor_pos - downstream
                win_end = anchor_pos + upstream
            else:
                anchor_pos = feat["start"]
                win_start = anchor_pos - upstream
                win_end = anchor_pos + downstream

            # Clamp to chromosome bounds (GFF3 1-based, so min is 1)
            chrom_len = len(self.genome[chrom])
            win_start = max(1, win_start)
            win_end = min(chrom_len, win_end)

            seq = extract_subsequence(self.genome, chrom, win_start, win_end, strand)

            region = {
                "region_id": f"{anchor_id}_window",
                "sequence": seq,
                "chrom": chrom,
                "start": win_start,
                "end": win_end,
                "strand": strand,
                "parent_id": anchor_id,
                "feature_type": f"{anchor_type}_window",
                "anchor_start": anchor_pos,
                "upstream": upstream,
                "downstream": downstream,
            }

            if dedup:
                if dedup.add(region["region_id"], seq, metadata=region):
                    regions.append(region)
            else:
                regions.append(region)

        return regions

    def write_fasta(self, regions: list[dict], path: str | Path) -> None:
        """Write extracted regions to a FASTA file."""
        sequences = {}
        for r in regions:
            header = (
                f"{r['region_id']} "
                f"loc={r['chrom']}:{r['start']}-{r['end']}:{r['strand']} "
                f"length={len(r['sequence'])}"
            )
            sequences[header] = r["sequence"]

        write_fasta(sequences, path)

    def write_metadata(self, regions: list[dict], path: str | Path) -> None:
        """Write region metadata to a TSV file."""
        path = Path(path)
        columns = [
            "region_id", "chrom", "start", "end", "strand",
            "parent_id", "feature_type", "length",
        ]

        with open(path, "w") as f:
            f.write("\t".join(columns) + "\n")
            for r in regions:
                row = [
                    r.get("region_id", ""),
                    r.get("chrom", ""),
                    str(r.get("start", "")),
                    str(r.get("end", "")),
                    r.get("strand", ""),
                    r.get("parent_id", ""),
                    r.get("feature_type", ""),
                    str(len(r.get("sequence", ""))),
                ]
                f.write("\t".join(row) + "\n")

    def _build_mrna_to_gene_map(self) -> dict[str, str]:
        """Build mRNA ID -> gene ID mapping from GFF3 features."""
        mrna_to_gene = {}
        for feat in self.features:
            if feat["type"] == "mRNA":
                mrna_id = feat["attributes"].get("ID", "")
                gene_id = feat["attributes"].get("Parent", "")
                if mrna_id and gene_id:
                    mrna_to_gene[mrna_id] = gene_id
        return mrna_to_gene
