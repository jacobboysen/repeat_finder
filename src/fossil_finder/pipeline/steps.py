"""Individual pipeline steps for TE fossil analysis.

Each step is a pure function: explicit inputs in, structured outputs out.
Steps do not depend on each other — the runner wires them together.
"""

from pathlib import Path

import pandas as pd

from fossil_finder.analysis.aggregation import aggregate_by_gene, compute_density
from fossil_finder.analysis.enrichment import test_gene_set_enrichment
from fossil_finder.analysis.families import (
    compute_class_distribution,
    compute_family_stats,
)
from fossil_finder.analysis.repeatmasker import (
    classify_hits,
    find_overlaps,
    parse_repeatmasker_out,
)
from fossil_finder.analysis.strand import (
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)
from fossil_finder.blast.dedup import HitDeduplicator
from fossil_finder.blast.filter import apply_filters
from fossil_finder.config.schema import GenomeConfig
from fossil_finder.io.blast import load_blast_results
from fossil_finder.regions.extractor import RegionExtractor


def step_extract_regions(
    config: GenomeConfig,
    feature_type: str,
    fasta_out: str | Path,
    metadata_out: str | Path,
    gene_ids: list[str] | None = None,
    deduplicate: bool = False,
    force: bool = True,
    base_dir: Path | None = None,
) -> list[dict]:
    """Extract genomic regions and write FASTA + metadata.

    Args:
        config: Validated GenomeConfig.
        feature_type: GFF3 feature type to extract.
        fasta_out: Path to write query FASTA.
        metadata_out: Path to write region metadata TSV.
        gene_ids: Optional gene ID filter.
        deduplicate: Collapse identical sequences.
        force: If False, skip extraction when output files exist.
        base_dir: Base directory for resolving config paths.

    Returns:
        List of region dicts (same as RegionExtractor.extract_features).
    """
    fasta_out = Path(fasta_out)
    metadata_out = Path(metadata_out)

    if not force and fasta_out.exists() and metadata_out.exists():
        return []  # Skip — files already exist

    extractor = RegionExtractor(config, base_dir=base_dir)
    regions = extractor.extract_features(
        feature_type, gene_ids=gene_ids, deduplicate=deduplicate,
    )

    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)

    extractor.write_fasta(regions, fasta_out)
    extractor.write_metadata(regions, metadata_out)

    return regions


def step_load_and_filter(
    blast_results: str | Path,
    max_evalue: float | None = None,
    min_pident: float | None = None,
    min_length: int | None = None,
) -> pd.DataFrame:
    """Load BLAST results and apply quality filters.

    Args:
        blast_results: Path to BLAST TSV output.
        max_evalue: Maximum e-value threshold.
        min_pident: Minimum percent identity threshold.
        min_length: Minimum alignment length threshold.

    Returns:
        Filtered BLAST results DataFrame.
    """
    df = load_blast_results(blast_results)
    return apply_filters(
        df,
        max_evalue=max_evalue,
        min_pident=min_pident,
        min_length=min_length,
    )


def step_deduplicate(
    df: pd.DataFrame,
    key_columns: list[str] | None = None,
    normalize_te_coords: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Deduplicate BLAST hits by coordinate key.

    Args:
        df: BLAST results (or annotated hits) DataFrame.
        key_columns: Columns forming the dedup key.
        normalize_te_coords: Normalize TE start/end before keying.

    Returns:
        Tuple of (deduplicated DataFrame, stats dict).
    """
    dedup = HitDeduplicator(
        key_columns=key_columns,
        normalize_te_coords=normalize_te_coords,
    )
    result = dedup.deduplicate(df)
    return result, dedup.stats


def step_aggregate(
    df: pd.DataFrame,
    query_to_gene: dict[str, str],
) -> pd.DataFrame:
    """Aggregate BLAST hits by gene and compute density.

    Args:
        df: BLAST results DataFrame.
        query_to_gene: Mapping from query ID to gene ID.

    Returns:
        Per-gene stats DataFrame with density column.
    """
    gene_stats = aggregate_by_gene(df, query_to_gene)
    return compute_density(gene_stats)


def step_strand_analysis(df: pd.DataFrame) -> dict:
    """Run strand bias analysis at all three levels.

    Args:
        df: BLAST results with gene_id and strand columns.

    Returns:
        Dict with 'gene', 'te_family', and 'genome' keys.
    """
    return {
        "gene": compute_gene_strand_bias(df),
        "te_family": compute_te_strand_bias(df),
        "genome": compute_genome_strand_bias(df),
    }


def step_family_analysis(
    df: pd.DataFrame,
    te_metadata: dict[str, dict] | None = None,
) -> dict:
    """Compute TE family statistics and optional class distribution.

    Args:
        df: BLAST results DataFrame.
        te_metadata: Optional TE ID -> metadata mapping.

    Returns:
        Dict with 'family_stats' DataFrame and optional 'class_distribution'.
    """
    result = {"family_stats": compute_family_stats(df)}
    if te_metadata is not None:
        result["class_distribution"] = compute_class_distribution(df, te_metadata)
    return result


def step_enrichment_analysis(
    gene_set: set[str],
    te_positive_genes: set[str],
    gene_densities: "pd.Series",
) -> dict:
    """Run enrichment analysis for a gene set.

    Args:
        gene_set: Gene IDs in the set of interest.
        te_positive_genes: All gene IDs with at least one TE hit.
        gene_densities: Series mapping gene_id -> density value.

    Returns:
        Dict with 'fisher' and 'mannwhitney' sub-dicts.
    """
    return test_gene_set_enrichment(gene_set, te_positive_genes, gene_densities)


def step_repeatmasker_overlap(
    blast_hits: pd.DataFrame,
    repeatmasker_path: str | Path,
    query_regions: pd.DataFrame,
) -> dict:
    """Classify BLAST hits as known (RepeatMasker) or novel.

    Args:
        blast_hits: BLAST results DataFrame.
        repeatmasker_path: Path to RepeatMasker .out file.
        query_regions: DataFrame with region_id, chrom, start, end.

    Returns:
        Dict with 'known' DataFrame, 'novel' DataFrame,
        and 'rm_stats' summary dict.
    """
    rm_regions = parse_repeatmasker_out(repeatmasker_path)
    overlaps = find_overlaps(rm_regions, query_regions)

    if overlaps.empty:
        return {
            "known": pd.DataFrame(),
            "novel": blast_hits.copy(),
            "rm_stats": {
                "total_rm_regions": len(rm_regions),
                "overlapping_regions": 0,
                "known_hits": 0,
                "novel_hits": len(blast_hits),
            },
        }

    known, novel = classify_hits(blast_hits, overlaps)

    return {
        "known": known,
        "novel": novel,
        "rm_stats": {
            "total_rm_regions": len(rm_regions),
            "overlapping_regions": len(overlaps),
            "known_hits": len(known),
            "novel_hits": len(novel),
        },
    }
