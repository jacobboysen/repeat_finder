"""Individual pipeline steps for TE fossil analysis.

Each step is a pure function: explicit inputs in, structured outputs out.
Steps do not depend on each other — the runner wires them together.
"""

from pathlib import Path

import pandas as pd

from fossil_finder.analysis.aggregation import aggregate_by_gene, compute_density
from fossil_finder.analysis.conservation import (
    compute_pident_conservation_correlation,
    hits_to_genomic_bed,
    score_with_bigwig,
    summarize_conservation_by_group,
)
from fossil_finder.analysis.enrichment import test_gene_set_enrichment
from fossil_finder.analysis.families import (
    compute_class_distribution,
    compute_family_stats,
)
from fossil_finder.analysis.multiplicity import compute_hit_multiplicity, compute_te_breadth
from fossil_finder.analysis.positional import (
    compute_end_bias,
    compute_positional_profile,
    compute_te_position,
    compute_utr_position,
)
from fossil_finder.analysis.quality_tiers import (
    assign_quality_tiers,
    compute_edit_stats,
    compute_tier_edit_summary,
    summarize_tiers,
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

    extractor = RegionExtractor(
        config, base_dir=base_dir,
        feature_types={feature_type},
    )
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


def step_quality_tiers(
    df: pd.DataFrame,
    strict_pident: float = 85,
    strict_length: int = 100,
    moderate_pident: float = 75,
    moderate_length: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assign quality tiers and compute edit stats.

    Args:
        df: BLAST results DataFrame.
        strict_pident: Min percent identity for strict tier.
        strict_length: Min alignment length for strict tier.
        moderate_pident: Min percent identity for moderate tier.
        moderate_length: Min alignment length for moderate tier.

    Returns:
        Tuple of (annotated DataFrame with tier+edit columns, tier summary DataFrame).
    """
    df = assign_quality_tiers(
        df,
        strict_pident=strict_pident,
        strict_length=strict_length,
        moderate_pident=moderate_pident,
        moderate_length=moderate_length,
    )
    df = compute_edit_stats(df)
    summary = summarize_tiers(df)
    return df, summary


def step_positional_analysis(df: pd.DataFrame) -> dict:
    """Compute positional profiles and end bias.

    Args:
        df: BLAST results with qstart, qlen, sstart, slen, strand columns.

    Returns:
        Dict with 'utr_profile', 'te_profile', and 'end_bias' keys.
    """
    df = compute_utr_position(df)
    df = compute_te_position(df)
    return {
        "utr_profile": compute_positional_profile(df, "normalized_utr_pos"),
        "te_profile": compute_positional_profile(df, "normalized_te_pos"),
        "end_bias": compute_end_bias(df),
    }


def step_multiplicity_analysis(
    df: pd.DataFrame,
    query_to_gene: dict[str, str] | None = None,
) -> dict:
    """Compute hit multiplicity and TE breadth statistics.

    Args:
        df: BLAST results DataFrame.
        query_to_gene: Optional mapping from query ID to gene ID.

    Returns:
        Dict with 'multiplicity' summary and 'te_breadth' DataFrame.
    """
    return {
        "multiplicity": compute_hit_multiplicity(df, query_to_gene),
        "te_breadth": compute_te_breadth(df, query_to_gene),
    }


def step_conservation_analysis(
    df: pd.DataFrame,
    regions: pd.DataFrame,
    bigwig_path: str | Path,
    tool_path: str | Path,
    sample_relaxed: int = 200_000,
) -> dict:
    """Score BLAST hits with phyloP conservation.

    Args:
        df: BLAST results with tier column.
        regions: Region metadata with region_id, chrom, start, end, strand.
        bigwig_path: Path to phyloP bigWig file.
        tool_path: Path to bigWigAverageOverBed binary.
        sample_relaxed: Max relaxed-tier hits to score (for performance).

    Returns:
        Dict with 'scores' DataFrame, 'by_tier' summary,
        'correlation' stats, and 'scored_df' (annotated hits).
    """
    bigwig_path = Path(bigwig_path)
    tool_path = Path(tool_path)

    # Select hits to score: all strict+moderate, sample relaxed
    if "tier" in df.columns:
        non_relaxed = df[df["tier"].isin(["strict", "moderate"])]
        relaxed = df[df["tier"] == "relaxed"]
        n_sample = min(sample_relaxed, len(relaxed))
        if n_sample > 0:
            sampled = relaxed.sample(n=n_sample, random_state=42)
            to_score = pd.concat([non_relaxed, sampled], ignore_index=True)
        else:
            to_score = non_relaxed.copy()
    else:
        to_score = df.copy()

    # Convert to genomic BED
    bed_df = hits_to_genomic_bed(to_score, regions)
    if bed_df.empty:
        return {"scores": pd.DataFrame(), "by_tier": pd.DataFrame(),
                "correlation": None, "scored_df": pd.DataFrame()}

    # Score with bigWig
    scores = score_with_bigwig(bed_df, bigwig_path, tool_path)

    # Merge scores back
    scored = bed_df.merge(
        scores[["hit_id", "mean0", "mean", "covered", "size"]],
        on="hit_id", how="left",
    )
    scored["phylop_mean"] = scored["mean0"]
    scored["phylop_coverage"] = scored["covered"] / scored["size"].clip(lower=1)

    result = {"scored_df": scored, "scores": scores}

    # Summarize by tier if available
    if "tier" in scored.columns:
        result["by_tier"] = summarize_conservation_by_group(scored, "tier")

    # Pident vs conservation correlation
    result["correlation"] = compute_pident_conservation_correlation(scored)

    return result
