"""Pipeline runner for TE fossil analysis.

Orchestrates extraction, BLAST, filtering, dedup, and analysis.
Returns a typed PipelineResult with all outputs.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from fossil_finder.config.schema import GenomeConfig
from fossil_finder.pipeline.steps import (
    step_aggregate,
    step_deduplicate,
    step_enrichment_analysis,
    step_extract_regions,
    step_family_analysis,
    step_load_and_filter,
    step_repeatmasker_overlap,
    step_strand_analysis,
)


@dataclass
class PipelineResult:
    """Typed container for all pipeline outputs."""

    blast_hits: pd.DataFrame | None = None
    gene_stats: pd.DataFrame | None = None
    strand_bias: dict | None = None
    family_stats: dict | None = None
    enrichment: dict = field(default_factory=dict)
    rm_overlap: dict | None = None

    def summary(self) -> dict:
        """Generate a summary dict of key metrics."""
        s: dict = {}

        if self.blast_hits is not None:
            s["total_blast_hits"] = len(self.blast_hits)

        if self.gene_stats is not None:
            s["n_genes_analyzed"] = len(self.gene_stats)
            if "hit_count" in self.gene_stats.columns:
                s["n_genes_with_hits"] = int(
                    (self.gene_stats["hit_count"] > 0).sum()
                )
            if "density" in self.gene_stats.columns and len(self.gene_stats) > 0:
                s["mean_density"] = float(self.gene_stats["density"].mean())

        if self.strand_bias and "genome" in self.strand_bias:
            s["genome_sense_pct"] = self.strand_bias["genome"].get("sense_pct", 0.0)

        if self.family_stats and "family_stats" in self.family_stats:
            s["n_te_families"] = len(self.family_stats["family_stats"])

        if self.rm_overlap and "rm_stats" in self.rm_overlap:
            s["known_hits"] = self.rm_overlap["rm_stats"].get("known_hits", 0)
            s["novel_hits"] = self.rm_overlap["rm_stats"].get("novel_hits", 0)

        s["n_gene_sets_tested"] = len(self.enrichment)

        return s


class PipelineRunner:
    """Orchestrate TE fossil analysis from config to results.

    Usage:
        runner = PipelineRunner(config, output_dir="results/")
        regions = runner.extract(feature_type="three_prime_UTR")
        result = runner.analyze(blast_results="blast_out.tsv",
                                query_to_gene=q2g_map)
    """

    def __init__(
        self,
        config: GenomeConfig,
        output_dir: str | Path = "output",
        base_dir: Path | None = None,
    ):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._base_dir = base_dir

    def extract(
        self,
        feature_type: str = "three_prime_UTR",
        gene_ids: list[str] | None = None,
        deduplicate: bool = False,
        force: bool = True,
    ) -> list[dict]:
        """Extract genomic regions for BLAST queries.

        Args:
            feature_type: GFF3 feature type to extract.
            gene_ids: Optional gene ID filter.
            deduplicate: Collapse identical sequences.
            force: Re-extract even if output exists.

        Returns:
            List of region dicts.
        """
        return step_extract_regions(
            config=self.config,
            feature_type=feature_type,
            fasta_out=self.output_dir / "regions.fa",
            metadata_out=self.output_dir / "regions.tsv",
            gene_ids=gene_ids,
            deduplicate=deduplicate,
            force=force,
            base_dir=self._base_dir,
        )

    def analyze(
        self,
        blast_results: str | Path,
        query_to_gene: dict[str, str],
        max_evalue: float | None = None,
        min_pident: float | None = None,
        min_length: int | None = None,
        gene_sets: dict[str, set[str]] | None = None,
        te_metadata: dict[str, dict] | None = None,
        repeatmasker_path: str | Path | None = None,
        query_regions: pd.DataFrame | None = None,
    ) -> PipelineResult:
        """Run full analysis on BLAST results.

        Args:
            blast_results: Path to BLAST TSV output.
            query_to_gene: Mapping from query ID to gene ID.
            max_evalue: E-value filter threshold.
            min_pident: Percent identity filter threshold.
            min_length: Alignment length filter threshold.
            gene_sets: Named gene sets for enrichment testing.
            te_metadata: TE ID -> metadata mapping for class distribution.
            repeatmasker_path: Path to RepeatMasker .out file.
            query_regions: DataFrame for RM overlap (region_id, chrom, start, end).

        Returns:
            PipelineResult with all analysis outputs.
        """
        result = PipelineResult()

        # Step 1: Load and filter
        df = step_load_and_filter(
            blast_results,
            max_evalue=max_evalue,
            min_pident=min_pident,
            min_length=min_length,
        )
        result.blast_hits = df

        # Step 2: Aggregate by gene
        gene_stats = step_aggregate(df, query_to_gene)
        result.gene_stats = gene_stats

        # Step 3: Strand analysis (needs gene_id column)
        if not df.empty:
            work = df.copy()
            work["gene_id"] = work["qseqid"].map(query_to_gene)
            work = work.dropna(subset=["gene_id"])
            result.strand_bias = step_strand_analysis(work)
        else:
            result.strand_bias = step_strand_analysis(df)

        # Step 4: TE family analysis
        result.family_stats = step_family_analysis(df, te_metadata=te_metadata)

        # Step 5: Enrichment testing (per gene set)
        if gene_sets and len(gene_stats) > 0:
            te_positive = set(
                gene_stats[gene_stats["hit_count"] > 0].index
            )
            for name, genes in gene_sets.items():
                result.enrichment[name] = step_enrichment_analysis(
                    gene_set=genes,
                    te_positive_genes=te_positive,
                    gene_densities=gene_stats["density"],
                )

        # Step 6: RepeatMasker overlap (optional)
        if repeatmasker_path and query_regions is not None:
            result.rm_overlap = step_repeatmasker_overlap(
                blast_hits=df,
                repeatmasker_path=repeatmasker_path,
                query_regions=query_regions,
            )

        # Save summary
        summary = result.summary()
        summary_path = self.output_dir / "summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        return result
