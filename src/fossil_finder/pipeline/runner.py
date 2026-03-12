"""Pipeline runner for TE fossil analysis.

Orchestrates extraction, BLAST, filtering, dedup, and analysis.
Returns a typed PipelineResult with all outputs.
"""

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from fossil_finder.blast.runner import BlastRunner
from fossil_finder.config.schema import GenomeConfig
from fossil_finder.repeatmasker.runner import RepeatMaskerRunner


def _json_safe(obj):
    """Handle non-JSON-serializable floats (inf, nan)."""
    if isinstance(obj, float):
        if math.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        if math.isnan(obj):
            return "NaN"
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")
from fossil_finder.pipeline.steps import (
    step_aggregate,
    step_conservation_analysis,
    step_deduplicate,
    step_enrichment_analysis,
    step_extract_regions,
    step_family_analysis,
    step_load_and_filter,
    step_motif_analysis,
    step_multiplicity_analysis,
    step_positional_analysis,
    step_quality_tiers,
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
    tier_summary: pd.DataFrame | None = None
    positional: dict | None = None
    multiplicity: dict | None = None
    conservation: dict | None = None
    motifs: dict | None = None

    def summary(self) -> dict:
        """Generate a summary dict of key metrics."""
        s: dict = {}

        if self.blast_hits is not None:
            s["total_blast_hits"] = len(self.blast_hits)
            if "tier" in self.blast_hits.columns:
                tier_counts = self.blast_hits["tier"].value_counts().to_dict()
                s["tier_counts"] = {k: int(v) for k, v in tier_counts.items()}

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

        if self.positional and "end_bias" in self.positional:
            s["end_bias_ratio"] = self.positional["end_bias"].get("end_ratio", 0.0)

        if self.conservation and "correlation" in self.conservation:
            corr = self.conservation["correlation"]
            if corr:
                s["pident_phylop_rho"] = corr.get("rho", 0.0)

        s["n_gene_sets_tested"] = len(self.enrichment)

        if self.motifs:
            top = self.motifs.get("top_motifs", [])
            s["n_top_motifs"] = len(top)

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
        deduplicate: bool = True,
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

    def blast(
        self,
        query: str | Path,
        database: str | Path,
        output: str | Path | None = None,
    ) -> Path:
        """Run BLAST search using the pipeline's BlastSpec configuration.

        Args:
            query: Path to query FASTA file.
            database: Path to BLAST database (without extension).
            output: Path for results TSV. Defaults to
                ``self.output_dir / "blast_results.tsv"``.

        Returns:
            Path to the BLAST results file.
        """
        if output is None:
            output = self.output_dir / "blast_results.tsv"
        runner = BlastRunner(self.config.blast)
        return runner.run(query, database, output)

    def repeatmasker(
        self,
        input_fasta: str | Path,
        output_dir: str | Path | None = None,
    ) -> Path:
        """Run RepeatMasker on input sequences.

        Produces a .out annotation file that can be passed to
        ``analyze(repeatmasker_path=...)`` for overlap classification.

        Args:
            input_fasta: Path to genome or region FASTA file.
            output_dir: Directory for RM output files. Defaults to
                ``self.output_dir / "repeatmasker"``.

        Returns:
            Path to the RepeatMasker .out file.
        """
        if output_dir is None:
            output_dir = self.output_dir / "repeatmasker"
        runner = RepeatMaskerRunner(self.config.repeatmasker)
        return runner.run(input_fasta, output_dir)

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
        phylop_bigwig: str | Path | None = None,
        bigwig_tool: str | Path | None = None,
        shuffled_kmer_counts: list[dict[str, int]] | None = None,
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
            phylop_bigwig: Path to phyloP bigWig file for conservation scoring.
            bigwig_tool: Path to bigWigAverageOverBed binary.

        Returns:
            PipelineResult with all analysis outputs.
        """
        result = PipelineResult()

        # Resolve RepeatMasker .out: explicit arg > config path > run de novo
        if repeatmasker_path is None and self.config.source.repeatmasker_out:
            rm_candidate = Path(self.config.source.repeatmasker_out)
            if self._base_dir and not rm_candidate.is_absolute():
                rm_candidate = self._base_dir / rm_candidate
            if rm_candidate.exists():
                repeatmasker_path = rm_candidate

        if repeatmasker_path is None:
            try:
                genome_fasta = Path(self.config.source.genome_fasta)
                if self._base_dir and not genome_fasta.is_absolute():
                    genome_fasta = self._base_dir / genome_fasta
                if genome_fasta.exists():
                    repeatmasker_path = self.repeatmasker(genome_fasta)
            except (RuntimeError, FileNotFoundError):
                pass  # RepeatMasker not installed — skip overlap step

        # Step 1: Load and filter
        df = step_load_and_filter(
            blast_results,
            max_evalue=max_evalue,
            min_pident=min_pident,
            min_length=min_length,
        )

        # Step 2: Quality tiers and edit stats
        df, tier_summary = step_quality_tiers(df)
        result.blast_hits = df
        result.tier_summary = tier_summary

        # Step 3: Aggregate by gene
        gene_stats = step_aggregate(df, query_to_gene)
        result.gene_stats = gene_stats

        # Step 4: Strand analysis (needs gene_id column)
        if not df.empty:
            work = df.copy()
            work["gene_id"] = work["qseqid"].map(query_to_gene)
            work = work.dropna(subset=["gene_id"])
            result.strand_bias = step_strand_analysis(work)
        else:
            result.strand_bias = step_strand_analysis(df)

        # Step 5: TE family analysis
        result.family_stats = step_family_analysis(df, te_metadata=te_metadata)

        # Step 6: Enrichment testing (per gene set)
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

        # Step 7: RepeatMasker overlap (optional)
        if repeatmasker_path and query_regions is not None:
            result.rm_overlap = step_repeatmasker_overlap(
                blast_hits=df,
                repeatmasker_path=repeatmasker_path,
                query_regions=query_regions,
            )

        # Step 8: Positional analysis
        if not df.empty:
            result.positional = step_positional_analysis(df)

        # Step 9: Multiplicity analysis
        if not df.empty:
            result.multiplicity = step_multiplicity_analysis(df, query_to_gene)

        # Step 10: Conservation analysis (optional)
        if phylop_bigwig and bigwig_tool and query_regions is not None:
            try:
                result.conservation = step_conservation_analysis(
                    df, query_regions, phylop_bigwig, bigwig_tool,
                )
            except (FileNotFoundError, RuntimeError):
                pass  # Conservation is optional — skip if tools missing

        # Step 11: Motif analysis
        if not df.empty:
            result.motifs = step_motif_analysis(
                df,
                shuffled_counts=shuffled_kmer_counts,
                query_to_gene=query_to_gene,
                gene_sets=gene_sets,
            )

        # ── Save all outputs ──────────────────────────────────
        self._save_results(result)

        return result

    def _save_results(self, result: PipelineResult) -> None:
        """Save all analysis outputs to flat files."""
        out = self.output_dir

        # BLAST hits (filtered)
        if result.blast_hits is not None and not result.blast_hits.empty:
            result.blast_hits.to_csv(
                out / "blast_hits_filtered.tsv", sep="\t", index=False,
            )

        # Per-gene stats
        if result.gene_stats is not None and not result.gene_stats.empty:
            result.gene_stats.to_csv(
                out / "gene_stats.tsv", sep="\t", index=True,
            )

        # Strand bias — mixed types: DataFrames + dict
        if result.strand_bias:
            strand_out = {}
            for level in ("gene", "te_family"):
                data = result.strand_bias.get(level)
                if data is not None and isinstance(data, pd.DataFrame) and not data.empty:
                    path = out / f"strand_bias_{level}.tsv"
                    data.to_csv(path, sep="\t", index=True)
                    strand_out[level] = str(path.name)
            genome = result.strand_bias.get("genome")
            if genome:
                strand_out["genome"] = genome
            with open(out / "strand_bias.json", "w") as f:
                json.dump(strand_out, f, indent=2)

        # TE family stats
        if result.family_stats:
            fs = result.family_stats.get("family_stats")
            if fs is not None and isinstance(fs, pd.DataFrame) and not fs.empty:
                fs.to_csv(out / "family_stats.tsv", sep="\t", index=True)
            cd = result.family_stats.get("class_distribution")
            if cd:
                with open(out / "class_distribution.json", "w") as f:
                    json.dump(cd, f, indent=2)

        # Enrichment results
        if result.enrichment:
            with open(out / "enrichment.json", "w") as f:
                json.dump(
                    result.enrichment, f, indent=2,
                    default=_json_safe,
                )

        # RepeatMasker overlap
        if result.rm_overlap:
            rm_stats = result.rm_overlap.get("rm_stats")
            if rm_stats:
                with open(out / "rm_overlap.json", "w") as f:
                    json.dump(rm_stats, f, indent=2)
            for key in ("known", "novel"):
                data = result.rm_overlap.get(key)
                if data is not None and isinstance(data, pd.DataFrame) and not data.empty:
                    data.to_csv(
                        out / f"rm_{key}_hits.tsv", sep="\t", index=False,
                    )

        # Tier summary
        if result.tier_summary is not None and not result.tier_summary.empty:
            result.tier_summary.to_csv(
                out / "tier_summary.tsv", sep="\t", index=False,
            )

        # Positional analysis
        if result.positional:
            utr = result.positional.get("utr_profile")
            if utr is not None and isinstance(utr, pd.DataFrame) and not utr.empty:
                utr.to_csv(out / "utr_profile.tsv", sep="\t", index=False)
            te = result.positional.get("te_profile")
            if te is not None and isinstance(te, pd.DataFrame) and not te.empty:
                te.to_csv(out / "te_profile.tsv", sep="\t", index=False)
            end_bias = result.positional.get("end_bias")
            if end_bias:
                with open(out / "end_bias.json", "w") as f:
                    json.dump(end_bias, f, indent=2, default=_json_safe)

        # Multiplicity analysis
        if result.multiplicity:
            te_breadth = result.multiplicity.get("te_breadth")
            if (
                te_breadth is not None
                and isinstance(te_breadth, pd.DataFrame)
                and not te_breadth.empty
            ):
                te_breadth.to_csv(
                    out / "te_breadth.tsv", sep="\t", index=False,
                )
            mult = result.multiplicity.get("multiplicity")
            if mult:
                with open(out / "multiplicity.json", "w") as f:
                    json.dump(mult, f, indent=2, default=_json_safe)

        # Conservation analysis
        if result.conservation:
            scored_df = result.conservation.get("scored_df")
            if (
                scored_df is not None
                and isinstance(scored_df, pd.DataFrame)
                and not scored_df.empty
            ):
                scored_df.to_csv(
                    out / "phylop_by_hit.tsv", sep="\t", index=False,
                )
            by_tier = result.conservation.get("by_tier")
            if (
                by_tier is not None
                and isinstance(by_tier, pd.DataFrame)
                and not by_tier.empty
            ):
                by_tier.to_csv(
                    out / "phylop_by_quality.tsv", sep="\t", index=False,
                )
            corr = result.conservation.get("correlation")
            if corr:
                with open(out / "quality_paradox_stats.json", "w") as f:
                    json.dump(corr, f, indent=2, default=_json_safe)

        # Motif analysis
        if result.motifs:
            enrichment = result.motifs.get("enrichment")
            if (
                enrichment is not None
                and isinstance(enrichment, pd.DataFrame)
                and not enrichment.empty
            ):
                enrichment.to_csv(
                    out / "motif_enrichment.tsv", sep="\t", index=False,
                )
            pos_profile = result.motifs.get("positional_profile")
            if (
                pos_profile is not None
                and isinstance(pos_profile, pd.DataFrame)
                and not pos_profile.empty
            ):
                pos_profile.to_csv(
                    out / "motif_positional_profile.tsv", sep="\t", index=False,
                )
            gs_assoc = result.motifs.get("gene_set_association")
            if (
                gs_assoc is not None
                and isinstance(gs_assoc, pd.DataFrame)
                and not gs_assoc.empty
            ):
                gs_assoc.to_csv(
                    out / "motif_gene_set_association.tsv", sep="\t", index=False,
                )
            top_motifs = result.motifs.get("top_motifs", [])
            kmer_counts = result.motifs.get("kmer_counts", {})
            with open(out / "motif_summary.json", "w") as f:
                json.dump(
                    {"top_motifs": top_motifs, "n_unique_kmers": len(kmer_counts)},
                    f, indent=2,
                )

        # Summary (always last — references other files)
        summary = result.summary()
        with open(out / "summary.json", "w") as f:
            json.dump(summary, f, indent=2)
