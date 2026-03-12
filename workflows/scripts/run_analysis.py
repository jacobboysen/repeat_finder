"""Bridge script: Snakemake -> PipelineRunner.analyze().

Called by Snakemake's script: directive. Reads the injected snakemake
object for input/output paths and config, then runs the analysis pipeline.

Also works standalone via argparse (for testing without Snakemake).
"""

from pathlib import Path

import pandas as pd

from fossil_finder.config.schema import GeneSetSpec, load_genome_config
from fossil_finder.pipeline.runner import PipelineRunner


def build_query_to_gene(
    regions_metadata: str | Path,
    annotation_gff: str | Path,
) -> dict[str, str]:
    """Build query_id -> gene_id mapping from regions.tsv and GFF3.

    regions.tsv parent_id is a transcript/mRNA ID (GFF3 Parent of a UTR).
    We resolve mRNA -> gene via the GFF3 hierarchy.
    """
    from fossil_finder.io.gff import parse_gff3

    df = pd.read_csv(regions_metadata, sep="\t")

    # Build mRNA_id -> gene_id lookup from GFF3
    features = parse_gff3(annotation_gff, feature_types={"mRNA"})
    mrna_to_gene: dict[str, str] = {}
    for feat in features:
        if feat["type"] == "mRNA":
            mrna_id = feat["attributes"].get("ID", "")
            gene_id = feat["attributes"].get("Parent", "")
            if mrna_id and gene_id:
                mrna_to_gene[mrna_id] = gene_id

    # Map region_id -> gene_id (fall back to parent_id if mRNA not found)
    return {
        row.region_id: mrna_to_gene.get(row.parent_id, row.parent_id)
        for _, row in df.iterrows()
    }


def load_gene_sets(
    gene_set_specs: dict[str, GeneSetSpec],
    base_dir: Path | None = None,
) -> dict[str, set[str]]:
    """Load gene set files into dict[str, set[str]].

    Each GeneSetSpec.genes points to a file with one gene ID per line.
    """
    result: dict[str, set[str]] = {}
    for name, spec in gene_set_specs.items():
        path = Path(spec.genes)
        if base_dir and not path.is_absolute():
            path = base_dir / path
        if not path.exists():
            continue
        genes = set()
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    genes.add(line.split()[0])
        result[name] = genes
    return result


def main():
    # Determine execution context: Snakemake or standalone
    try:
        sm = snakemake  # noqa: F821 — injected by Snakemake script: directive
        config_path = sm.params.genome_config
        blast_results = sm.input.blast_results
        regions_metadata = sm.input.regions_metadata
        rm_out = sm.input.rm_out
        output_dir = str(Path(sm.output.summary).parent)
        base_dir = sm.config.get("base_dir", None)
        feature_type = sm.config.get("feature_type", "three_prime_UTR")
        max_evalue = sm.config.get("max_evalue", None)
        min_pident = sm.config.get("min_pident", None)
        min_length = sm.config.get("min_length", None)
    except NameError:
        import argparse
        parser = argparse.ArgumentParser(description="Run fossil_finder analysis")
        parser.add_argument("--config", required=True, help="GenomeConfig YAML path")
        parser.add_argument("--blast-results", required=True)
        parser.add_argument("--regions-metadata", required=True)
        parser.add_argument("--rm-out", default=None, help="RepeatMasker .out path")
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--base-dir", default=None)
        parser.add_argument("--feature-type", default="three_prime_UTR")
        parser.add_argument("--max-evalue", type=float, default=None)
        parser.add_argument("--min-pident", type=float, default=None)
        parser.add_argument("--min-length", type=int, default=None)
        args = parser.parse_args()
        config_path = args.config
        blast_results = args.blast_results
        regions_metadata = args.regions_metadata
        rm_out = args.rm_out
        output_dir = args.output_dir
        base_dir = args.base_dir
        feature_type = args.feature_type
        max_evalue = args.max_evalue
        min_pident = args.min_pident
        min_length = args.min_length

    base_dir_path = Path(base_dir) if base_dir else None

    # Load config
    config = load_genome_config(config_path)

    # Resolve annotation GFF path for mRNA -> gene mapping
    annotation_gff = config.source.annotation_gff
    if base_dir_path and not Path(annotation_gff).is_absolute():
        annotation_gff = str(base_dir_path / annotation_gff)

    # Build query -> gene mapping from regions metadata + GFF3
    query_to_gene = build_query_to_gene(regions_metadata, annotation_gff)

    # Load query regions for RM overlap
    query_regions = pd.read_csv(regions_metadata, sep="\t")

    # Load gene sets from config
    gene_sets = load_gene_sets(config.gene_sets, base_dir=base_dir_path)

    # Run analysis
    runner = PipelineRunner(
        config=config,
        output_dir=output_dir,
        base_dir=base_dir_path,
    )

    result = runner.analyze(
        blast_results=blast_results,
        query_to_gene=query_to_gene,
        max_evalue=max_evalue,
        min_pident=min_pident,
        min_length=min_length,
        gene_sets=gene_sets if gene_sets else None,
        repeatmasker_path=rm_out,
        query_regions=query_regions,
    )

    n_genes = len(result.gene_stats) if result.gene_stats is not None else 0
    print(f"Analysis complete. {n_genes} genes analyzed.")
    print(f"Summary saved to {Path(output_dir) / 'summary.json'}")


# Snakemake script: directive injects `snakemake` into globals but does
# not set __name__ to "__main__". Check for both execution contexts.
if "snakemake" in dir():
    main()
elif __name__ == "__main__":
    main()
