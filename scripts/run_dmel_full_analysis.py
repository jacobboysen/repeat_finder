#!/usr/bin/env python3
"""Full dmel TE fossil analysis: 3'UTR + 5'UTR, evalue=10, all gene sets.

Runs the complete fossil_finder pipeline on Drosophila melanogaster dm6:
  1. Extract all UTR regions (deduplicated)
  2. BLAST vs TE consensus database (evalue=10, dust=no)
  3. Build region→gene mapping for all extracted UTRs
  4. Run analysis with enrichment testing for 5 gene sets

USAGE:
    conda run -n fossil-finder python scripts/run_dmel_full_analysis.py

OUTPUTS:
    results/dmel_3utr_e10/  — 3'UTR analysis (gene_stats.tsv, enrichment.json, etc.)
    results/dmel_5utr_e10/  — 5'UTR analysis (same structure)
"""

import json
import shutil
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# BLAST+ fallback for conda envs
_BLAST_FALLBACK = Path.home() / "miniconda3/envs/bioinformatics-program/bin"


def find_blastn():
    """Find blastn binary."""
    blastn = shutil.which("blastn")
    if blastn:
        return blastn
    fallback = _BLAST_FALLBACK / "blastn"
    if fallback.exists():
        return str(fallback)
    return None


def load_gene_list(path):
    """Load gene IDs from a text file (one per line)."""
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.add(line)
    return genes


def build_query_to_gene(regions, gff_features):
    """Build region_id → gene_id mapping from extracted regions + GFF.

    Regions have parent_id (transcript IDs, possibly comma-separated).
    GFF mRNA features map transcript_id → gene_id.
    """
    # Build transcript → gene mapping
    t2g = {}
    for feat in gff_features:
        if feat["type"] == "mRNA":
            tid = feat["attributes"].get("ID", "")
            gid = feat["attributes"].get("Parent", "")
            if tid and gid:
                t2g[tid] = gid

    # Map region_id → gene_id
    q2g = {}
    for r in regions:
        region_id = r["region_id"]
        parents = [p.strip() for p in r["parent_id"].split(",")]
        for parent in parents:
            gene_id = t2g.get(parent, parent)
            if gene_id.startswith("FBgn"):
                q2g[region_id] = gene_id
                break

    return q2g


def run_utr_analysis(feature_type, output_name, config, te_db, gene_sets):
    """Run full extraction → BLAST → analysis for one UTR type."""
    from fossil_finder.blast.runner import BlastRunner
    from fossil_finder.io.gff import parse_gff3
    from fossil_finder.pipeline.runner import PipelineRunner

    output_dir = RESULTS_DIR / output_name
    print(f"\n{'='*70}")
    print(f"  {feature_type} analysis → {output_dir}")
    print(f"{'='*70}")

    runner = PipelineRunner(config=config, output_dir=output_dir)

    # Step 1: Extract regions
    t0 = time.time()
    regions = runner.extract(feature_type=feature_type, deduplicate=True)
    t_extract = time.time() - t0
    print(f"  Extracted {len(regions)} deduplicated {feature_type} regions  [{t_extract:.1f}s]")

    if len(regions) == 0:
        print(f"  ERROR: No regions extracted for {feature_type}!")
        return None

    # Step 2: BLAST
    t0 = time.time()
    blast_runner = BlastRunner(config.blast)
    blast_out = output_dir / "blast_results.tsv"
    blast_runner.run(
        query=output_dir / "regions.fa",
        database=te_db,
        output=blast_out,
    )
    n_hits = sum(1 for _ in open(blast_out))
    t_blast = time.time() - t0
    print(f"  BLAST: {n_hits} hits  [{t_blast:.1f}s]")

    # Step 3: Build query → gene mapping
    t0 = time.time()
    gff_path = config.source.annotation_gff
    if not Path(gff_path).is_absolute():
        gff_path = DATA_DIR / gff_path
    features = parse_gff3(gff_path, feature_types={feature_type, "mRNA"})
    q2g = build_query_to_gene(regions, features)
    t_map = time.time() - t0
    n_genes = len(set(q2g.values()))
    print(f"  Mapped {len(q2g)} regions → {n_genes} genes  [{t_map:.1f}s]")

    # Step 4: Analysis with enrichment
    t0 = time.time()
    result = runner.analyze(
        blast_results=blast_out,
        query_to_gene=q2g,
        gene_sets=gene_sets,
    )
    t_analysis = time.time() - t0

    n_analyzed = len(result.gene_stats) if result.gene_stats is not None else 0
    n_with_hits = int((result.gene_stats["hit_count"] > 0).sum()) if result.gene_stats is not None else 0
    n_families = len(result.family_stats["family_stats"]) if result.family_stats else 0

    print(f"  Analysis: {n_analyzed} genes ({n_with_hits} with hits), "
          f"{n_families} TE families  [{t_analysis:.1f}s]")

    if result.strand_bias and "genome" in result.strand_bias:
        gb = result.strand_bias["genome"]
        print(f"  Strand bias: {gb['sense_pct']:.1f}% sense "
              f"({gb['sense_hits']} sense, {gb['antisense_hits']} antisense)")

    if result.enrichment:
        print(f"  Enrichment tests: {len(result.enrichment)} gene sets")
        for name, enr in result.enrichment.items():
            fisher = enr.get("fisher", {})
            mw = enr.get("mannwhitney", {})
            print(f"    {name:20s}  OR={fisher.get('odds_ratio', 0):.2f}  "
                  f"p_fisher={fisher.get('p_value', 1):.4f}  "
                  f"p_mw={mw.get('p_value', 1):.4f}")

    # List output files
    print(f"\n  Output files:")
    total_bytes = 0
    for f in sorted(output_dir.iterdir()):
        sz = f.stat().st_size
        total_bytes += sz
        print(f"    {f.name:35s} {sz:>12,} bytes")
    print(f"    {'TOTAL':35s} {total_bytes:>12,} bytes")

    total_time = t_extract + t_blast + t_analysis + t_map
    print(f"\n  Total time: {total_time:.1f}s  "
          f"(extract={t_extract:.1f}s  blast={t_blast:.1f}s  "
          f"map={t_map:.1f}s  analysis={t_analysis:.1f}s)")

    return result


def main():
    print("FOSSIL FINDER — Full dmel TE fossil analysis")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Project root: {PROJECT_ROOT}")

    # Pre-flight checks
    blastn = find_blastn()
    if not blastn:
        print("ABORT: blastn not found in PATH or conda envs")
        return 1

    coding_gff = DATA_DIR / "references" / "dmel_annotation_coding.gff3"
    if not coding_gff.exists():
        print(f"ABORT: {coding_gff} not found")
        return 1

    te_db = DATA_DIR / "blastdb" / "dmel_te_combined"
    if not (DATA_DIR / "blastdb" / "dmel_te_combined.nin").exists():
        print("ABORT: dmel_te_combined BLAST DB not found")
        return 1

    genome_fasta = DATA_DIR / "references" / "dmel_genome.fasta"
    if not genome_fasta.exists():
        print("ABORT: dmel_genome.fasta not found")
        return 1

    print(f"BLAST binary: {blastn}")
    print(f"Genome: {genome_fasta}")
    print(f"Annotation: {coding_gff}")
    print(f"TE database: {te_db}")

    # Load config
    from fossil_finder.config.schema import GenomeConfig

    config = GenomeConfig(
        genome={
            "species": "Drosophila melanogaster",
            "assembly": "dm6",
            "chromosomes": ["2L", "2R", "3L", "3R", "X", "4"],
        },
        source={
            "adapter": "flybase",
            "genome_fasta": str(genome_fasta),
            "annotation_gff": str(coding_gff),
            "te_consensus": str(DATA_DIR / "references" / "dmel_te_consensus.fasta"),
        },
        blast={
            "program": blastn,
            "word_size": 7, "gapopen": 2, "gapextend": 1,
            "penalty": -1, "reward": 1, "dust": False, "evalue": 10,
        },
    )

    # Load gene sets for enrichment testing
    gene_list_dir = DATA_DIR / "gene_lists"
    gene_sets = {}
    for name, filename in [
        ("germ_plasm", "germ_plasm_fbgn_ids.txt"),
        ("housekeeping", "housekeeping_fbgn_ids.txt"),
        ("somatic", "somatic_fbgn_ids.txt"),
        ("cleared", "cleared_fbgn_ids.txt"),
        ("adult", "adult_fbgn_ids.txt"),
    ]:
        path = gene_list_dir / filename
        if path.exists():
            gene_sets[name] = load_gene_list(path)
            print(f"  Gene set '{name}': {len(gene_sets[name])} genes")
        else:
            print(f"  WARNING: {path} not found, skipping")

    # Run 3'UTR analysis
    result_3utr = run_utr_analysis(
        "three_prime_UTR", "dmel_3utr_e10", config, te_db, gene_sets,
    )

    # Run 5'UTR analysis
    result_5utr = run_utr_analysis(
        "five_prime_UTR", "dmel_5utr_e10", config, te_db, gene_sets,
    )

    # Final summary
    print(f"\n{'='*70}")
    print("  ANALYSIS COMPLETE")
    print(f"{'='*70}")
    for label, result in [("3'UTR", result_3utr), ("5'UTR", result_5utr)]:
        if result:
            s = result.summary()
            print(f"\n  {label}:")
            for k, v in s.items():
                print(f"    {k}: {v}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
