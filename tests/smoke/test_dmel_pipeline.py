#!/usr/bin/env python3
"""SMOKE TEST: End-to-end dmel TE fossil pipeline.

PURPOSE
-------
Validate that fossil_finder can run a complete analysis on real
Drosophila melanogaster data: GFF3 filtering → region extraction →
BLAST → analysis. NOT a unit test — exercises the real data path.

LOCATION
--------
tests/smoke/ — explicitly separated from unit tests to prevent
accumulation of one-off scripts in scripts/.

USAGE
-----
    conda run -n fossil-finder python tests/smoke/test_dmel_pipeline.py

EXPECTS
-------
- BLAST+ in PATH (fossil-finder env has it)
- Real dmel data in data/references/ and data/blastdb/
- Pre-filtered GFF: data/references/dmel_annotation_coding.gff3
- fossil_finder package importable (pip install -e . or sys.path hack)

OUTPUTS
-------
All outputs go to a temp directory that is printed at the start and
cleaned up automatically. Nothing is written to results/.
"""

import json
import shutil
import sys
import tempfile
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

DATA_DIR = PROJECT_ROOT / "data"

# If BLAST+ isn't in PATH, check the bioinformatics-program conda env
_BLAST_FALLBACK = Path.home() / "miniconda3/envs/bioinformatics-program/bin"

# 5 germ plasm genes — the core of the TE fossil hypothesis
TEST_GENES = [
    "FBgn0002962",  # nanos
    "FBgn0003015",  # oskar
    "FBgn0016053",  # pgc
    "FBgn0005695",  # gcl (germ cell-less)
    "FBgn0283442",  # cycB (cyclin B)
]


class Timer:
    """Simple timing context manager."""

    def __init__(self, label):
        self.label = label
        self.elapsed = 0.0

    def __enter__(self):
        print(f"\n{'='*60}")
        print(f"  {self.label}")
        print(f"{'='*60}")
        self._start = time.time()
        return self

    def __exit__(self, *args):
        self.elapsed = time.time() - self._start
        print(f"  [{self.elapsed:.1f}s]")


def filter_gff_for_genes(gff_path, gene_ids, output_path):
    """Stream-filter a GFF3 to features related to specific genes.

    Two-pass approach (handles FlyBase's comma-separated Parent fields):
      Pass 1: Collect transcript IDs (mRNA) parented by target genes
      Pass 2: Write gene/mRNA/subfeature lines that match

    Returns (n_features_written, transcript_ids).
    """
    gene_set = set(gene_ids)
    transcript_ids = set()

    # Pass 1: Find transcripts for our genes
    print("  Pass 1: Collecting transcript IDs...")
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or "\t" not in line:
                continue
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = parts[8]
            # Check if Parent contains any of our gene IDs
            for gene_id in gene_set:
                if f"Parent={gene_id}" in attrs:
                    # Extract the transcript ID
                    for pair in attrs.split(";"):
                        if pair.startswith("ID="):
                            transcript_ids.add(pair.split("=", 1)[1].strip())
                            break
                    break

    print(f"  Found {len(transcript_ids)} transcripts for {len(gene_set)} genes")
    for tid in sorted(transcript_ids)[:5]:
        print(f"    {tid}")

    # Pass 2: Write filtered GFF
    print("  Pass 2: Writing filtered GFF...")
    kept = 0
    target_ids = gene_set | transcript_ids

    with open(gff_path) as fin, open(output_path, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue

            ftype = parts[2]
            attrs = parts[8]
            keep = False

            if ftype == "gene":
                # Keep gene lines matching our gene IDs
                for gid in gene_set:
                    if f"ID={gid}" in attrs:
                        keep = True
                        break
            elif ftype == "mRNA":
                # Keep mRNA lines parented by our genes
                for gid in gene_set:
                    if f"Parent={gid}" in attrs:
                        keep = True
                        break
            else:
                # Keep subfeatures (UTR, exon, CDS) parented by our transcripts
                # Handle comma-separated Parent values
                for tid in target_ids:
                    if tid in attrs:
                        keep = True
                        break

            if keep:
                fout.write(line)
                kept += 1

    print(f"  Wrote {kept} features")
    return kept, transcript_ids


def build_query_to_gene(features, transcript_ids, gene_set):
    """Build region_id → gene_id mapping for FlyBase UTRs.

    FlyBase UTRs often lack ID attributes and have comma-separated parents.
    Region IDs are the extractor's fallback: "chrom:start-end".
    """
    # transcript → gene mapping
    t2g = {}
    for feat in features:
        if feat["type"] == "mRNA":
            tid = feat["attributes"].get("ID", "")
            gid = feat["attributes"].get("Parent", "")
            if tid and gid:
                t2g[tid] = gid

    # UTR region_id → gene_id
    q2g = {}
    for feat in features:
        if feat["type"] == "three_prime_UTR":
            # Region ID: either ID attribute or fallback
            region_id = feat["attributes"].get(
                "ID", f"{feat['seqid']}:{feat['start']}-{feat['end']}"
            )
            # Handle comma-separated parents
            parents = feat["attributes"].get("Parent", "").split(",")
            for parent in parents:
                parent = parent.strip()
                gene_id = t2g.get(parent, parent)
                if gene_id in gene_set:
                    q2g[region_id] = gene_id
                    break

    return q2g


def main():
    issues = []  # Collect findings

    print("FOSSIL FINDER SMOKE TEST — dmel end-to-end")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Test genes: {', '.join(TEST_GENES)}")

    # Pre-flight checks
    if not (DATA_DIR / "references" / "dmel_genome.fasta").exists():
        print("ABORT: dmel_genome.fasta not found")
        return 1
    coding_gff = DATA_DIR / "references" / "dmel_annotation_coding.gff3"
    if not coding_gff.exists():
        print("ABORT: dmel_annotation_coding.gff3 not found")
        print("  Generate with: grep -E '^#|gene|mRNA|UTR|CDS|exon' dmel-all-r6.66.gff > dmel_annotation_coding.gff3")
        return 1
    blastn_bin = shutil.which("blastn")
    if not blastn_bin:
        # Check conda env fallback
        fallback = _BLAST_FALLBACK / "blastn"
        if fallback.exists():
            blastn_bin = str(fallback)
            print(f"  Using BLAST+ from {_BLAST_FALLBACK}")
        else:
            print("ABORT: blastn not in PATH and not found in conda envs")
            return 1

    te_db = DATA_DIR / "blastdb" / "dmel_te_combined"
    if not (DATA_DIR / "blastdb" / "dmel_te_combined.nin").exists():
        print("ABORT: dmel_te_combined BLAST DB not found")
        return 1

    tmpdir = Path(tempfile.mkdtemp(prefix="fossil_smoke_"))
    output_dir = tmpdir / "output"
    output_dir.mkdir()
    print(f"Temp dir: {tmpdir}")

    try:
        # ── Step 1: Extract 3'UTRs ─────────────────────────────
        with Timer("Load genome + GFF + extract 3'UTR regions") as t:
            from fossil_finder.config.schema import GenomeConfig
            from fossil_finder.pipeline.runner import PipelineRunner

            config = GenomeConfig(
                genome={
                    "species": "Drosophila melanogaster",
                    "assembly": "dm6",
                    "chromosomes": ["2L", "2R", "3L", "3R", "X", "4"],
                },
                source={
                    "adapter": "flybase",
                    "genome_fasta": str(DATA_DIR / "references" / "dmel_genome.fasta"),
                    "annotation_gff": str(coding_gff),
                    "te_consensus": str(DATA_DIR / "references" / "dmel_te_consensus.fasta"),
                },
                blast={
                    "program": blastn_bin,
                    "word_size": 7, "gapopen": 2, "gapextend": 1,
                    "penalty": -1, "reward": 1, "dust": False, "evalue": 10,
                },
            )

            runner = PipelineRunner(config=config, output_dir=output_dir)
            regions = runner.extract(
                feature_type="three_prime_UTR",
                gene_ids=TEST_GENES,
            )
            print(f"  Extracted {len(regions)} 3'UTR regions")

            if len(regions) == 0:
                print("ERROR: No regions extracted!")
                issues.append("EXTRACTION: 0 regions — check GFF feature hierarchy")
                return 1

            for r in regions[:5]:
                print(f"    {r['region_id']:30s}  {r['chrom']}:{r['start']}-{r['end']:>10}  "
                      f"{len(r['sequence']):>5}bp  {r['strand']}")
        extraction_time = t.elapsed

        # ── Step 2: Run BLAST ───────────────────────────────────
        with Timer("BLAST query regions vs TE database") as t:
            from fossil_finder.blast.runner import BlastRunner

            blast_runner = BlastRunner(config.blast)
            blast_out = output_dir / "blast_results.tsv"

            blast_runner.run(
                query=output_dir / "regions.fa",
                database=te_db,
                output=blast_out,
            )

            n_hits = sum(1 for _ in open(blast_out))
            print(f"  BLAST hits: {n_hits}")

            if n_hits == 0:
                issues.append("BLAST: 0 hits — check database or parameters")
        blast_time = t.elapsed

        # ── Step 3: Build query→gene mapping ────────────────────
        with Timer("Build query → gene mapping") as t:
            from fossil_finder.io.gff import parse_gff3

            features = parse_gff3(coding_gff, feature_types={"three_prime_UTR", "mRNA"})
            gene_set = set(TEST_GENES)
            # Collect transcript IDs for our genes
            transcript_ids = set()
            for feat in features:
                if feat["type"] == "mRNA":
                    parent = feat["attributes"].get("Parent", "")
                    if parent in gene_set:
                        tid = feat["attributes"].get("ID", "")
                        if tid:
                            transcript_ids.add(tid)
            query_to_gene = build_query_to_gene(features, transcript_ids, gene_set)

            print(f"  Mapped {len(query_to_gene)} queries → genes")
            for qid, gid in list(query_to_gene.items())[:5]:
                print(f"    {qid:30s} → {gid}")

            if len(query_to_gene) == 0:
                issues.append("MAPPING: 0 query→gene mappings — analysis will be empty")

        # ── Step 4: Run analysis ────────────────────────────────
        with Timer("Full analysis pipeline") as t:
            result = runner.analyze(
                blast_results=blast_out,
                query_to_gene=query_to_gene,
            )

            n_genes = len(result.gene_stats) if result.gene_stats is not None else 0
            n_blast = len(result.blast_hits) if result.blast_hits is not None else 0

            print(f"\n  Genes analyzed:    {n_genes}")
            print(f"  BLAST hits loaded: {n_blast}")

            if result.gene_stats is not None and len(result.gene_stats) > 0:
                print(f"\n  Per-gene stats:")
                for gene_id, row in result.gene_stats.iterrows():
                    print(f"    {gene_id:20s}  hits={int(row.get('hit_count', 0)):3d}  "
                          f"bp={int(row.get('hit_bp', 0)):5d}  "
                          f"density={row.get('density', 0):.4f}")

            if result.strand_bias:
                gb = result.strand_bias.get("genome", {})
                print(f"\n  Genome strand bias: sense={gb.get('sense_pct', 0):.1f}%  "
                      f"n_sense={gb.get('n_sense', 0)}  n_anti={gb.get('n_antisense', 0)}")

            if result.family_stats and "family_stats" in result.family_stats:
                fs = result.family_stats["family_stats"]
                n_fam = len(fs) if hasattr(fs, '__len__') else 0
                print(f"\n  TE families detected: {n_fam}")
                if hasattr(fs, 'head'):
                    print(fs.head(5).to_string())
        analysis_time = t.elapsed

        # ── Step 5: Verify outputs ─────────────────────────────
        with Timer("Verify outputs"):
            summary_path = output_dir / "summary.json"
            assert summary_path.exists(), "summary.json not created!"
            summary = json.loads(summary_path.read_text())
            print(f"  Summary JSON:")
            print(f"  {json.dumps(summary, indent=4)}")

            expected_files = ["regions.fa", "regions.tsv", "blast_results.tsv", "summary.json"]
            for fname in expected_files:
                fpath = output_dir / fname
                if fpath.exists():
                    size = fpath.stat().st_size
                    print(f"  ✓ {fname:25s} ({size:,} bytes)")
                else:
                    print(f"  ✗ {fname:25s} MISSING")
                    issues.append(f"OUTPUT: {fname} missing")

        # ── Report ──────────────────────────────────────────────
        print(f"\n{'='*60}")
        print("  TIMING SUMMARY")
        print(f"{'='*60}")
        print(f"  Extraction:    {extraction_time:6.1f}s  ({len(regions)} regions)")
        print(f"  BLAST:         {blast_time:6.1f}s  ({n_hits} hits)")
        print(f"  Analysis:      {analysis_time:6.1f}s  ({n_genes} genes)")
        total = extraction_time + blast_time + analysis_time
        print(f"  TOTAL:         {total:6.1f}s")

        print(f"\n{'='*60}")
        if issues:
            print("  SMOKE TEST COMPLETED WITH ISSUES")
            for issue in issues:
                print(f"  ⚠  {issue}")
        else:
            print("  SMOKE TEST PASSED")
        print(f"{'='*60}")

        return 1 if any("ERROR" in i or "ABORT" in i for i in issues) else 0

    finally:
        # Clean up
        shutil.rmtree(tmpdir, ignore_errors=True)
        print(f"\nCleaned up {tmpdir}")


if __name__ == "__main__":
    sys.exit(main())
