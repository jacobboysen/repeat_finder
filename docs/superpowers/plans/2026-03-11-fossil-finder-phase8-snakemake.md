> **SUPERSEDED** by `docs/superpowers/plans/2026-03-12-snakemake-workflow.md` which adds
> multi-genome support, conditional RepeatMasker, and the bridge script pattern.

# Fossil Finder Phase 8: Snakemake Workflow

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the placeholder Snakefile with a production Snakemake workflow that orchestrates the full TE fossil pipeline — from region extraction through BLAST to analysis — using the `fossil_finder` Python package as the computational backend.

**Architecture:** A single `Snakefile` with ~7 rules forming a linear DAG: `extract_regions` → `make_blast_db` → `run_blast` → `analyze`. Config is loaded from the same `GenomeConfig` YAML that the Python package uses. Each rule produces file-based outputs, enabling Snakemake's automatic resume-from-failure. A thin Python wrapper script (`workflow/scripts/run_analysis.py`) bridges the Snakemake file-based world with `PipelineRunner.analyze()`.

**Tech Stack:** Snakemake >=8.0, fossil_finder package, BLAST+

---

## File Structure (Phase 8)

```
workflows/
├── Snakefile                    # Main workflow definition (~7 rules)
├── scripts/
│   └── run_analysis.py          # Python bridge: Snakemake → PipelineRunner.analyze()
tests/
├── test_workflow/
│   ├── __init__.py
│   └── test_snakemake.py        # Dry-run + integration tests (6 tests)
```

**Key design decisions:**

1. **Single Snakefile, not modular rules** — the pipeline is linear with ~7 steps. Splitting into multiple rule files adds complexity for no benefit at this scale.

2. **Config reuse** — the Snakemake workflow reads the same `GenomeConfig` YAML as the Python API. No duplicate config format.

3. **`run_analysis.py` bridge script** — Snakemake rules are best kept to shell commands and `script:` directives. Complex Python logic (aggregation, enrichment, strand analysis) stays in the package; the bridge script just calls `PipelineRunner.analyze()` and writes outputs.

4. **File-based checkpointing** — Snakemake tracks file timestamps natively. No custom skip logic needed (unlike `PipelineRunner.extract(force=...)`).

5. **Thread-aware** — BLAST rule uses `threads:` directive mapped to `BlastSpec.num_threads`.

6. **No cluster config** — HPC profiles (SLURM, SGE) are user-provided via `--profile`. The workflow stays portable.

---

## Chunk 1: Snakefile + Bridge Script

### Task 1: Snakemake workflow and analysis bridge

**Files:**
- Modify: `workflows/Snakefile`
- Create: `workflows/scripts/run_analysis.py`
- Create: `tests/test_workflow/__init__.py`
- Create: `tests/test_workflow/test_snakemake.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_workflow/__init__.py` (empty) and `tests/test_workflow/test_snakemake.py`:

```python
"""Tests for the Snakemake workflow.

Tests use --dry-run (no BLAST+ needed) unless BLAST+ is available.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

WORKFLOWS_DIR = Path(__file__).resolve().parents[2] / "workflows"
SNAKEFILE = WORKFLOWS_DIR / "Snakefile"


def snakemake_available() -> bool:
    return shutil.which("snakemake") is not None


def blast_available() -> bool:
    return shutil.which("blastn") is not None


@pytest.mark.skipif(not snakemake_available(), reason="snakemake not installed")
class TestSnakemakeDryRun:
    """Validate workflow DAG without executing rules."""

    def test_snakefile_parses(self, test_data_dir, tmp_path):
        """Snakefile loads without syntax errors."""
        config_path = test_data_dir / "mini_genome_config.yaml"
        result = subprocess.run(
            [
                "snakemake", "-s", str(SNAKEFILE),
                "--dry-run", "--quiet",
                "--config",
                f"genome_config={config_path}",
                f"output_dir={tmp_path / 'output'}",
                f"base_dir={test_data_dir}",
            ],
            capture_output=True, text=True, cwd=str(WORKFLOWS_DIR),
        )
        assert result.returncode == 0, f"Snakefile parse failed:\n{result.stderr}"

    def test_dag_has_expected_rules(self, test_data_dir, tmp_path):
        """DAG contains the core rules."""
        config_path = test_data_dir / "mini_genome_config.yaml"
        result = subprocess.run(
            [
                "snakemake", "-s", str(SNAKEFILE),
                "--list-rules",
                "--config",
                f"genome_config={config_path}",
                f"output_dir={tmp_path / 'output'}",
                f"base_dir={test_data_dir}",
            ],
            capture_output=True, text=True, cwd=str(WORKFLOWS_DIR),
        )
        rules = result.stdout.strip().split("\n")
        expected = ["extract_regions", "make_blast_db", "run_blast", "analyze"]
        for rule_name in expected:
            assert any(rule_name in r for r in rules), f"Missing rule: {rule_name}"

    def test_config_validation(self, test_data_dir, tmp_path):
        """Workflow fails gracefully with missing config."""
        result = subprocess.run(
            [
                "snakemake", "-s", str(SNAKEFILE),
                "--dry-run", "--quiet",
                "--config",
                f"genome_config=/nonexistent/config.yaml",
                f"output_dir={tmp_path / 'output'}",
            ],
            capture_output=True, text=True, cwd=str(WORKFLOWS_DIR),
        )
        # Should fail (config file missing)
        assert result.returncode != 0


@pytest.mark.skipif(not snakemake_available(), reason="snakemake not installed")
class TestBridgeScript:
    def test_bridge_script_exists(self):
        bridge = WORKFLOWS_DIR / "scripts" / "run_analysis.py"
        assert bridge.exists()

    def test_bridge_script_imports(self):
        """Bridge script has valid Python syntax and imports."""
        bridge = WORKFLOWS_DIR / "scripts" / "run_analysis.py"
        result = subprocess.run(
            ["python", "-c", f"import ast; ast.parse(open('{bridge}').read())"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0


@pytest.mark.skipif(
    not (snakemake_available() and blast_available()),
    reason="snakemake + BLAST+ required",
)
class TestSnakemakeIntegration:
    def test_full_pipeline_mini_genome(self, test_data_dir, tmp_path):
        """End-to-end run on synthetic mini genome."""
        config_path = test_data_dir / "mini_genome_config.yaml"
        result = subprocess.run(
            [
                "snakemake", "-s", str(SNAKEFILE),
                "--cores", "1",
                "--config",
                f"genome_config={config_path}",
                f"output_dir={tmp_path / 'output'}",
                f"base_dir={test_data_dir}",
            ],
            capture_output=True, text=True, cwd=str(WORKFLOWS_DIR),
            timeout=120,
        )
        assert result.returncode == 0, f"Snakemake failed:\n{result.stderr}"

        # Check expected outputs exist
        out = tmp_path / "output"
        assert (out / "regions.fa").exists()
        assert (out / "regions.tsv").exists()
        assert (out / "blast_results.tsv").exists()
        assert (out / "summary.json").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_workflow/test_snakemake.py -v`
Expected: FAIL — rules not defined yet (dry-run will fail or tests skip)

- [ ] **Step 3: Write the bridge script**

Create `workflows/scripts/run_analysis.py`:

```python
"""Bridge script: Snakemake → PipelineRunner.analyze().

Called by Snakemake's `script:` directive. Reads snakemake object
for input/output paths and config, then runs the analysis pipeline.
"""

import sys
from pathlib import Path

import pandas as pd

# Ensure the package is importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from fossil_finder.config.schema import load_genome_config
from fossil_finder.io.gff import parse_gff3
from fossil_finder.pipeline.runner import PipelineRunner


def build_query_to_gene(gff_path: str | Path, feature_type: str) -> dict[str, str]:
    """Build query_id → gene_id mapping from GFF3 annotation.

    For feature-based queries (e.g., three_prime_UTR), the query ID
    in the FASTA header matches the GFF3 feature ID, and we map that
    back to the parent gene via mRNA → gene relationships.
    """
    features = parse_gff3(gff_path)

    # Build mRNA → gene mapping
    mrna_to_gene: dict[str, str] = {}
    for feat in features:
        if feat["type"] == "mRNA":
            mrna_id = feat["attributes"].get("ID", "")
            gene_id = feat["attributes"].get("Parent", "")
            if mrna_id and gene_id:
                mrna_to_gene[mrna_id] = gene_id

    # Build feature_id → gene_id mapping
    query_to_gene: dict[str, str] = {}
    for feat in features:
        if feat["type"] == feature_type:
            feat_id = feat["attributes"].get("ID", "")
            parent = feat["attributes"].get("Parent", "")
            gene_id = mrna_to_gene.get(parent, parent)
            if feat_id and gene_id:
                query_to_gene[feat_id] = gene_id

    return query_to_gene


def main():
    # When called via Snakemake script: directive, `snakemake` object is available
    # When called standalone, read from command-line args
    try:
        # Snakemake script: directive injects `snakemake` into globals
        sm = snakemake  # noqa: F821
        config_path = sm.config["genome_config"]
        blast_results = sm.input.blast_results
        regions_metadata = sm.input.regions_metadata
        output_dir = str(Path(sm.output.summary).parent)
        base_dir = sm.config.get("base_dir", None)
        feature_type = sm.config.get("feature_type", "three_prime_UTR")
        max_evalue = sm.config.get("max_evalue", None)
        min_pident = sm.config.get("min_pident", None)
    except NameError:
        # Standalone mode (for testing)
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--config", required=True)
        parser.add_argument("--blast-results", required=True)
        parser.add_argument("--regions-metadata", required=True)
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--base-dir", default=None)
        parser.add_argument("--feature-type", default="three_prime_UTR")
        parser.add_argument("--max-evalue", type=float, default=None)
        parser.add_argument("--min-pident", type=float, default=None)
        args = parser.parse_args()
        config_path = args.config
        blast_results = args.blast_results
        regions_metadata = args.regions_metadata
        output_dir = args.output_dir
        base_dir = args.base_dir
        feature_type = args.feature_type
        max_evalue = args.max_evalue
        min_pident = args.min_pident

    # Load config
    base_dir_path = Path(base_dir) if base_dir else None
    config = load_genome_config(config_path)

    # Build query -> gene mapping from annotation
    gff_path = config.source.annotation_gff
    if base_dir_path and not Path(gff_path).is_absolute():
        gff_path = str(base_dir_path / gff_path)
    query_to_gene = build_query_to_gene(gff_path, feature_type)

    # Load query regions for RepeatMasker overlap
    query_regions = pd.read_csv(regions_metadata, sep="\t")

    # Run analysis
    runner = PipelineRunner(
        config=config,
        output_dir=output_dir,
        base_dir=base_dir_path,
    )

    # Check for RepeatMasker annotations
    rm_path = config.source.repeatmasker_out
    if rm_path and base_dir_path and not Path(rm_path).is_absolute():
        rm_path = str(base_dir_path / rm_path)

    result = runner.analyze(
        blast_results=blast_results,
        query_to_gene=query_to_gene,
        max_evalue=max_evalue,
        min_pident=min_pident,
        repeatmasker_path=rm_path if rm_path else None,
        query_regions=query_regions if rm_path else None,
    )

    print(f"Analysis complete. {len(result.gene_stats)} genes analyzed.")
    print(f"Summary saved to {Path(output_dir) / 'summary.json'}")


# Called unconditionally: Snakemake script: directive does not set
# __name__ to "__main__". The try/except inside main() handles
# both Snakemake and standalone contexts.
main()
```

- [ ] **Step 4: Write the Snakefile**

Update `workflows/Snakefile`:

```python
"""TE Fossil Finder — Snakemake Workflow.

Usage:
    snakemake --cores 4 --config genome_config=path/to/config.yaml output_dir=results/

Config keys:
    genome_config: Path to GenomeConfig YAML (required)
    output_dir:    Output directory (default: "output")
    base_dir:      Base directory for resolving relative paths in config
    feature_type:  GFF3 feature type to extract (default: "three_prime_UTR")
    max_evalue:    E-value filter (optional)
    min_pident:    Percent identity filter (optional)
"""

import sys
from pathlib import Path

# Ensure fossil_finder is importable
sys.path.insert(0, str(Path(workflow.basedir).parent / "src"))

from fossil_finder.config.schema import load_genome_config


# ---------- Config ----------

GENOME_CONFIG = config.get("genome_config")
if not GENOME_CONFIG:
    raise ValueError("Must provide --config genome_config=path/to/config.yaml")

OUTPUT_DIR = Path(config.get("output_dir", "output"))
BASE_DIR = Path(config["base_dir"]) if config.get("base_dir") else None
FEATURE_TYPE = config.get("feature_type", "three_prime_UTR")

# Load and validate genome config
cfg = load_genome_config(GENOME_CONFIG)
blast_spec = cfg.blast

# Resolve paths relative to base_dir
def resolve_path(p):
    """Resolve a path, optionally relative to BASE_DIR."""
    path = Path(p)
    if not path.is_absolute() and BASE_DIR:
        path = BASE_DIR / path
    return str(path)

GENOME_FASTA = resolve_path(cfg.source.genome_fasta)
ANNOTATION_GFF = resolve_path(cfg.source.annotation_gff)
TE_DB_SOURCE = resolve_path(cfg.source.te_consensus or cfg.source.te_instances or "")

if not TE_DB_SOURCE or not Path(TE_DB_SOURCE).exists():
    raise ValueError(
        f"TE database not found: {TE_DB_SOURCE}. "
        "Set source.te_consensus or source.te_instances in config."
    )


# ---------- Targets ----------

rule all:
    input:
        OUTPUT_DIR / "summary.json",


# ---------- Rules ----------

rule extract_regions:
    """Extract query regions (e.g., 3'UTRs) from genome."""
    input:
        genome=GENOME_FASTA,
        annotation=ANNOTATION_GFF,
    output:
        fasta=OUTPUT_DIR / "regions.fa",
        metadata=OUTPUT_DIR / "regions.tsv",
    params:
        feature_type=FEATURE_TYPE,
    run:
        from fossil_finder.pipeline.steps import step_extract_regions
        step_extract_regions(
            config=cfg,
            feature_type=params.feature_type,
            fasta_out=output.fasta,
            metadata_out=output.metadata,
            base_dir=BASE_DIR,
        )


rule make_blast_db:
    """Build BLAST nucleotide database from TE sequences."""
    input:
        fasta=TE_DB_SOURCE,
    output:
        nhr=OUTPUT_DIR / "blastdb" / "te_db.nhr",
        nin=OUTPUT_DIR / "blastdb" / "te_db.nin",
        nsq=OUTPUT_DIR / "blastdb" / "te_db.nsq",
    params:
        db_prefix=str(OUTPUT_DIR / "blastdb" / "te_db"),
    shell:
        """
        mkdir -p {OUTPUT_DIR}/blastdb
        makeblastdb \
            -in {input.fasta} \
            -dbtype nucl \
            -out {params.db_prefix} \
            -parse_seqids
        """


rule run_blast:
    """Run BLAST search: query regions vs TE database."""
    input:
        query=OUTPUT_DIR / "regions.fa",
        db_nhr=OUTPUT_DIR / "blastdb" / "te_db.nhr",
    output:
        tsv=OUTPUT_DIR / "blast_results.tsv",
    params:
        db_prefix=str(OUTPUT_DIR / "blastdb" / "te_db"),
    threads: blast_spec.num_threads
    run:
        from fossil_finder.blast.runner import BlastRunner
        runner = BlastRunner(blast_spec)
        runner.run(
            query=input.query,
            database=params.db_prefix,
            output=output.tsv,
        )


rule analyze:
    """Run full analysis pipeline on BLAST results."""
    input:
        blast_results=OUTPUT_DIR / "blast_results.tsv",
        regions_metadata=OUTPUT_DIR / "regions.tsv",
        annotation=ANNOTATION_GFF,
    output:
        summary=OUTPUT_DIR / "summary.json",
    script:
        "scripts/run_analysis.py"
```

- [ ] **Step 5: Run dry-run test**

Run: `pytest tests/test_workflow/test_snakemake.py::TestSnakemakeDryRun -v`
Expected: Tests that can run without BLAST+ should PASS (or skip if snakemake not installed)

- [ ] **Step 6: Run full test suite**

Run: `pytest tests/ --tb=short -q`
Expected: ALL PASS (workflow tests may skip if snakemake/BLAST+ not installed)

- [ ] **Step 7: Commit**

```bash
git add workflows/Snakefile workflows/scripts/run_analysis.py \
  tests/test_workflow/__init__.py tests/test_workflow/test_snakemake.py
git commit -m "feat: Snakemake workflow — extract → makeblastdb → BLAST → analyze"
```

---

## Chunk 2: Version Bump + Documentation

### Task 2: Version bump and dependency sync

**Files:**
- Modify: `src/fossil_finder/__init__.py`
- Modify: `pyproject.toml`
- Modify: `requirements.txt`

- [ ] **Step 1: Bump version**

Update `src/fossil_finder/__init__.py`:

```python
"""Fossil Finder: Multi-genome TE fossil mining and regulatory analysis framework."""

__version__ = "0.8.0"
```

Update `pyproject.toml` version to match:

```toml
version = "0.8.0"
```

Update `requirements.txt` snakemake pin (line 8):

```
snakemake>=8.0
```

- [ ] **Step 2: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

- [ ] **Step 3: Verify Snakemake dry-run works**

Run: `cd workflows && snakemake --dry-run --config genome_config=../tests/data/mini_genome_config.yaml output_dir=/tmp/test_output base_dir=../tests/data`
Expected: Shows planned execution of 4 rules (extract_regions, make_blast_db, run_blast, analyze)

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/__init__.py pyproject.toml requirements.txt
git commit -m "chore: bump version to 0.8.0 for Phase 8"
```

---

## Summary

| Task | What | Tests |
|------|------|-------|
| 1 | Snakefile + bridge script + workflow tests | 6 (3 dry-run, 2 bridge, 1 integration) |
| 2 | Version bump to 0.8.0 | 0 |

**Total new tests:** ~6 (most skip gracefully without Snakemake/BLAST+)

**Snakemake DAG:**
```
extract_regions ──→ run_blast ──→ analyze ──→ summary.json
                       ↑
make_blast_db ─────────┘
```

**What this enables:**
- `snakemake --cores 4 --config genome_config=dmel.yaml` — runs the full pipeline
- `snakemake --cores 4 --forcerun analyze` — re-runs only analysis (reuses BLAST cache)
- `snakemake --profile slurm` — runs on an HPC cluster
- `snakemake --dry-run` — shows what would run without executing
- File-based checkpointing — if BLAST completes but analysis fails, only analysis re-runs

**What's NOT included:**
- Shuffled-sequence control rules (future)
- Multi-genome comparison rules (future)
- Report generation rule (depends on Phase 7)
- Cluster profile templates (user-provided)
