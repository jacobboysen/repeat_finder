# Snakemake Workflow Design — Phase 6

## Goal

Replace the placeholder Snakefile with a production Snakemake workflow that orchestrates the full TE fossil pipeline — extraction, BLAST database creation, BLAST search, RepeatMasker annotation, and analysis — using the `fossil_finder` package as the computational backend. Support single-genome and multi-genome invocations with the same workflow.

## Architecture

A single `Snakefile` with 5 rules forming a DAG. Config is loaded from the same `GenomeConfig` YAML the Python package uses. Each rule produces file-based outputs, enabling Snakemake's automatic resume-from-failure. A bridge script connects the Snakemake file-based world with `PipelineRunner.analyze()`.

Multi-genome support is achieved via a `{genome}` wildcard — Snakemake expands the DAG per genome automatically. Single-genome invocations are the same code path (list of one).

**Prerequisite:** The `fossil_finder` package must be installed (`uv pip install -e .`). Rules import from the package directly — no `sys.path` manipulation. This is already documented in the project's environment setup.

## DAG

```
extract_regions ──→ run_blast ──→ analyze ──→ summary.json
                       ↑              ↑
make_blast_db ─────────┘              │
                                      │
run_repeatmasker (conditional) ───────┘
```

Per genome. Multi-genome = the DAG is expanded N times.

## Config Interface

```bash
# Single genome
snakemake --cores 4 \
  --config genome_config=config/dmel_r6.66.yaml output_dir=results/

# Multi-genome
snakemake --cores 8 \
  --config genome_configs='["config/dmel.yaml","config/dsim.yaml"]' \
  output_dir=results/

# HPC cluster (user-provided profile)
snakemake --profile slurm --config genome_config=config/dmel.yaml
```

### Config keys

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `genome_config` | Yes (unless `genome_configs`) | — | Path to single GenomeConfig YAML |
| `genome_configs` | No | — | JSON list of GenomeConfig YAML paths |
| `output_dir` | No | `"output"` | Base output directory |
| `base_dir` | No | None | Base dir for resolving relative paths in config |
| `feature_type` | No | `"three_prime_UTR"` | GFF3 feature type to extract |
| `max_evalue` | No | None | E-value filter for analysis |
| `min_pident` | No | None | Percent identity filter for analysis |
| `min_length` | No | None | Alignment length filter for analysis |

### Output directory structure

Single genome: `{output_dir}/{assembly}/` (assembly name from config)
Multi-genome: `{output_dir}/{assembly1}/`, `{output_dir}/{assembly2}/`, etc.

Both modes use the same subdirectory structure — single-genome is just multi-genome with one item.

## Rules

### 1. `extract_regions`

Extracts query regions from the genome FASTA using GFF3 annotation.

- **Input:** genome FASTA, annotation GFF3
- **Output:** `{genome}/regions.fa`, `{genome}/regions.tsv`
- **Implementation:** `run:` block calling `step_extract_regions()`
- **Note:** Does not pass `force=` — defaults to `False`, which is correct for Snakemake (let Snakemake handle re-execution via file timestamps).

### 2. `make_blast_db`

Builds a BLAST nucleotide database from TE consensus/instance sequences.

- **Input:** TE FASTA (from `config.source.te_consensus` or `config.source.te_instances`)
- **Output:** `{genome}/blastdb/te_db.{nhr,nin,nsq}`
- **Implementation:** `shell:` calling `makeblastdb`
- **Note:** Each genome gets its own DB copy even if the TE source is shared. This is a minor inefficiency traded for simpler wildcard logic.

### 3. `run_blast`

Runs BLAST search of extracted regions against the TE database.

- **Input:** `regions.fa`, blast database
- **Output:** `{genome}/blast_results.tsv`
- **Implementation:** `run:` block calling `BlastRunner.run()`
- **Threads:** `blast_spec.num_threads` via Snakemake `threads:` directive

### 4. `run_repeatmasker` (conditional)

Runs RepeatMasker de novo on the genome FASTA. Only fires when no pre-existing `.out` file is configured.

- **Input:** genome FASTA
- **Output:** `{genome}/repeatmasker/{genome_fasta_name}.out`
- **Implementation:** `run:` block calling `RepeatMaskerRunner.run()`
- **Condition:** See "RepeatMasker Resolution" below.

### 5. `analyze`

Runs the full analysis pipeline on BLAST results.

- **Input:** `blast_results.tsv`, `regions.tsv`, RepeatMasker `.out` (resolved conditionally)
- **Output:** `summary.json` (plus all analysis TSVs/JSONs written by `PipelineRunner._save_results()`)
- **Implementation:** `script:` directive calling `workflows/scripts/run_analysis.py`
- **Params:** Per-genome config path and base_dir passed via `params:` (because `script:` runs in a separate Python process and cannot access Snakefile-scope variables like `GENOMES`)

## Bridge Script (`workflows/scripts/run_analysis.py`)

The bridge script connects Snakemake's file-based world to `PipelineRunner.analyze()`.

**Import strategy:** Imports `fossil_finder` as a normal installed package. No `sys.path` manipulation — `uv pip install -e .` is a prerequisite.

**Snakemake integration:** When called via `script:` directive, the `snakemake` object is injected into globals. The script reads:
- `snakemake.input.*` for file paths (blast_results, regions_metadata, rm_out)
- `snakemake.output.*` for output paths
- `snakemake.params.*` for per-genome config path
- `snakemake.config` for global settings (feature_type, filters)

**Responsibilities:**
1. Load `GenomeConfig` from the per-genome config path (passed via `params.genome_config`)
2. Build `query_to_gene` mapping from GFF3 annotation
3. Load `query_regions` DataFrame from `regions.tsv`
4. Load gene sets from `config.gene_sets` — read each gene list file and convert to `dict[str, set[str]]`
5. Call `PipelineRunner.analyze()` with all inputs, **always passing `repeatmasker_path`** explicitly from `snakemake.input.rm_out` (bypasses the internal RM fallback)

**Standalone mode:** Also works via argparse for testing without Snakemake.

## RepeatMasker Resolution

The `analyze` rule uses a Python input function to resolve its RM `.out` input:

```python
def resolve_rm_input(wildcards):
    cfg = GENOMES[wildcards.genome]
    rm_path = cfg.source.repeatmasker_out
    if rm_path:
        resolved = Path(rm_path)
        if BASE_DIR and not resolved.is_absolute():
            resolved = BASE_DIR / resolved
        if resolved.exists():
            return str(resolved)
    # Fall through: request the output of run_repeatmasker rule
    genome_fasta_name = Path(cfg.source.genome_fasta).name
    return f"{OUTPUT_DIR}/{wildcards.genome}/repeatmasker/{genome_fasta_name}.out"
```

**Behavior:**
1. `config.source.repeatmasker_out` is set AND the file exists on disk → use it directly. `run_repeatmasker` rule never fires.
2. Path is set but file doesn't exist → treat as missing, fall through to de novo run (user typo will surface as a clear Snakemake "missing input" error if RM isn't installed).
3. Path is not set → `run_repeatmasker` rule fires on the genome FASTA.

**Relationship to PipelineRunner.analyze() internal RM logic:** The bridge script always passes `repeatmasker_path=` explicitly, so the internal 3-tier fallback in `PipelineRunner.analyze()` is bypassed. The Snakemake DAG owns RM resolution. The internal fallback remains for programmatic users who call the Python API directly.

## Multi-Genome Implementation

At Snakefile parse time:

1. Accept `genome_config` (single string) or `genome_configs` (JSON list of strings)
2. If single, wrap in a list
3. For each config path, load it via `load_genome_config()` and extract `genome.assembly` as the genome key
4. Store two dicts:
   - `GENOMES: dict[str, GenomeConfig]` — assembly → loaded config
   - `GENOME_CONFIG_PATHS: dict[str, str]` — assembly → original YAML path
5. All rules use `{genome}` wildcard, resolved from `GENOMES.keys()`
6. `rule all` uses `expand("{output_dir}/{genome}/summary.json", genome=GENOMES.keys())`
7. Rules that need per-genome config access (analyze) pass the config path via `params:`:
   ```python
   params:
       genome_config=lambda wc: GENOME_CONFIG_PATHS[wc.genome]
   ```

Single-genome invocations are identical — one item in the list, one wildcard value.

## File Structure

```
workflows/
├── Snakefile                       # Config parsing, 5 rules, input functions (~120 lines)
├── scripts/
│   └── run_analysis.py             # Bridge: Snakemake → PipelineRunner.analyze()
tests/
├── test_workflow/
│   ├── __init__.py
│   └── test_snakemake.py           # Dry-run, bridge validation, integration tests
```

## Testing Strategy

Tests skip gracefully when Snakemake or BLAST+ are not installed.

| Test | Requires | What it validates |
|------|----------|-------------------|
| `test_snakefile_parses` | Snakemake | Snakefile loads without syntax errors |
| `test_dag_has_expected_rules` | Snakemake | DAG contains all 5 rules |
| `test_config_validation` | Snakemake | Fails gracefully with missing config |
| `test_bridge_script_exists` | — | Bridge script is present |
| `test_bridge_script_valid_python` | — | Bridge script parses without errors |
| `test_full_pipeline_mini_genome` | Snakemake + BLAST+ | End-to-end on synthetic test data |

## What's NOT Included (Deferred)

- **Shuffled-sequence controls** — will be added as additional rules later (shuffle → BLAST → compare). The DAG grows without touching existing rules.
- **Report generation** — depends on Phase 7 (CLI + reporting)
- **NTv3 scoring rule** — depends on Phase 8
- **Cluster profile templates** — users provide their own via `--profile`. Snakemake docs cover this well.

## Dependencies

- Snakemake >=8.0 (already in `pyproject.toml` optional deps)
- `fossil_finder` package installed (`uv pip install -e ".[workflow]"`)
- BLAST+ (for `makeblastdb` and `blastn`)
- RepeatMasker (optional — only needed if no pre-existing `.out`)
