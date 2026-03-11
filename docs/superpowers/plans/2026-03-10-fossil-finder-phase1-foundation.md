# Fossil Finder Phase 1: Package Foundation

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create the installable `fossil_finder` Python package with config-driven genome specification, generic I/O modules, and a synthetic test genome — the foundation everything else builds on.

**Architecture:** Replace the flat `scripts/` + `sys.path.insert` approach with a proper `src/fossil_finder/` package using modern Python packaging (`pyproject.toml`). All organism-specific values move into YAML config validated by Pydantic. Generic I/O modules (BLAST, FASTA, GFF3) are ported first since they have zero Dmel coupling. BED I/O deferred to Phase 3 (region extraction) where it's first needed.

**Tech Stack:** Python 3.11, Pydantic v2, pytest, pyproject.toml (PEP 621), conda/pip hybrid install

**Phases Overview:**
- **Phase 1** (this plan): Package skeleton, config schema, generic I/O, test fixtures
- **Phase 2**: Adapter layer (FlyBase, Ensembl, Custom) + organism-specific parsing
- **Phase 3**: Region extraction (generic from GFF + genome)
- **Phase 4**: BLAST search + deduplication
- **Phase 5**: Core analysis modules (enrichment, strand, density, families)
- **Phase 6**: Snakemake workflow DAG
- **Phase 7**: Functional integration + reporting
- **Phase 8**: NTv3 scoring pipeline
- **Phase 9**: Legacy validation + archive old scripts

---

## File Structure (Phase 1)

```
fossil_finder/                        # Project root (renamed from repeat_finder later or aliased)
├── pyproject.toml                    # Package definition (replaces setup.py)
├── environment.yml                   # Updated conda env (Python 3.11)
├── src/
│   └── fossil_finder/
│       ├── __init__.py               # Package version + top-level exports
│       ├── config/
│       │   ├── __init__.py
│       │   ├── schema.py             # Pydantic models for genome config
│       │   └── genomes/
│       │       ├── dmel_r6.66.yaml   # Drosophila melanogaster config
│       │       └── _template.yaml    # Annotated template for new genomes
│       └── io/
│           ├── __init__.py
│           ├── blast.py              # BLAST I/O (port of blast_io.py)
│           ├── fasta.py              # Generic FASTA parsing
│           └── gff.py                # GFF3 parsing (new, standards-based)
├── tests/
│   ├── conftest.py                   # Shared fixtures
│   ├── data/                         # Test fixture data
│   │   ├── mini_genome.fasta         # Synthetic ~1kb genome (2 chromosomes)
│   │   ├── mini_tes.fasta            # 3 synthetic TE consensus seqs
│   │   ├── mini_annotation.gff3      # Gene models for mini genome
│   │   ├── mini_blast_results.tsv    # Pre-computed BLAST output
│   │   └── mini_genome_config.yaml   # Config for test genome
│   ├── test_config/
│   │   └── test_schema.py
│   └── test_io/
│       ├── test_blast.py
│       ├── test_fasta.py
│       └── test_gff.py
└── scripts/                          # Existing scripts (untouched in Phase 1)
```

**Key decisions:**
- `src/` layout prevents accidental imports of uninstalled code
- `fossil_finder` as package name (importable as `from fossil_finder.io import blast`)
- Existing `scripts/` untouched — old code keeps working throughout refactor
- `pyproject.toml` replaces `setup.py` (PEP 621 standard)

---

## Chunk 1: Package Skeleton + Python Upgrade

### Task 1: Create pyproject.toml and package structure

**Files:**
- Create: `pyproject.toml`
- Create: `src/fossil_finder/__init__.py`
- Create: `src/fossil_finder/config/__init__.py`
- Create: `src/fossil_finder/io/__init__.py`
- Modify: `environment.yml`
- Delete (later): `setup.py` (replaced by pyproject.toml)

**Context:** The existing `setup.py` references a `src/` layout but `src/` doesn't exist. We create the real package structure. The old `setup.py` pointed to `bioinformatics-program` — we replace it entirely.

- [ ] **Step 1: Create the src directory structure**

```bash
mkdir -p src/fossil_finder/config/genomes
mkdir -p src/fossil_finder/io
```

- [ ] **Step 2: Create `pyproject.toml`**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "fossil-finder"
version = "0.1.0"
description = "Multi-genome TE fossil mining and regulatory analysis framework"
requires-python = ">=3.11"
license = "MIT"
dependencies = [
    "numpy>=1.24",
    "pandas>=2.0",
    "biopython>=1.82",
    "pydantic>=2.0",
    "pyyaml>=6.0",
    "requests>=2.28",
    "tqdm>=4.64",
]

[project.optional-dependencies]
analysis = [
    "scikit-learn>=1.3",
    "scipy>=1.10",
    "statsmodels>=0.14",
    "matplotlib>=3.7",
    "seaborn>=0.12",
]
scoring = [
    "torch>=2.0",
    "transformers>=4.55",
    "pyfaidx>=0.7",
]
workflow = [
    "snakemake>=8.0",
]
dev = [
    "pytest>=7.0",
    "pytest-cov>=4.0",
]
all = ["fossil-finder[analysis,scoring,workflow,dev]"]

[tool.hatch.build.targets.wheel]
packages = ["src/fossil_finder"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
```

Note: Core deps are minimal (numpy, pandas, biopython, pydantic, pyyaml, requests, tqdm). Heavy deps (torch, snakemake, matplotlib) are optional extras. This means `pip install fossil-finder` works on any machine; `pip install fossil-finder[all]` gets everything.

- [ ] **Step 3: Create `src/fossil_finder/__init__.py`**

```python
"""Fossil Finder: Multi-genome TE fossil mining and regulatory analysis framework."""

__version__ = "0.1.0"
```

- [ ] **Step 4: Create empty `__init__.py` files for subpackages**

```python
# src/fossil_finder/config/__init__.py
# src/fossil_finder/io/__init__.py
# Both empty for now — populated as modules are added.
```

- [ ] **Step 5: Update `environment.yml` for Python 3.11**

```yaml
name: fossil-finder
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - numpy>=1.24
  - pandas>=2.0
  - biopython>=1.82
  - scikit-learn>=1.3
  - matplotlib>=3.7
  - seaborn>=0.12
  - jupyter
  - snakemake>=8.0
  - blast
  - pyyaml>=6.0
  - pip
  - pip:
      - pydantic>=2.0
      - requests>=2.28
      - tqdm>=4.64
      - pytest>=7.0
      - pytest-cov>=4.0
```

Note: Environment name changes from `bioinformatics-program` to `fossil-finder`. Old env still works for legacy scripts.

- [ ] **Step 6: Verify the package installs**

```bash
cd /Users/jacobboysen/git_repos/repeat_finder
pip install -e ".[dev]"
python -c "import fossil_finder; print(fossil_finder.__version__)"
```

Expected: prints `0.1.0`

- [ ] **Step 7: Verify old scripts still work**

```bash
cd /Users/jacobboysen/git_repos/repeat_finder
python -c "import sys; sys.path.insert(0, 'scripts'); from utils.blast_io import BLAST_COLUMNS; print(len(BLAST_COLUMNS))"
```

Expected: prints `17`. Old import path still works — we haven't broken anything.

- [ ] **Step 8: Rename old setup.py to avoid conflicts**

```bash
mv setup.py setup.py.legacy
```

The old `setup.py` defined package `bioinformatics-program` with `src/` layout but no `src/` existed. `pyproject.toml` now owns the package definition. Keeping the file as `.legacy` preserves history without confusing build tools.

- [ ] **Step 9: Commit**

```bash
git add pyproject.toml src/ environment.yml setup.py.legacy
git rm setup.py
git commit -m "feat: initialize fossil_finder package with pyproject.toml and Python 3.11

Replaces setup.py with modern PEP 621 pyproject.toml.
Package installs as fossil-finder with optional dependency groups.
Old scripts/ untouched — both import paths coexist during migration.
Old setup.py preserved as setup.py.legacy for reference."
```

---

### Task 2: Create test infrastructure

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/__init__.py`
- Create: `tests/test_config/__init__.py`
- Create: `tests/test_io/__init__.py`

- [ ] **Step 1: Create test directory structure**

```bash
mkdir -p tests/test_config tests/test_io tests/data
```

- [ ] **Step 2: Create `tests/conftest.py` with path fixtures**

```python
"""Shared test fixtures for fossil_finder."""

from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir() -> Path:
    """Path to test fixture data."""
    return Path(__file__).parent / "data"


@pytest.fixture
def mini_genome_fasta(test_data_dir) -> Path:
    """Path to synthetic mini genome FASTA."""
    return test_data_dir / "mini_genome.fasta"


@pytest.fixture
def mini_tes_fasta(test_data_dir) -> Path:
    """Path to synthetic mini TE consensus FASTA."""
    return test_data_dir / "mini_tes.fasta"


@pytest.fixture
def mini_annotation_gff(test_data_dir) -> Path:
    """Path to synthetic mini genome GFF3 annotation."""
    return test_data_dir / "mini_annotation.gff3"


@pytest.fixture
def mini_blast_results(test_data_dir) -> Path:
    """Path to pre-computed BLAST results for mini genome."""
    return test_data_dir / "mini_blast_results.tsv"


@pytest.fixture
def mini_genome_config(test_data_dir) -> Path:
    """Path to genome config YAML for mini test genome."""
    return test_data_dir / "mini_genome_config.yaml"
```

- [ ] **Step 3: Create empty `__init__.py` files**

```python
# tests/__init__.py
# tests/test_config/__init__.py
# tests/test_io/__init__.py
# All empty.
```

- [ ] **Step 4: Verify pytest discovers tests**

```bash
pytest --collect-only 2>&1 | head -5
```

Expected: `no tests ran` (no test files yet, but no errors)

- [ ] **Step 5: Commit**

```bash
git add tests/
git commit -m "test: add pytest infrastructure with fixture paths"
```

---

## Chunk 2: Genome Config Schema

### Task 3: Define Pydantic config models

**Files:**
- Create: `src/fossil_finder/config/schema.py`
- Test: `tests/test_config/test_schema.py`

**Context:** This is the central abstraction — every organism-specific value that was hardcoded across 50+ files becomes a field in this schema. The schema validates config at load time so errors surface immediately, not 3 hours into a pipeline run.

- [ ] **Step 1: Write the failing test for config loading**

`tests/test_config/test_schema.py`:

```python
"""Tests for genome configuration schema."""

import pytest
import yaml

from fossil_finder.config.schema import (
    GenomeConfig,
    load_genome_config,
)


MINIMAL_CONFIG = {
    "genome": {
        "species": "Testus synthetica",
        "assembly": "test_v1",
    },
    "source": {
        "adapter": "custom",
        "genome_fasta": "/path/to/genome.fasta",
        "annotation_gff": "/path/to/annotation.gff3",
    },
}


class TestGenomeConfig:
    def test_minimal_config_loads(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.genome.species == "Testus synthetica"
        assert config.genome.assembly == "test_v1"
        assert config.source.adapter == "custom"

    def test_missing_required_fields_raises(self):
        from pydantic import ValidationError
        with pytest.raises(ValidationError):
            GenomeConfig(**{"genome": {"species": "Test"}})

    def test_chromosomes_default_to_auto(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.genome.chromosomes is None  # None means auto-detect

    def test_optional_fields_default_to_none(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.source.te_instances is None
        assert config.source.phylop_bigwig is None
        assert config.source.synteny is None

    def test_blast_params_have_defaults(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.blast.word_size == 7
        assert config.blast.dust is False
        assert config.blast.evalue == 0.001

    def test_gene_sets_optional(self):
        config = GenomeConfig(**MINIMAL_CONFIG)
        assert config.gene_sets == {}

    def test_full_config_with_gene_sets(self):
        full = {
            **MINIMAL_CONFIG,
            "gene_sets": {
                "target_genes": {
                    "description": "Genes of interest",
                    "genes": "/path/to/gene_ids.txt",
                }
            },
        }
        config = GenomeConfig(**full)
        assert "target_genes" in config.gene_sets
        assert config.gene_sets["target_genes"].description == "Genes of interest"

    def test_adapter_must_be_valid_choice(self):
        from pydantic import ValidationError
        bad = {**MINIMAL_CONFIG}
        bad["source"] = {**bad["source"], "adapter": "nonexistent_db"}
        with pytest.raises(ValidationError):
            GenomeConfig(**bad)


class TestLoadGenomeConfig:
    def test_load_from_yaml_file(self, tmp_path):
        config_path = tmp_path / "test_config.yaml"
        config_path.write_text(yaml.dump(MINIMAL_CONFIG))
        config = load_genome_config(config_path)
        assert config.genome.species == "Testus synthetica"

    def test_load_nonexistent_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_genome_config("/nonexistent/config.yaml")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
pytest tests/test_config/test_schema.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.config.schema'`

- [ ] **Step 3: Implement the config schema**

`src/fossil_finder/config/schema.py`:

```python
"""Genome configuration schema for fossil_finder.

Every organism-specific value lives here. Scripts never hardcode
species names, chromosome lists, file paths, or ID formats.
"""

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, model_validator


class GenomeSpec(BaseModel):
    """Core genome identity and structure."""

    species: str
    assembly: str
    release: str | None = None
    chromosomes: list[str] | None = None  # None = auto-detect from FASTA index
    ucsc_prefix: str | None = None  # e.g. "chr" for UCSC bigWig queries


class SyntenySpec(BaseModel):
    """Synteny/conservation comparison species."""

    maf_dir: str | None = None
    species: list[str] = []


class FunctionalSpec(BaseModel):
    """Optional functional annotation data sources."""

    expression: str | None = None
    go_annotations: str | None = None
    localization: str | None = None
    gene_groups: str | None = None


class SourceSpec(BaseModel):
    """Data source configuration — where reference files come from."""

    adapter: Literal["flybase", "ensembl", "ncbi", "custom"]
    genome_fasta: str
    annotation_gff: str
    gene_id_prefix: str | None = None  # e.g. "FBgn", "ENSMUSG"
    transcript_id_prefix: str | None = None

    # TE sources (at least one should be provided for TE analysis)
    te_instances: str | None = None
    te_consensus: str | None = None  # "dfam:Species_name" or path

    # Optional organism-specific data
    gene_symbol_map: str | None = None
    repeatmasker_out: str | None = None
    phylop_bigwig: str | None = None
    synteny: SyntenySpec | None = None
    functional: FunctionalSpec | None = None


class BlastSpec(BaseModel):
    """BLAST search parameters.

    Defaults based on DUST_FILTERING_ANALYSIS.md (2026-01-22):
    dust=no with stringent e-value captures 52% more high-quality hits
    and higher real/shuffled enrichment ratios than dust=yes.
    """

    word_size: int = 7
    gapopen: int = 2
    gapextend: int = 1
    penalty: int = -1
    reward: int = 1
    dust: bool = False
    evalue: float = 0.001
    max_target_seqs: int = 1000
    num_threads: int = 4


class ScoringSpec(BaseModel):
    """ML/LM scoring configuration."""

    model: str = "InstaDeepAI/NTv3_100M_post"
    species_key: str | None = None  # auto-detected if None
    window_size: int = 32768


class GeneSetSpec(BaseModel):
    """A named gene set for comparative analysis."""

    description: str = ""
    genes: str  # path to gene ID list file
    tier: int = 1


class GenomeConfig(BaseModel):
    """Top-level genome configuration.

    Validated at load time — errors surface before the pipeline runs,
    not 3 hours in.
    """

    genome: GenomeSpec
    source: SourceSpec
    blast: BlastSpec = BlastSpec()
    scoring: ScoringSpec = ScoringSpec()
    gene_sets: dict[str, GeneSetSpec] = {}

    @model_validator(mode="after")
    def warn_no_te_source(self):
        if not self.source.te_instances and not self.source.te_consensus:
            import warnings
            warnings.warn(
                "No TE source configured (te_instances or te_consensus). "
                "TE fossil analysis requires at least one TE database.",
                UserWarning,
                stacklevel=2,
            )
        return self


def load_genome_config(path: str | Path) -> GenomeConfig:
    """Load and validate a genome config from YAML file.

    Args:
        path: Path to genome config YAML.

    Returns:
        Validated GenomeConfig instance.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        pydantic.ValidationError: If config is invalid.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path) as f:
        raw = yaml.safe_load(f)

    return GenomeConfig(**raw)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_config/test_schema.py -v
```

Expected: All 10 tests PASS (8 in TestGenomeConfig + 2 in TestLoadGenomeConfig).

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/config/schema.py tests/test_config/test_schema.py
git commit -m "feat: genome config schema with Pydantic validation

Config-driven genome specification replaces all hardcoded species values.
Supports flybase/ensembl/ncbi/custom adapters with optional data sources.
BLAST params default to optimized values (word_size=7, dust=yes)."
```

---

### Task 4: Create Dmel config and annotated template

**Files:**
- Create: `src/fossil_finder/config/genomes/dmel_r6.66.yaml`
- Create: `src/fossil_finder/config/genomes/_template.yaml`
- Test: extend `tests/test_config/test_schema.py`

- [ ] **Step 1: Write test that Dmel config validates**

Append to `tests/test_config/test_schema.py`:

```python
class TestDmelConfig:
    def test_dmel_config_validates(self):
        dmel_path = (
            Path(__file__).parent.parent.parent
            / "src" / "fossil_finder" / "config" / "genomes" / "dmel_r6.66.yaml"
        )
        config = load_genome_config(dmel_path)
        assert config.genome.species == "Drosophila melanogaster"
        assert config.genome.assembly == "dm6"
        assert config.source.adapter == "flybase"
        assert "2L" in config.genome.chromosomes
        assert config.blast.word_size == 7
        assert config.blast.dust is False
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_config/test_schema.py::TestDmelConfig -v
```

Expected: FAIL — file not found.

- [ ] **Step 3: Create `src/fossil_finder/config/genomes/dmel_r6.66.yaml`**

```yaml
# Drosophila melanogaster (dm6 / FlyBase r6.66)
# Reference configuration for the original TE fossil pipeline.

genome:
  species: "Drosophila melanogaster"
  assembly: "dm6"
  release: "r6.66"
  chromosomes: ["2L", "2R", "3L", "3R", "X", "4"]
  ucsc_prefix: "chr"

source:
  adapter: "flybase"

  # Reference files (relative to project data/ directory)
  genome_fasta: "references/dmel_genome.fasta"
  annotation_gff: "references/dmel-all-r6.66.gff"
  gene_id_prefix: "FBgn"
  transcript_id_prefix: "FBtr"

  # TE databases
  te_instances: "references/dmel_te_flybase.fasta"
  te_consensus: "references/dmel_te_consensus.fasta"

  # Organism-specific data
  gene_symbol_map: "references/fbgn_to_symbol.tsv"
  repeatmasker_out: "references/dm6.fa.out"
  phylop_bigwig: "references/dm6.phyloP27way.bw"

  synteny:
    maf_dir: "references/maf/"
    species: ["droSim1", "droYak2", "droEre2", "droSec1"]

  functional:
    expression: "annotations/raw/gene_rpkm_report.tsv"
    go_annotations: "annotations/raw/gene_association.gaf"
    localization: "annotations/raw/flyfish_localization.csv"
    gene_groups: "annotations/raw/gene_group_data.tsv"

blast:
  word_size: 7
  gapopen: 2
  gapextend: 1
  penalty: -1
  reward: 1
  dust: false
  evalue: 0.001

scoring:
  model: "InstaDeepAI/NTv3_100M_post"
  species_key: "fruit_fly"
  window_size: 32768

gene_sets:
  germ_plasm:
    description: "Germ plasm-localized mRNAs"
    genes: "gene_lists/germ_plasm_fbgn_ids.txt"
    tier: 1
  housekeeping:
    description: "Housekeeping controls"
    genes: "gene_lists/housekeeping_fbgn_ids.txt"
    tier: 2
  somatic:
    description: "Somatically-localized genes"
    genes: "gene_lists/somatic_fbgn_ids.txt"
    tier: 2
  cleared:
    description: "Posteriorly-cleared transcripts"
    genes: "gene_lists/cleared_fbgn_ids.txt"
    tier: 2
  adult:
    description: "Adult-expressed genes"
    genes: "gene_lists/adult_fbgn_ids.txt"
    tier: 2
```

- [ ] **Step 4: Create `src/fossil_finder/config/genomes/_template.yaml`**

```yaml
# Genome Configuration Template for fossil_finder
#
# Copy this file, rename to <species>_<assembly>.yaml, and fill in your values.
# Fields marked REQUIRED must be provided. All others are optional and will
# gracefully degrade (analyses that need missing data will be skipped).
#
# Paths are relative to the project's data/ directory unless absolute.

genome:
  species: ""              # REQUIRED. Binomial name, e.g. "Apis mellifera"
  assembly: ""             # REQUIRED. Assembly ID, e.g. "Amel_HAv3.1"
  release: null            # Database release version, e.g. "r6.66", "112"
  chromosomes: null        # List of chromosome names to analyze, or null for auto-detect
                           # Example: ["1", "2", "3", "X"]
  ucsc_prefix: null        # Prefix for UCSC tool compatibility, e.g. "chr"

source:
  adapter: "custom"        # REQUIRED. One of: flybase, ensembl, ncbi, custom
                           #   flybase  - FlyBase GFF + ID formats (Drosophila)
                           #   ensembl  - Ensembl GFF3 + BioMart (broad species)
                           #   ncbi     - NCBI RefSeq GFF3 + gene IDs
                           #   custom   - Standard GFF3 with minimal assumptions

  genome_fasta: ""         # REQUIRED. Genome assembly FASTA (local path or URL)
  annotation_gff: ""       # REQUIRED. Gene annotations in GFF3 format
  gene_id_prefix: null     # Gene ID prefix for validation, e.g. "FBgn", "ENSMUSG"
  transcript_id_prefix: null

  # TE databases (provide at least one for TE analysis)
  te_instances: null       # Known TE instances FASTA (if available for your species)
  te_consensus: null       # TE consensus library. Can be:
                           #   - Local FASTA path
                           #   - "dfam:Species_name" to query Dfam
                           #   - null (run de novo TE discovery later)

  # Optional reference data
  gene_symbol_map: null    # TSV mapping gene IDs to symbols
  repeatmasker_out: null   # RepeatMasker .out file for this assembly
  phylop_bigwig: null      # phyloP conservation scores (BigWig format)

  synteny: null            # Set to enable cross-species conservation analysis:
  #  maf_dir: "path/to/maf/"
  #  species: ["species1", "species2"]

  functional: null         # Set to enable functional data integration:
  #  expression: "path/to/rpkm_or_tpm.tsv"
  #  go_annotations: "path/to/gene_association.gaf"
  #  localization: null    # Spatial expression data (organism-specific)
  #  gene_groups: null     # Gene group definitions

# BLAST parameters (defaults are optimized for diverged TE detection)
# dust=false with strict evalue per DUST_FILTERING_ANALYSIS.md findings
blast:
  word_size: 7
  gapopen: 2
  gapextend: 1
  penalty: -1
  reward: 1
  dust: false
  evalue: 0.001

# Gene sets for comparative analysis (optional)
gene_sets: {}
#  my_target_genes:
#    description: "Genes of interest"
#    genes: "gene_lists/target_gene_ids.txt"
#    tier: 1
```

- [ ] **Step 5: Run test to verify it passes**

```bash
pytest tests/test_config/test_schema.py -v
```

Expected: All tests PASS (including new TestDmelConfig).

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/config/genomes/
git commit -m "feat: add Dmel config and annotated template for new genomes

dmel_r6.66.yaml encodes all previously-hardcoded values.
_template.yaml is the starting point for adding any new genome."
```

---

## Chunk 3: Generic I/O Modules

### Task 5: Port BLAST I/O (already 100% generic)

**Files:**
- Create: `src/fossil_finder/io/blast.py`
- Test: `tests/test_io/test_blast.py`
- Create: `tests/data/mini_blast_results.tsv`

**Context:** `scripts/utils/blast_io.py` is already fully generic — no Dmel coupling. This is a clean copy-and-modernize: update typing syntax, improve docstrings, keep identical behavior. The test validates against a small fixture file.

- [ ] **Step 1: Create BLAST test fixture**

`tests/data/mini_blast_results.tsv` — 5 representative BLAST hits (17-column TSV, no header):

```tsv
gene1_utr	TE_gypsy1	78.5	120	20	3	50	170	1	120	1.5e-10	85.2	500	5000	ATCGATCG	ATCGATCG	plus
gene1_utr	TE_copia1	65.3	85	25	4	200	284	300	216	3.2e-5	42.1	500	4500	GCTAGCTA	GCTAGCTA	minus
gene2_utr	TE_gypsy1	82.1	200	30	5	10	210	50	250	2.1e-15	120.5	800	5000	TTAATTAA	TTAATTAA	plus
gene2_utr	TE_pogo1	71.0	60	15	2	300	360	100	160	0.005	30.8	800	3000	AAGGCCTT	AAGGCCTT	plus
gene3_utr	TE_jockey1	90.2	150	12	1	100	250	1	150	8.3e-20	155.0	600	4000	CCGGAATT	CCGGAATT	plus
```

- [ ] **Step 2: Write failing tests**

`tests/test_io/test_blast.py`:

```python
"""Tests for BLAST I/O module."""

import pandas as pd
import pytest

from fossil_finder.io.blast import (
    BLAST_COLUMNS,
    classify_strand,
    load_blast_results,
    parse_blast_line,
    iter_blast_results,
    summarize_blast_results,
)


class TestClassifyStrand:
    def test_plus_strand(self):
        assert classify_strand(1, 100) == "plus"

    def test_minus_strand(self):
        assert classify_strand(100, 1) == "minus"

    def test_equal_positions(self):
        # Edge case: same position
        assert classify_strand(50, 50) == "minus"


class TestBlastColumns:
    def test_column_count(self):
        assert len(BLAST_COLUMNS) == 17

    def test_required_columns_present(self):
        for col in ["qseqid", "sseqid", "pident", "evalue", "bitscore"]:
            assert col in BLAST_COLUMNS


class TestLoadBlastResults:
    def test_load_17col_file(self, mini_blast_results):
        df = load_blast_results(mini_blast_results)
        assert len(df) == 5
        assert list(df.columns) == BLAST_COLUMNS

    def test_filter_by_evalue(self, mini_blast_results):
        df = load_blast_results(mini_blast_results, max_evalue=1e-10)
        assert all(df["evalue"] <= 1e-10)

    def test_filter_by_min_length(self, mini_blast_results):
        df = load_blast_results(mini_blast_results, min_length=100)
        assert all(df["length"] >= 100)

    def test_nonexistent_file_returns_empty(self, tmp_path):
        df = load_blast_results(tmp_path / "nonexistent.tsv")
        assert df.empty
        assert list(df.columns) == BLAST_COLUMNS

    def test_empty_file_returns_empty(self, tmp_path):
        empty = tmp_path / "empty.tsv"
        empty.write_text("")
        df = load_blast_results(empty)
        assert df.empty


class TestParseBlastLine:
    def test_parse_complete_line(self):
        line = "q1\ts1\t80.0\t100\t15\t2\t1\t100\t1\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG\tplus"
        result = parse_blast_line(line)
        assert result["qseqid"] == "q1"
        assert result["pident"] == 80.0
        assert result["length"] == 100
        assert result["strand"] == "plus"

    def test_parse_16col_line_adds_strand(self):
        line = "q1\ts1\t80.0\t100\t15\t2\t1\t100\t200\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG"
        result = parse_blast_line(line)
        assert result["strand"] == "minus"  # sstart=200 > send=100


class TestIterBlastResults:
    def test_iterates_all_lines(self, mini_blast_results):
        hits = list(iter_blast_results(mini_blast_results))
        assert len(hits) == 5

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            list(iter_blast_results(tmp_path / "nope.tsv"))


class TestSummarizeBlastResults:
    def test_summary_stats(self, mini_blast_results):
        df = load_blast_results(mini_blast_results)
        summary = summarize_blast_results(df)
        assert summary["total_hits"] == 5
        assert summary["unique_queries"] == 3
        assert summary["unique_subjects"] == 4
        assert summary["strand_plus"] == 4
        assert summary["strand_minus"] == 1

    def test_empty_dataframe_summary(self):
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        summary = summarize_blast_results(df)
        assert summary["total_hits"] == 0
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
pytest tests/test_io/test_blast.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.io.blast'`

- [ ] **Step 4: Port blast_io.py to new location**

`src/fossil_finder/io/blast.py`:

Copy the full contents of `scripts/utils/blast_io.py` with these changes:
1. Module docstring updated to reference `fossil_finder`
2. `from typing import Union` replaced with `str | Path` union syntax (Python 3.11)
3. All logic identical — this is a straight port

```python
"""BLAST I/O utilities for fossil_finder.

Provides consistent parsing of BLAST tabular output.
All TSV files use the 17-column format (see BLAST_COLUMNS).
"""

from pathlib import Path

import pandas as pd


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq", "strand",
]

BLAST_COLUMNS_NO_STRAND = BLAST_COLUMNS[:-1]


def classify_strand(sstart: int, send: int) -> str:
    """Classify hit strand based on subject coordinates.

    When sstart > send in BLAST output, the hit is on the minus strand.
    """
    return "plus" if sstart < send else "minus"


def load_blast_results(
    results_file: str | Path,
    add_strand: bool = True,
    min_length: int | None = None,
    min_pident: float | None = None,
    max_evalue: float | None = None,
) -> pd.DataFrame:
    """Load BLAST results from TSV file.

    Auto-detects 16 vs 17 column format. Adds strand classification
    if missing. Returns empty DataFrame (with correct columns) for
    missing or empty files.
    """
    results_file = Path(results_file)

    if not results_file.exists() or results_file.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    with open(results_file) as f:
        first_line = f.readline()
    num_cols = len(first_line.strip().split("\t"))

    if num_cols == 17:
        df = pd.read_csv(results_file, sep="\t", names=BLAST_COLUMNS)
    elif num_cols == 16:
        df = pd.read_csv(results_file, sep="\t", names=BLAST_COLUMNS_NO_STRAND)
        if add_strand:
            df["strand"] = df.apply(
                lambda row: classify_strand(row["sstart"], row["send"]), axis=1
            )
    else:
        basic_cols = BLAST_COLUMNS[: min(num_cols, len(BLAST_COLUMNS))]
        df = pd.read_csv(results_file, sep="\t", names=basic_cols)
        if add_strand and "sstart" in df.columns and "send" in df.columns:
            df["strand"] = df.apply(
                lambda row: classify_strand(row["sstart"], row["send"]), axis=1
            )

    if min_length is not None and "length" in df.columns:
        df = df[df["length"] >= min_length]
    if min_pident is not None and "pident" in df.columns:
        df = df[df["pident"] >= min_pident]
    if max_evalue is not None and "evalue" in df.columns:
        df = df[df["evalue"] <= max_evalue]

    return df


def parse_blast_line(line: str) -> dict:
    """Parse a single BLAST TSV line into a dictionary.

    Useful for streaming large files without loading into memory.
    """
    parts = line.strip().split("\t")
    result = {}

    for i, col in enumerate(BLAST_COLUMNS):
        if i >= len(parts):
            result[col] = None
            continue

        value = parts[i]
        if col in ("pident", "evalue", "bitscore"):
            result[col] = float(value) if value else 0.0
        elif col in (
            "length", "mismatch", "gapopen", "qstart", "qend",
            "sstart", "send", "qlen", "slen",
        ):
            result[col] = int(value) if value else 0
        else:
            result[col] = value

    if "strand" not in result or result["strand"] is None:
        if result.get("sstart") is not None and result.get("send") is not None:
            result["strand"] = classify_strand(result["sstart"], result["send"])

    return result


def iter_blast_results(results_file: str | Path):
    """Iterate over BLAST results line by line.

    Yields dict per hit. Useful for large files that don't fit in memory.
    Raises FileNotFoundError if file doesn't exist (fail-fast).
    """
    results_file = Path(results_file)
    if not results_file.exists():
        raise FileNotFoundError(f"BLAST results file not found: {results_file}")

    with open(results_file) as f:
        for line in f:
            if line.strip():
                yield parse_blast_line(line)


def summarize_blast_results(df: pd.DataFrame) -> dict:
    """Generate summary statistics for BLAST results DataFrame."""
    if df.empty:
        return {
            "total_hits": 0, "unique_queries": 0, "unique_subjects": 0,
            "mean_pident": 0, "mean_length": 0, "mean_evalue": 0,
            "strand_plus": 0, "strand_minus": 0,
        }

    summary = {
        "total_hits": len(df),
        "unique_queries": df["qseqid"].nunique() if "qseqid" in df.columns else 0,
        "unique_subjects": df["sseqid"].nunique() if "sseqid" in df.columns else 0,
        "mean_pident": df["pident"].mean() if "pident" in df.columns else 0,
        "mean_length": df["length"].mean() if "length" in df.columns else 0,
        "mean_evalue": df["evalue"].mean() if "evalue" in df.columns else 0,
    }

    if "strand" in df.columns:
        strand_counts = df["strand"].value_counts()
        summary["strand_plus"] = strand_counts.get("plus", 0)
        summary["strand_minus"] = strand_counts.get("minus", 0)
    else:
        summary["strand_plus"] = 0
        summary["strand_minus"] = 0

    return summary
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/test_io/test_blast.py -v
```

Expected: All 16 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/io/blast.py tests/test_io/test_blast.py tests/data/mini_blast_results.tsv
git commit -m "feat: port BLAST I/O to fossil_finder.io.blast

Direct port from scripts/utils/blast_io.py (already 100% generic).
Updated to Python 3.11 syntax (str | Path unions).
16 tests covering load, parse, iterate, summarize, and edge cases."
```

---

### Task 6: Create generic FASTA parser

**Files:**
- Create: `src/fossil_finder/io/fasta.py`
- Test: `tests/test_io/test_fasta.py`
- Create: `tests/data/mini_genome.fasta`
- Create: `tests/data/mini_tes.fasta`

**Context:** The existing `data_loaders.py` has two FASTA parsers: `parse_fasta()` (generic) and `parse_fasta_by_parent()` (FlyBase-specific header format). We port `parse_fasta()` and create a new, more flexible header parser that works with any key=value header format.

- [ ] **Step 1: Create test FASTA fixtures**

`tests/data/mini_genome.fasta` — A synthetic 1000bp "genome" with 2 chromosomes. Each chromosome has a known sequence with planted TE-like regions:

```fasta
>chr1 Synthetic chromosome 1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATATATATATATATATATATATATATATATATATATATATATATATATATAT
CCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCC
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
ATATATATGCGCGCGCATATATATGCGCGCGCATATATATGCGCGCGCATAT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
>chr2 Synthetic chromosome 2
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAAT
CCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATATATATCCCCGGGGATATATATCCCCGGGGATATATATCCCCGGGGATAT
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
```

`tests/data/mini_tes.fasta` — 3 synthetic TE consensus sequences:

```fasta
>TE_LTR1 class=LTR; family=gypsy; length=200
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATATATATATATATATATATATATATATATATATATATATATATATATATAT
CCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCC
>TE_DNA1 class=DNA; family=mariner; length=150
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAAT
CCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGG
>TE_LINE1 class=LINE; family=jockey; length=180
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATATATATGCGCGCGCATATATATGCGCGCGCATATATATGCGCGCGCATAT
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

- [ ] **Step 2: Write failing tests**

`tests/test_io/test_fasta.py`:

```python
"""Tests for FASTA I/O module."""

import pytest

from fossil_finder.io.fasta import (
    parse_fasta,
    parse_fasta_headers,
    iter_fasta,
    write_fasta,
)


class TestParseFasta:
    def test_parse_sequences(self, mini_genome_fasta):
        seqs = parse_fasta(mini_genome_fasta)
        assert len(seqs) == 2
        assert "chr1" in seqs
        assert "chr2" in seqs

    def test_sequences_are_strings(self, mini_genome_fasta):
        seqs = parse_fasta(mini_genome_fasta)
        for seq in seqs.values():
            assert isinstance(seq, str)
            assert "\n" not in seq  # No line breaks in sequence

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_fasta(tmp_path / "nope.fasta")

    def test_parse_te_fasta(self, mini_tes_fasta):
        seqs = parse_fasta(mini_tes_fasta)
        assert len(seqs) == 3
        assert "TE_LTR1" in seqs
        assert "TE_DNA1" in seqs
        assert "TE_LINE1" in seqs


class TestParseFastaHeaders:
    def test_extract_key_value_attributes(self, mini_tes_fasta):
        headers = parse_fasta_headers(mini_tes_fasta)
        assert headers["TE_LTR1"]["class"] == "LTR"
        assert headers["TE_LTR1"]["family"] == "gypsy"

    def test_raw_description_preserved(self, mini_tes_fasta):
        headers = parse_fasta_headers(mini_tes_fasta)
        assert "description" in headers["TE_LTR1"]

    def test_plain_headers(self, mini_genome_fasta):
        headers = parse_fasta_headers(mini_genome_fasta)
        # chr1 has description "Synthetic chromosome 1" but no key=value pairs
        assert headers["chr1"]["description"] == "Synthetic chromosome 1"


class TestIterFasta:
    def test_yields_id_and_sequence(self, mini_genome_fasta):
        records = list(iter_fasta(mini_genome_fasta))
        assert len(records) == 2
        seq_id, seq = records[0]
        assert seq_id == "chr1"
        assert isinstance(seq, str)


class TestWriteFasta:
    def test_roundtrip(self, mini_tes_fasta, tmp_path):
        seqs = parse_fasta(mini_tes_fasta)
        out_path = tmp_path / "out.fasta"
        write_fasta(seqs, out_path)
        reloaded = parse_fasta(out_path)
        assert seqs == reloaded

    def test_line_wrapping(self, tmp_path):
        seqs = {"test": "A" * 200}
        out_path = tmp_path / "wrapped.fasta"
        write_fasta(seqs, out_path, line_width=80)
        text = out_path.read_text()
        lines = text.strip().split("\n")
        assert lines[0] == ">test"
        assert len(lines[1]) == 80  # First sequence line wrapped at 80
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
pytest tests/test_io/test_fasta.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.io.fasta'`

- [ ] **Step 4: Implement generic FASTA module**

`src/fossil_finder/io/fasta.py`:

```python
"""Generic FASTA I/O for fossil_finder.

Parses any FASTA file without organism-specific assumptions.
Handles key=value header attributes flexibly.
"""

import re
from pathlib import Path
from collections.abc import Iterator


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Parse FASTA file into {sequence_id: sequence} dict.

    Sequence ID is the first whitespace-delimited token after '>'.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def parse_fasta_headers(path: str | Path) -> dict[str, dict]:
    """Parse FASTA headers, extracting key=value attributes.

    Returns {seq_id: {"description": "...", "key1": "val1", ...}}.
    Handles both 'key=value;' and 'key=value' formats.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    headers: dict[str, dict] = {}

    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            parts = line.strip()[1:].split(None, 1)
            seq_id = parts[0]
            raw_desc = parts[1] if len(parts) > 1 else ""

            attrs: dict[str, str] = {"description": raw_desc}

            # Extract key=value pairs (with optional trailing semicolons)
            for match in re.finditer(r"(\w+)=([^;\s]+)", raw_desc):
                attrs[match.group(1)] = match.group(2)

            headers[seq_id] = attrs

    return headers


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Iterate over FASTA records yielding (seq_id, sequence) tuples.

    Memory-efficient for large files — only one record in memory at a time.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    current_id: str | None = None
    current_seq: list[str] = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    yield current_id, "".join(current_seq)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            yield current_id, "".join(current_seq)


def write_fasta(
    sequences: dict[str, str],
    path: str | Path,
    line_width: int = 80,
) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: {seq_id: sequence} dict.
        path: Output file path.
        line_width: Characters per sequence line (default 80, FASTA standard).
    """
    path = Path(path)
    with open(path, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/test_io/test_fasta.py -v
```

Expected: All 10 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/io/fasta.py tests/test_io/test_fasta.py tests/data/mini_genome.fasta tests/data/mini_tes.fasta
git commit -m "feat: generic FASTA I/O with flexible header parsing

parse_fasta() for simple ID→sequence loading.
parse_fasta_headers() extracts key=value attrs from any header format.
iter_fasta() for memory-efficient streaming of large files.
write_fasta() for roundtrip output with configurable line width."
```

---

### Task 7: Create GFF3 parser

**Files:**
- Create: `src/fossil_finder/io/gff.py`
- Test: `tests/test_io/test_gff.py`
- Create: `tests/data/mini_annotation.gff3`

**Context:** This is NEW code, not a port. The existing codebase parses GFF ad-hoc in each extraction script with FlyBase-specific assumptions. We create a standards-compliant GFF3 parser that works with any GFF3 file (FlyBase, Ensembl, NCBI, custom).

- [ ] **Step 1: Create GFF3 test fixture**

`tests/data/mini_annotation.gff3`:

```gff3
##gff-version 3
##sequence-region chr1 1 520
##sequence-region chr2 1 520
chr1	test	gene	1	300	.	+	.	ID=gene001;Name=testgene1
chr1	test	mRNA	1	300	.	+	.	ID=mRNA001;Parent=gene001;Name=testgene1-RA
chr1	test	exon	1	100	.	+	.	ID=exon001;Parent=mRNA001
chr1	test	CDS	20	90	.	+	0	ID=cds001;Parent=mRNA001
chr1	test	exon	200	300	.	+	.	ID=exon002;Parent=mRNA001
chr1	test	CDS	200	280	.	+	2	ID=cds002;Parent=mRNA001
chr1	test	five_prime_UTR	1	19	.	+	.	ID=utr5_001;Parent=mRNA001
chr1	test	three_prime_UTR	281	300	.	+	.	ID=utr3_001;Parent=mRNA001
chr1	test	gene	350	520	.	-	.	ID=gene002;Name=testgene2
chr1	test	mRNA	350	520	.	-	.	ID=mRNA002;Parent=gene002;Name=testgene2-RA
chr1	test	exon	350	520	.	-	.	ID=exon003;Parent=mRNA002
chr1	test	CDS	370	500	.	-	0	ID=cds003;Parent=mRNA002
chr1	test	three_prime_UTR	350	369	.	-	.	ID=utr3_002;Parent=mRNA002
chr1	test	five_prime_UTR	501	520	.	-	.	ID=utr5_002;Parent=mRNA002
chr2	test	gene	50	400	.	+	.	ID=gene003;Name=testgene3
chr2	test	mRNA	50	400	.	+	.	ID=mRNA003;Parent=gene003;Name=testgene3-RA
chr2	test	exon	50	150	.	+	.	ID=exon004;Parent=mRNA003
chr2	test	exon	250	400	.	+	.	ID=exon005;Parent=mRNA003
chr2	test	CDS	70	140	.	+	0	ID=cds004;Parent=mRNA003
chr2	test	CDS	250	380	.	+	1	ID=cds005;Parent=mRNA003
chr2	test	five_prime_UTR	50	69	.	+	.	ID=utr5_003;Parent=mRNA003
chr2	test	three_prime_UTR	381	400	.	+	.	ID=utr3_003;Parent=mRNA003
```

- [ ] **Step 2: Write failing tests**

`tests/test_io/test_gff.py`:

```python
"""Tests for GFF3 I/O module."""

import pytest

from fossil_finder.io.gff import (
    parse_gff3,
    iter_gff3,
    get_features_by_type,
    get_children,
    get_gene_to_transcripts,
)


class TestParseGff3:
    def test_loads_features(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        assert len(features) > 0

    def test_feature_has_required_fields(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        f = features[0]
        assert "seqid" in f
        assert "type" in f
        assert "start" in f
        assert "end" in f
        assert "strand" in f
        assert "attributes" in f

    def test_coordinates_are_integers(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        for f in features:
            assert isinstance(f["start"], int)
            assert isinstance(f["end"], int)

    def test_skips_comment_lines(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        for f in features:
            assert not f["seqid"].startswith("#")

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_gff3(tmp_path / "nope.gff3")


class TestGetFeaturesByType:
    def test_get_genes(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        genes = get_features_by_type(features, "gene")
        assert len(genes) == 3

    def test_get_utrs(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        utr3s = get_features_by_type(features, "three_prime_UTR")
        assert len(utr3s) == 3

    def test_get_nonexistent_type(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        enhancers = get_features_by_type(features, "enhancer")
        assert enhancers == []


class TestGetChildren:
    def test_get_exons_of_mrna(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        children = get_children(features, "mRNA001")
        child_types = {c["type"] for c in children}
        assert "exon" in child_types
        assert "CDS" in child_types

    def test_get_mrnas_of_gene(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        children = get_children(features, "gene001")
        assert len(children) == 1
        assert children[0]["type"] == "mRNA"


class TestGetGeneToTranscripts:
    def test_maps_genes_to_transcripts(self, mini_annotation_gff):
        features = parse_gff3(mini_annotation_gff)
        mapping = get_gene_to_transcripts(features)
        assert "gene001" in mapping
        assert "mRNA001" in mapping["gene001"]
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
pytest tests/test_io/test_gff.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.io.gff'`

- [ ] **Step 4: Implement GFF3 parser**

`src/fossil_finder/io/gff.py`:

```python
"""Standards-compliant GFF3 parser for fossil_finder.

Parses any valid GFF3 file (FlyBase, Ensembl, NCBI, custom).
No organism-specific assumptions — feature types and attributes
are treated generically.

GFF3 spec: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
"""

from pathlib import Path
from urllib.parse import unquote


def _parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 attribute column (col 9) into key-value dict."""
    attrs = {}
    if attr_string == ".":
        return attrs
    for pair in attr_string.split(";"):
        pair = pair.strip()
        if "=" not in pair:
            continue
        key, value = pair.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def parse_gff3(path: str | Path) -> list[dict]:
    """Parse GFF3 file into list of feature dicts.

    Each feature dict has: seqid, source, type, start (int), end (int),
    score, strand, phase, attributes (dict).

    Coordinates are 1-based inclusive (GFF3 standard).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")

    features = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            feature = {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": parts[6],
                "phase": None if parts[7] == "." else int(parts[7]),
                "attributes": _parse_attributes(parts[8]),
            }
            features.append(feature)

    return features


def iter_gff3(path: str | Path):
    """Iterate over GFF3 features one at a time (memory-efficient)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GFF3 file not found: {path}")

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            yield {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": parts[6],
                "phase": None if parts[7] == "." else int(parts[7]),
                "attributes": _parse_attributes(parts[8]),
            }


def get_features_by_type(features: list[dict], feature_type: str) -> list[dict]:
    """Filter features to those matching a specific type (e.g. 'gene', 'exon')."""
    return [f for f in features if f["type"] == feature_type]


def get_children(features: list[dict], parent_id: str) -> list[dict]:
    """Get all features whose Parent attribute matches the given ID."""
    return [
        f for f in features
        if parent_id in f["attributes"].get("Parent", "").split(",")
    ]


def get_gene_to_transcripts(features: list[dict]) -> dict[str, list[str]]:
    """Build mapping from gene IDs to their transcript (mRNA) IDs.

    Looks for mRNA features whose Parent is a gene.
    """
    genes = {
        f["attributes"]["ID"]: []
        for f in features
        if f["type"] == "gene" and "ID" in f["attributes"]
    }

    for f in features:
        if f["type"] == "mRNA" and "Parent" in f["attributes"]:
            parent = f["attributes"]["Parent"]
            if parent in genes and "ID" in f["attributes"]:
                genes[parent].append(f["attributes"]["ID"])

    return genes
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
pytest tests/test_io/test_gff.py -v
```

Expected: All 11 tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/io/gff.py tests/test_io/test_gff.py tests/data/mini_annotation.gff3
git commit -m "feat: standards-compliant GFF3 parser (organism-agnostic)

Parses any valid GFF3 (FlyBase, Ensembl, NCBI, custom).
No hardcoded feature types or attribute assumptions.
Supports parent-child traversal for gene→mRNA→exon hierarchies."
```

---

## Chunk 4: Synthetic Test Genome Config

### Task 8: Create test genome config and verify full integration

**Files:**
- Create: `tests/data/mini_genome_config.yaml`
- Test: `tests/test_config/test_schema.py` (extend with integration test)

**Context:** Wire everything together — config loads, points to test fixtures, and the I/O modules can read those fixtures. This validates the full Phase 1 stack.

- [ ] **Step 1: Create test genome config**

`tests/data/mini_genome_config.yaml`:

```yaml
genome:
  species: "Testus synthetica"
  assembly: "test_v1"
  chromosomes: ["chr1", "chr2"]

source:
  adapter: "custom"
  genome_fasta: "mini_genome.fasta"
  annotation_gff: "mini_annotation.gff3"
  te_consensus: "mini_tes.fasta"

blast:
  word_size: 7
  dust: false

gene_sets: {}
```

- [ ] **Step 2: Write integration test**

Append to `tests/test_config/test_schema.py`:

```python
class TestIntegration:
    """Verify config + I/O modules work together."""

    def test_config_points_to_readable_fixtures(self, mini_genome_config, test_data_dir):
        config = load_genome_config(mini_genome_config)

        from fossil_finder.io.fasta import parse_fasta
        from fossil_finder.io.gff import parse_gff3

        # Resolve paths relative to test data dir
        genome_path = test_data_dir / config.source.genome_fasta
        gff_path = test_data_dir / config.source.annotation_gff
        te_path = test_data_dir / config.source.te_consensus

        seqs = parse_fasta(genome_path)
        assert len(seqs) == len(config.genome.chromosomes)

        features = parse_gff3(gff_path)
        assert len(features) > 0

        tes = parse_fasta(te_path)
        assert len(tes) == 3
```

- [ ] **Step 3: Run full test suite**

```bash
pytest tests/ -v
```

Expected: All tests PASS across config and I/O modules.

- [ ] **Step 4: Commit**

```bash
git add tests/data/mini_genome_config.yaml tests/test_config/test_schema.py
git commit -m "test: integration test validates config + I/O pipeline end-to-end

Synthetic test genome config points to FASTA and GFF3 fixtures.
Verifies config loading → file parsing works as a complete stack."
```

---

### Task 9: Update package exports and documentation

**Files:**
- Modify: `src/fossil_finder/io/__init__.py`
- Modify: `src/fossil_finder/config/__init__.py`

- [ ] **Step 1: Set up package exports**

`src/fossil_finder/io/__init__.py`:

```python
"""I/O modules for reading and writing bioinformatics file formats."""

from .blast import BLAST_COLUMNS, load_blast_results, classify_strand
from .fasta import parse_fasta, parse_fasta_headers, iter_fasta, write_fasta
from .gff import parse_gff3, get_features_by_type, get_children
```

`src/fossil_finder/config/__init__.py`:

```python
"""Genome configuration management."""

from .schema import GenomeConfig, load_genome_config
```

- [ ] **Step 2: Verify imports work from top-level**

```bash
python -c "
from fossil_finder.config import GenomeConfig, load_genome_config
from fossil_finder.io import BLAST_COLUMNS, parse_fasta, parse_gff3
print('All imports OK')
print(f'BLAST columns: {len(BLAST_COLUMNS)}')
"
```

Expected: `All imports OK` and `BLAST columns: 17`

- [ ] **Step 3: Run full test suite one final time**

```bash
pytest tests/ -v --tb=short
```

Expected: All tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/io/__init__.py src/fossil_finder/config/__init__.py
git commit -m "feat: package exports for config and I/O modules

from fossil_finder.config import GenomeConfig, load_genome_config
from fossil_finder.io import BLAST_COLUMNS, parse_fasta, parse_gff3"
```

---

## Phase 1 Completion Checklist

After all tasks, verify:

- [ ] `pip install -e ".[dev]"` succeeds
- [ ] `python -c "import fossil_finder"` works
- [ ] `pytest tests/ -v` — all tests pass
- [ ] Old scripts still work: `cd scripts && python -c "from utils.blast_io import BLAST_COLUMNS"`
- [ ] Dmel config validates: `python -c "from fossil_finder.config import load_genome_config; load_genome_config('src/fossil_finder/config/genomes/dmel_r6.66.yaml')"`

---

## Phases 2-9 Overview (to be planned in detail separately)

### Phase 2: Adapter Layer
- Define `GenomeAdapter` ABC in `src/fossil_finder/adapters/base.py`
- Port FlyBase-specific code from `data_loaders.py` + `annotation_loaders.py` → `adapters/flybase.py`
- Stub `adapters/ensembl.py` and `adapters/custom.py`
- Test: Each adapter loads its format correctly

### Phase 3: Region Extraction
- Generic extractor: (genome FASTA + GFF + feature type) → (sequences FASTA + metadata TSV)
- Port `extract_promoters.py`, `extract_exons.py`, etc. as thin wrappers
- Chromosomes come from config, not hardcoded sets
- Port `shuffle_sequences.py` (already generic)

### Phase 4: BLAST Search + Deduplication
- Port `blast_runner.py` — DB paths from config
- Merge `deduplicate_te_hits.py` + `deduplicate_exon_te_hits.py` into one module
- Port `te_parameter_sweep.py`

### Phase 5: Core Analysis
- Port analysis modules: enrichment, strand bias, density, TE families
- All read config for species-specific parameters
- Conservation/synteny: optional, skip if not configured

### Phase 6: Snakemake Workflow
- Write rules for each pipeline stage
- `snakemake --configfile genomes/dmel_r6.66.yaml` runs full pipeline
- Each rule maps to one `fossil_finder` function

### Phase 7: Functional Integration + Reporting
- Expression, GO, localization — adapter-driven
- HTML report generation
- Figure generation

### Phase 8: NTv3 Scoring
- Port scoring pipeline (already has species detection)
- Integrate with config for model/species selection

### Phase 9: Legacy Validation + Archive
- Run new pipeline against Dmel, compare outputs to existing results
- Archive old scripts to `scripts/legacy/`
- Remove v1 scripts (keep only v2 logic)
- Update CLAUDE.md and FILE_MAP.md
