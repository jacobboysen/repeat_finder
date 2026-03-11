# Fossil Finder Phase 2: Adapter Layer & TE Classification

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a pluggable adapter layer that encapsulates all organism-specific data loading, plus a configurable TE classification module — making the pipeline work with any genome, not just Dmel.

**Architecture:** An abstract `GenomeAdapter` base class defines the contract for organism-specific operations (gene ID loading, symbol mapping, expression, GO annotations, localization, TE metadata). Concrete adapters (`FlyBaseAdapter`, `CustomAdapter`) implement the interface. A factory function `get_adapter(config)` returns the right adapter based on `config.source.adapter`. The TE domain classifier is ported from `scripts/utils/te_domain_classifier.py` with Dfam-style `#class/family` notation as primary classification and configurable organism-specific family name lists as fallback.

**Tech Stack:** Python 3.11, Pydantic v2, pytest, abstract base classes

**Source files being refactored:**
- `scripts/utils/data_loaders.py` (405 lines) — 67% Dmel-specific
- `scripts/utils/annotation_loaders.py` (650 lines) — 62% Dmel-specific
- `scripts/utils/te_domain_classifier.py` (273 lines) — domain boundaries generic, family lists Dmel-biased

**Phases Overview:**
- **Phase 1** (done): Package skeleton, config schema, generic I/O, test fixtures
- **Phase 2** (this plan): Adapter layer + TE classification
- **Phase 3**: Region extraction (generic from GFF + genome)
- **Phase 4**: BLAST search + deduplication
- **Phase 5**: Core analysis modules
- **Phase 6**: Snakemake workflow DAG
- **Phase 7**: Functional integration + reporting
- **Phase 8**: NTv3 scoring pipeline
- **Phase 9**: Legacy validation + archive old scripts

---

## File Structure (Phase 2)

```
src/fossil_finder/
├── adapters/
│   ├── __init__.py              # Exports: get_adapter, GenomeAdapter
│   ├── base.py                  # GenomeAdapter ABC (8 abstract methods)
│   ├── flybase.py               # FlyBaseAdapter (ports data_loaders + annotation_loaders)
│   └── custom.py                # CustomAdapter (minimal: GFF3 + FASTA, no database)
├── te/
│   ├── __init__.py              # Exports: classify_te_domain, infer_te_class, TEClassifier
│   ├── taxonomy.py              # TE class inference (Dfam notation + configurable families)
│   └── classifier.py            # TE domain position classification (from te_domain_classifier.py)
├── config/
│   └── schema.py                # Modify: add TETaxonomySpec to GenomeConfig
└── ... (existing, unchanged)

tests/
├── data/
│   ├── mini_gene_list.txt           # 5 gene IDs (ID001-ID005)
│   ├── mini_symbol_map.tsv          # 5 rows: gene_id → symbol
│   ├── mini_expression.tsv          # 5 genes × 3 tissues, generic TSV matrix
│   ├── mini_gene_association.gaf    # 8 GO annotations for 3 genes
│   ├── mini_gene_groups.tsv         # 2 groups with 3+2 members
│   ├── mini_flyfish.csv             # 4 genes with localization patterns
│   ├── mini_te_consensus_fb.fasta   # 3 TEs with FlyBase-style headers (name=...)
│   └── mini_genome_config.yaml      # Modify: add functional data paths
├── test_adapters/
│   ├── __init__.py
│   ├── test_base.py                 # ABC contract tests (4 tests)
│   ├── test_flybase.py              # FlyBaseAdapter tests (14 tests)
│   └── test_custom.py               # CustomAdapter tests (6 tests)
├── test_te/
│   ├── __init__.py
│   ├── test_taxonomy.py             # TE class inference tests (10 tests)
│   └── test_classifier.py           # TE domain classification tests (8 tests)
└── conftest.py                      # Modify: add Phase 2 fixtures
```

---

## Chunk 1: GenomeAdapter ABC + Factory

### Task 1: GenomeAdapter Abstract Base Class

**Files:**
- Create: `src/fossil_finder/adapters/__init__.py`
- Create: `src/fossil_finder/adapters/base.py`
- Test: `tests/test_adapters/__init__.py`
- Test: `tests/test_adapters/test_base.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_adapters/test_base.py
"""Tests for GenomeAdapter ABC contract."""

import pytest

from fossil_finder.adapters.base import GenomeAdapter


class TestGenomeAdapterABC:
    def test_cannot_instantiate_abc_directly(self):
        """ABC should not be instantiable without implementing all methods."""
        with pytest.raises(TypeError, match="abstract"):
            GenomeAdapter.__init__  # just verify it's defined
            # Actually try to instantiate
            GenomeAdapter(config=None)

    def test_concrete_subclass_must_implement_all_methods(self):
        """Partial implementation should also fail."""
        class PartialAdapter(GenomeAdapter):
            def load_gene_ids(self, path):
                return []

        with pytest.raises(TypeError, match="abstract"):
            PartialAdapter(config=None)

    def test_complete_subclass_can_instantiate(self):
        """A complete implementation should work."""
        class FullAdapter(GenomeAdapter):
            def load_gene_ids(self, path):
                return []
            def load_gene_id_symbol_map(self):
                return {}
            def load_expression(self):
                return None
            def load_go_annotations(self):
                return None
            def load_gene_groups(self):
                return None
            def load_localization(self):
                return None
            def load_te_metadata(self, path):
                return {}
            def parse_fasta_metadata(self, header):
                return {}

        adapter = FullAdapter(config=None)
        assert adapter is not None
        assert adapter.config is None

    def test_abstract_methods_list_is_complete(self):
        """Verify all 8 abstract methods are declared."""
        import inspect
        abstract_methods = {
            name for name, _ in inspect.getmembers(GenomeAdapter)
            if getattr(getattr(GenomeAdapter, name, None), '__isabstractmethod__', False)
        }
        expected = {
            'load_gene_ids', 'load_gene_id_symbol_map',
            'load_expression', 'load_go_annotations',
            'load_gene_groups', 'load_localization',
            'load_te_metadata', 'parse_fasta_metadata',
        }
        assert abstract_methods == expected
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_base.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.adapters'`

- [ ] **Step 3: Write the GenomeAdapter ABC**

```python
# src/fossil_finder/adapters/base.py
"""Abstract base class for genome-specific data adapters.

Each supported data source (FlyBase, Ensembl, NCBI, custom) provides a
concrete adapter that knows how to parse that source's file formats and
ID conventions. The rest of the pipeline uses only the adapter interface,
never organism-specific parsing directly.
"""

from abc import ABC, abstractmethod
from pathlib import Path

from fossil_finder.config.schema import GenomeConfig


class GenomeAdapter(ABC):
    """Interface for organism-specific data loading.

    Adapters are constructed with a validated GenomeConfig and provide
    uniform access to gene IDs, symbols, expression, GO annotations,
    gene groups, localization, and TE metadata.
    """

    def __init__(self, config: GenomeConfig | None):
        self.config = config

    @abstractmethod
    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load gene identifiers from a gene list file.

        Args:
            path: Path to gene list (one ID per line, or TSV with IDs).

        Returns:
            List of gene identifiers matching this organism's ID format.
        """

    @abstractmethod
    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """Load gene ID → symbol mapping from the configured source.

        Returns:
            Dict mapping gene IDs (e.g. FBgn0000003) to symbols (e.g. Act5C).
            Empty dict if no symbol map is configured.
        """

    @abstractmethod
    def load_expression(self) -> dict[str, dict[str, float]] | None:
        """Load expression data (e.g. RPKM/TPM matrix).

        Returns:
            Dict mapping gene_id → {tissue_name: value}, or None if
            expression data is not configured.
        """

    @abstractmethod
    def load_go_annotations(self) -> dict[str, list[dict]] | None:
        """Load Gene Ontology annotations.

        Returns:
            Dict mapping gene_id → list of GO annotation dicts,
            or None if not configured. Each annotation dict has keys:
            go_id, aspect, evidence, symbol, qualifier.
        """

    @abstractmethod
    def load_gene_groups(self) -> dict[str, list[str]] | None:
        """Load gene groups / pathway memberships.

        Returns:
            Dict mapping group_name → list of gene IDs,
            or None if not configured.
        """

    @abstractmethod
    def load_localization(self) -> dict[str, list[str]] | None:
        """Load subcellular localization data.

        Returns:
            Dict mapping gene_id → list of localization patterns,
            or None if not configured.
        """

    @abstractmethod
    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE database metadata from a FASTA file.

        Returns:
            Dict mapping te_id → {name, te_class, length}.
        """

    @abstractmethod
    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract structured metadata from a FASTA header line.

        Different databases encode metadata differently in headers
        (FlyBase uses key=value; pairs, Ensembl uses space-separated
        key:value, NCBI uses pipe-delimited fields).

        Args:
            header: Full FASTA header line (including '>').

        Returns:
            Dict of extracted key-value metadata.
        """
```

```python
# src/fossil_finder/adapters/__init__.py
"""Genome adapter layer — pluggable organism-specific data loading."""

from .base import GenomeAdapter

__all__ = ["GenomeAdapter", "get_adapter"]


def get_adapter(config) -> GenomeAdapter:
    """Factory: return the appropriate adapter for a genome config.

    Args:
        config: A GenomeConfig instance (config.source.adapter selects the adapter).

    Returns:
        Concrete GenomeAdapter subclass instance.

    Raises:
        ValueError: If the adapter name is not recognized.
    """
    from .custom import CustomAdapter
    from .flybase import FlyBaseAdapter

    _registry: dict[str, type[GenomeAdapter]] = {
        "flybase": FlyBaseAdapter,
        "custom": CustomAdapter,
        # "ensembl": EnsemblAdapter,  # Phase 9+
        # "ncbi": NCBIAdapter,        # Phase 9+
    }

    adapter_name = config.source.adapter
    cls = _registry.get(adapter_name)
    if cls is None:
        raise ValueError(
            f"Unknown adapter '{adapter_name}'. "
            f"Available: {sorted(_registry.keys())}"
        )
    return cls(config)
```

Create empty `__init__.py` files:
```python
# tests/test_adapters/__init__.py
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_base.py -v`
Expected: 4 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/__init__.py src/fossil_finder/adapters/base.py \
       tests/test_adapters/__init__.py tests/test_adapters/test_base.py
git commit -m "feat: add GenomeAdapter ABC and adapter factory"
```

---

### Task 2: Test Fixtures for Phase 2

**Files:**
- Create: `tests/data/mini_gene_list.txt`
- Create: `tests/data/mini_symbol_map.tsv`
- Create: `tests/data/mini_expression.tsv`
- Create: `tests/data/mini_gene_association.gaf`
- Create: `tests/data/mini_gene_groups.tsv`
- Create: `tests/data/mini_flyfish.csv`
- Create: `tests/data/mini_te_consensus_fb.fasta`
- Modify: `tests/conftest.py`

- [ ] **Step 1: Create all test fixture files**

```
# tests/data/mini_gene_list.txt
# Test gene list
FBgn0000001
FBgn0000002
FBgn0000003
FBgn0000004
FBgn0000005
```

```
# tests/data/mini_symbol_map.tsv
FBgn0000001	geneA
FBgn0000002	geneB
FBgn0000003	geneC
FBgn0000004	geneD
FBgn0000005	geneE
```

```
# tests/data/mini_expression.tsv
# Generic expression matrix (gene_id, gene_symbol, then tissue columns)
#gene_primary_id	gene_symbol	gene_fullname	gene_type	ovary	testis	embryo
FBgn0000001	geneA	Gene A protein	protein_coding_gene	125.4	10.2	55.0
FBgn0000002	geneB	Gene B protein	protein_coding_gene	0.5	200.3	15.1
FBgn0000003	geneC	Gene C protein	protein_coding_gene	50.0	50.0	50.0
FBgn0000004	geneD	Gene D protein	protein_coding_gene	0.0	0.0	300.5
FBgn0000005	geneE	Gene E protein	protein_coding_gene	75.2	75.2	75.2
```

```
# tests/data/mini_gene_association.gaf
!gaf-version: 2.2
!generated-by: test
FB	FBgn0000001	geneA		GO:0003674	FB:FBrf0000001	IDA		F	Gene A	geneA	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000001	geneA		GO:0006468	FB:FBrf0000001	IDA		P	Gene A	geneA	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000002	geneB		GO:0005634	FB:FBrf0000002	IEA		C	Gene B	geneB	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000002	geneB	NOT	GO:0005737	FB:FBrf0000002	IEA		C	Gene B	geneB	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000003	geneC		GO:0003700	FB:FBrf0000003	IDA		F	Gene C	geneC	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000003	geneC		GO:0006355	FB:FBrf0000003	IMP		P	Gene C	geneC	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000003	geneC		GO:0005634	FB:FBrf0000003	IDA		C	Gene C	geneC	gene	taxon:7227	20240101	FlyBase
FB	FBgn0000004	geneD		GO:0016020	FB:FBrf0000004	IEA		C	Gene D	geneD	gene	taxon:7227	20240101	FlyBase
```

```
# tests/data/mini_gene_groups.tsv
#FB_group_id	FB_group_symbol	FB_group_name	Parent_FB_group_id	Parent_FB_group_symbol	Group_member_FB_gene_id	Group_member_FB_gene_symbol
FBgg0000001	RIBOSOMAL	Ribosomal proteins			FBgn0000001	geneA
FBgg0000001	RIBOSOMAL	Ribosomal proteins			FBgn0000002	geneB
FBgg0000001	RIBOSOMAL	Ribosomal proteins			FBgn0000003	geneC
FBgg0000002	KINASE	Protein kinases			FBgn0000003	geneC
FBgg0000002	KINASE	Protein kinases			FBgn0000004	geneD
```

```csv
# tests/data/mini_flyfish.csv
Probe,gene,FBgn ID,stage,term
RE00001,geneA,FBcl0000001,stage 1-3,posterior
RE00002,geneA,FBcl0000001,stage 4-5,pole plasm
RE00003,geneB,FBcl0000002,stage 1-3,apical
RE00004,geneC,FBcl0000003,stage 1-3,cytoplasmic
RE00005,geneC,FBcl0000003,stage 4-5,nuclear
RE00006,geneD,FBcl0000004,stage 1-3,NA
```

```
# tests/data/mini_te_consensus_fb.fasta
>FBte0000001 name=gypsy12; class=LTR; family=Gypsy; length=200
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
ATATATATATATATATATATATATATATATATATATATATATATATATATAT
CCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCC
>FBte0000002 name=mariner2; class=DNA; family=Tc1-Mariner; length=150
GGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCC
TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAAT
CCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGG
>FBte0000003 name=jockey3; class=LINE; family=Jockey; length=180
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATATATATGCGCGCGCATATATATGCGCGCGCATATATATGCGCGCGCATAT
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

- [ ] **Step 2: Update conftest.py with Phase 2 fixtures**

Add to `tests/conftest.py`:

```python
@pytest.fixture
def mini_gene_list(test_data_dir) -> Path:
    """Path to synthetic mini gene list (FBgn IDs)."""
    return test_data_dir / "mini_gene_list.txt"


@pytest.fixture
def mini_symbol_map(test_data_dir) -> Path:
    """Path to synthetic gene ID → symbol mapping TSV."""
    return test_data_dir / "mini_symbol_map.tsv"


@pytest.fixture
def mini_expression(test_data_dir) -> Path:
    """Path to synthetic expression matrix TSV."""
    return test_data_dir / "mini_expression.tsv"


@pytest.fixture
def mini_gaf(test_data_dir) -> Path:
    """Path to synthetic GAF file (GO annotations)."""
    return test_data_dir / "mini_gene_association.gaf"


@pytest.fixture
def mini_gene_groups(test_data_dir) -> Path:
    """Path to synthetic gene group data TSV."""
    return test_data_dir / "mini_gene_groups.tsv"


@pytest.fixture
def mini_flyfish(test_data_dir) -> Path:
    """Path to synthetic FlyFISH localization CSV."""
    return test_data_dir / "mini_flyfish.csv"


@pytest.fixture
def mini_te_consensus_fb(test_data_dir) -> Path:
    """Path to synthetic TE consensus FASTA with FlyBase-style headers."""
    return test_data_dir / "mini_te_consensus_fb.fasta"


@pytest.fixture
def flybase_config(test_data_dir) -> "GenomeConfig":
    """A GenomeConfig wired to FlyBase adapter with test fixture paths."""
    from fossil_finder.config.schema import GenomeConfig

    return GenomeConfig(
        genome={
            "species": "Testus synthetica",
            "assembly": "test_v1",
            "chromosomes": ["chr1", "chr2"],
        },
        source={
            "adapter": "flybase",
            "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
            "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            "gene_id_prefix": "FBgn",
            "transcript_id_prefix": "FBtr",
            "te_consensus": str(test_data_dir / "mini_te_consensus_fb.fasta"),
            "gene_symbol_map": str(test_data_dir / "mini_symbol_map.tsv"),
            "functional": {
                "expression": str(test_data_dir / "mini_expression.tsv"),
                "go_annotations": str(test_data_dir / "mini_gene_association.gaf"),
                "gene_groups": str(test_data_dir / "mini_gene_groups.tsv"),
                "localization": str(test_data_dir / "mini_flyfish.csv"),
            },
        },
    )


@pytest.fixture
def custom_config(test_data_dir) -> "GenomeConfig":
    """A GenomeConfig wired to custom adapter (no functional annotations)."""
    from fossil_finder.config.schema import GenomeConfig

    return GenomeConfig(
        genome={
            "species": "Testus synthetica",
            "assembly": "test_v1",
            "chromosomes": ["chr1", "chr2"],
        },
        source={
            "adapter": "custom",
            "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
            "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            "gene_id_prefix": "GENE",
            "te_consensus": str(test_data_dir / "mini_tes.fasta"),
        },
    )
```

- [ ] **Step 3: Commit**

```bash
git add tests/data/mini_gene_list.txt tests/data/mini_symbol_map.tsv \
       tests/data/mini_expression.tsv tests/data/mini_gene_association.gaf \
       tests/data/mini_gene_groups.tsv tests/data/mini_flyfish.csv \
       tests/data/mini_te_consensus_fb.fasta tests/conftest.py
git commit -m "feat: add Phase 2 test fixtures for adapter layer"
```

---

### Task 3: FlyBaseAdapter — Gene ID Loading + Symbol Mapping

**Files:**
- Create: `src/fossil_finder/adapters/flybase.py`
- Test: `tests/test_adapters/test_flybase.py`

**Ported from:** `scripts/utils/data_loaders.py:18-101` (`load_gene_list`, `load_gene_list_with_symbols`), `scripts/utils/annotation_loaders.py:18-126` (`load_fbgn_to_symbol`, `build_symbol_to_fbgn_map`, `build_comprehensive_symbol_map`)

- [ ] **Step 1: Write the failing tests for gene ID loading and symbol mapping**

```python
# tests/test_adapters/test_flybase.py
"""Tests for FlyBaseAdapter."""

import pytest

from fossil_finder.adapters.flybase import FlyBaseAdapter


class TestFlyBaseGeneIDs:
    def test_load_gene_ids_simple_list(self, flybase_config, mini_gene_list):
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(mini_gene_list)
        assert ids == ["FBgn0000001", "FBgn0000002", "FBgn0000003",
                        "FBgn0000004", "FBgn0000005"]

    def test_load_gene_ids_skips_comments_and_blanks(self, flybase_config, tmp_path):
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("# comment\n\nFBgn0000001\nFBgn0000002\n\n")
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["FBgn0000001", "FBgn0000002"]

    def test_load_gene_ids_extracts_from_tsv(self, flybase_config, tmp_path):
        gene_file = tmp_path / "genes.tsv"
        gene_file.write_text("geneA\tFBgn0000001\nFBgn0000002\tgeneB\n")
        adapter = FlyBaseAdapter(flybase_config)
        ids = adapter.load_gene_ids(gene_file)
        assert "FBgn0000001" in ids
        assert "FBgn0000002" in ids


class TestFlyBaseSymbolMap:
    def test_load_symbol_map(self, flybase_config, mini_symbol_map):
        # Point config to our fixture
        adapter = FlyBaseAdapter(flybase_config)
        mapping = adapter.load_gene_id_symbol_map()
        assert mapping["FBgn0000001"] == "geneA"
        assert mapping["FBgn0000005"] == "geneE"
        assert len(mapping) == 5

    def test_symbol_map_returns_empty_when_unconfigured(self, custom_config):
        """CustomAdapter used here to verify behavior with no symbol map."""
        from fossil_finder.adapters.custom import CustomAdapter
        adapter = CustomAdapter(custom_config)
        mapping = adapter.load_gene_id_symbol_map()
        assert mapping == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseGeneIDs -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.adapters.flybase'`

- [ ] **Step 3: Write FlyBaseAdapter with gene ID loading and symbol mapping**

```python
# src/fossil_finder/adapters/flybase.py
"""FlyBase genome adapter.

Handles FlyBase-specific file formats:
- Gene IDs with FBgn prefix
- FlyBase symbol maps (TSV: FBgn → symbol)
- FlyBase RNA-Seq RPKM reports and matrix format
- FlyBase GAF (GO annotations)
- FlyBase gene group data
- FlyFISH localization data
- FlyBase TE FASTA with name=...; headers
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

from fossil_finder.adapters.base import GenomeAdapter
from fossil_finder.config.schema import GenomeConfig


class FlyBaseAdapter(GenomeAdapter):
    """Adapter for FlyBase data sources (Drosophila)."""

    def __init__(self, config: GenomeConfig):
        super().__init__(config)
        self._gene_id_prefix = config.source.gene_id_prefix or "FBgn"

    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load FBgn IDs from a gene list file.

        Handles simple one-ID-per-line files and TSV files where the
        FBgn ID may be in any column.
        """
        path = Path(path)
        prefix = self._gene_id_prefix
        genes = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "\t" in line:
                    for part in line.split("\t"):
                        if part.startswith(prefix):
                            genes.append(part)
                            break
                elif line.startswith(prefix):
                    genes.append(line)

        return genes

    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """Load FBgn → symbol mapping from configured gene_symbol_map."""
        if not self.config.source.gene_symbol_map:
            return {}

        path = Path(self.config.source.gene_symbol_map)
        if not path.exists():
            return {}

        prefix = self._gene_id_prefix
        mapping = {}

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    fbgn = parts[0].strip()
                    symbol = parts[1].strip()
                    if fbgn.startswith(prefix):
                        mapping[fbgn] = symbol

        return mapping

    # Remaining methods stubbed — implemented in subsequent tasks
    def load_expression(self):
        raise NotImplementedError

    def load_go_annotations(self):
        raise NotImplementedError

    def load_gene_groups(self):
        raise NotImplementedError

    def load_localization(self):
        raise NotImplementedError

    def load_te_metadata(self, path):
        raise NotImplementedError

    def parse_fasta_metadata(self, header):
        raise NotImplementedError
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseGeneIDs tests/test_adapters/test_flybase.py::TestFlyBaseSymbolMap -v`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/flybase.py tests/test_adapters/test_flybase.py
git commit -m "feat: FlyBaseAdapter gene ID loading and symbol mapping"
```

---

### Task 4: FlyBaseAdapter — Expression Data

**Files:**
- Modify: `src/fossil_finder/adapters/flybase.py`
- Modify: `tests/test_adapters/test_flybase.py`

**Ported from:** `scripts/utils/annotation_loaders.py:129-272` (`load_rnaseq_expression`, `load_rnaseq_matrix`, `_clean_tissue_name`)

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_adapters/test_flybase.py`:

```python
class TestFlyBaseExpression:
    def test_load_expression_matrix_format(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        expr = adapter.load_expression()
        assert expr is not None
        assert "FBgn0000001" in expr
        # Check tissue values from mini_expression.tsv
        assert abs(expr["FBgn0000001"]["ovary"] - 125.4) < 0.01
        assert abs(expr["FBgn0000002"]["testis"] - 200.3) < 0.01

    def test_load_expression_gene_count(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        expr = adapter.load_expression()
        assert len(expr) == 5

    def test_load_expression_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_expression() is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseExpression -v`
Expected: FAIL — `NotImplementedError`

- [ ] **Step 3: Implement load_expression and _clean_tissue_name**

Replace the `load_expression` stub in `flybase.py` with:

```python
    def load_expression(self) -> dict[str, dict[str, float]] | None:
        """Load FlyBase RNA-Seq expression data (matrix format).

        Supports the gene_rpkm_matrix format: genes as rows, tissues as columns.
        First columns are metadata (gene_primary_id, gene_symbol, etc.),
        followed by tissue/sample RPKM columns.
        """
        func = self.config.source.functional
        if not func or not func.expression:
            return None

        path = Path(func.expression)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        expression = {}
        metadata_cols = {
            "gene_primary_id", "gene_symbol", "gene_fullname", "gene_type",
            "#gene_primary_id",
        }

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            header = None
            fbgn_col = 0
            tissue_names: list[str] = []
            tissue_indices: list[int] = []

            for row in reader:
                if not row or row[0].startswith("##"):
                    continue

                if header is None:
                    header = row
                    for i, col in enumerate(header):
                        col_clean = col.lstrip("#")
                        if col_clean in metadata_cols:
                            if "primary_id" in col_clean or "fbgn" in col_clean.lower():
                                fbgn_col = i
                            continue
                        tissue_indices.append(i)
                        tissue_names.append(self._clean_tissue_name(col))
                    continue

                if len(row) <= fbgn_col:
                    continue

                gene_id = row[fbgn_col]
                if not gene_id.startswith(prefix):
                    continue

                expression[gene_id] = {}
                for t_idx, col_idx in enumerate(tissue_indices):
                    if col_idx < len(row) and t_idx < len(tissue_names):
                        try:
                            expression[gene_id][tissue_names[t_idx]] = float(row[col_idx])
                        except ValueError:
                            expression[gene_id][tissue_names[t_idx]] = 0.0

        return expression

    @staticmethod
    def _clean_tissue_name(name: str) -> str:
        """Clean FlyBase tissue/sample names for use as dictionary keys.

        Handles FlyBase column naming conventions:
        - FlyAtlas2: RNA-Seq_Profile_FlyAtlas2_Adult_Female_Ovary_(FBlc...)
        - modENCODE: mE_mRNA_A_MateF_4d_ovary_(FBlc...)
        - BCM: BCM_1_E2-4hr_(FBlc...)
        """
        # Remove FlyBase library IDs
        name = re.sub(r"\(FBlc\d+\)", "", name).strip("_")

        # Handle FlyAtlas2 format
        if "FlyAtlas2" in name:
            name = re.sub(r"RNA-Seq_Profile_FlyAtlas2_", "", name)
            return re.sub(r"\s+", "_", name.replace("_", " ").lower())

        # Handle modENCODE format
        if name.startswith("mE_mRNA_"):
            name = name.replace("mE_mRNA_", "")
            for old, new in [
                ("_A_", "_adult_"), ("_L3_", "_larva3_"),
                ("_L1_", "_larva1_"), ("_L2_", "_larva2_"),
                ("MateF", "mated_female"), ("MateM", "mated_male"),
                ("VirF", "virgin_female"), ("AdF", "adult_female"),
                ("AdM", "adult_male"),
            ]:
                name = name.replace(old, new)

        # Handle BCM format
        if name.startswith("BCM_"):
            name = name.replace("BCM_1_", "")
            for old, new in [
                ("E2-4hr", "embryo_2_4hr"), ("E14-16hr", "embryo_14_16hr"),
                ("E2-16hr", "embryo_2_16hr"),
            ]:
                name = name.replace(old, new)

        # Remove common prefixes
        for prefix in [r"^BCM_HGSC_", r"^modENCODE_", r"^Knoblich_mRNA_"]:
            name = re.sub(prefix, "", name)

        # Normalize
        name = name.lower()
        name = re.sub(r"[,;:\s]+", "_", name)
        name = re.sub(r"-+", "_", name)
        name = re.sub(r"_+", "_", name)
        return name.strip("_")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseExpression -v`
Expected: 3 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/flybase.py tests/test_adapters/test_flybase.py
git commit -m "feat: FlyBaseAdapter expression data loading with tissue name cleaning"
```

---

### Task 5: FlyBaseAdapter — GO Annotations + Gene Groups

**Files:**
- Modify: `src/fossil_finder/adapters/flybase.py`
- Modify: `tests/test_adapters/test_flybase.py`

**Ported from:** `scripts/utils/annotation_loaders.py:343-478` (`load_go_annotations`, `load_gene_groups`)

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_adapters/test_flybase.py`:

```python
class TestFlyBaseGO:
    def test_load_go_annotations(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        go = adapter.load_go_annotations()
        assert go is not None
        # FBgn0000001 has 2 annotations (F and P)
        assert len(go["FBgn0000001"]) == 2
        aspects = {a["aspect"] for a in go["FBgn0000001"]}
        assert aspects == {"F", "P"}

    def test_go_skips_negative_annotations(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        go = adapter.load_go_annotations()
        # FBgn0000002 has 1 positive + 1 NOT annotation → only 1 kept
        assert len(go["FBgn0000002"]) == 1

    def test_go_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_go_annotations() is None


class TestFlyBaseGeneGroups:
    def test_load_gene_groups(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        groups = adapter.load_gene_groups()
        assert groups is not None
        assert "RIBOSOMAL" in groups
        assert len(groups["RIBOSOMAL"]) == 3
        assert "KINASE" in groups
        assert len(groups["KINASE"]) == 2

    def test_gene_groups_returns_none_when_unconfigured(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "flybase",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        adapter = FlyBaseAdapter(config)
        assert adapter.load_gene_groups() is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseGO tests/test_adapters/test_flybase.py::TestFlyBaseGeneGroups -v`
Expected: FAIL — `NotImplementedError`

- [ ] **Step 3: Implement load_go_annotations and load_gene_groups**

Replace the stubs in `flybase.py`:

```python
    def load_go_annotations(self) -> dict[str, list[dict]] | None:
        """Load GO annotations from FlyBase GAF file.

        GAF 2.2 format. Filters to gene IDs matching the configured
        prefix and skips negative (NOT) annotations.
        """
        func = self.config.source.functional
        if not func or not func.go_annotations:
            return None

        path = Path(func.go_annotations)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        annotations: dict[str, list[dict]] = defaultdict(list)

        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("!"):
                    continue

                parts = line.split("\t")
                if len(parts) < 9:
                    continue

                gene_id = parts[1]
                qualifier = parts[3]

                if "NOT" in qualifier.upper():
                    continue
                if not gene_id.startswith(prefix):
                    continue

                annotations[gene_id].append({
                    "go_id": parts[4],
                    "aspect": parts[8] if len(parts) > 8 else "",
                    "evidence": parts[6] if len(parts) > 6 else "",
                    "symbol": parts[2],
                    "qualifier": qualifier,
                })

        return dict(annotations)

    def load_gene_groups(self) -> dict[str, list[str]] | None:
        """Load FlyBase gene group data.

        Format: FB_group_id, FB_group_symbol, ..., Group_member_FB_gene_id, ...
        """
        func = self.config.source.functional
        if not func or not func.gene_groups:
            return None

        path = Path(func.gene_groups)
        if not path.exists():
            return None

        prefix = self._gene_id_prefix
        groups: dict[str, list[str]] = defaultdict(list)

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#") or "FB_group" in row[0]:
                    continue
                if len(row) >= 7:
                    group_symbol = row[1]
                    member_id = row[5]
                    if member_id.startswith(prefix):
                        groups[group_symbol].append(member_id)

        return dict(groups)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseGO tests/test_adapters/test_flybase.py::TestFlyBaseGeneGroups -v`
Expected: 5 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/flybase.py tests/test_adapters/test_flybase.py
git commit -m "feat: FlyBaseAdapter GO annotations and gene groups loading"
```

---

### Task 6: FlyBaseAdapter — Localization + TE Metadata + FASTA Headers

**Files:**
- Modify: `src/fossil_finder/adapters/flybase.py`
- Modify: `tests/test_adapters/test_flybase.py`

**Ported from:** `scripts/utils/annotation_loaders.py:481-607` (`load_flyfish_localization`), `scripts/utils/data_loaders.py:188-309` (`parse_te_name`, `load_te_database`, `parse_fasta_by_parent`)

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_adapters/test_flybase.py`:

```python
class TestFlyBaseLocalization:
    def test_load_localization(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        loc = adapter.load_localization()
        assert loc is not None
        # geneA has "posterior" and "pole plasm"
        assert "geneA" in loc
        assert "posterior" in loc["geneA"]
        assert "pole plasm" in loc["geneA"]

    def test_localization_skips_na(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        loc = adapter.load_localization()
        # geneD has "NA" — should be skipped
        assert "geneD" not in loc

    def test_localization_returns_none_when_unconfigured(self, custom_config):
        from fossil_finder.adapters.custom import CustomAdapter
        adapter = CustomAdapter(custom_config)
        assert adapter.load_localization() is None


class TestFlyBaseTEMetadata:
    def test_load_te_metadata(self, flybase_config, mini_te_consensus_fb):
        adapter = FlyBaseAdapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)
        assert len(te_meta) == 3
        assert te_meta["FBte0000001"]["name"] == "gypsy12"
        assert te_meta["FBte0000002"]["name"] == "mariner2"
        assert te_meta["FBte0000003"]["name"] == "jockey3"

    def test_te_metadata_has_class_and_length(self, flybase_config, mini_te_consensus_fb):
        adapter = FlyBaseAdapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)
        for te_id, info in te_meta.items():
            assert "name" in info
            assert "te_class" in info
            assert "length" in info
            assert info["length"] > 0


class TestFlyBaseFastaMetadata:
    def test_parse_flybase_header(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        header = ">FBte0000001 name=gypsy12; class=LTR; family=Gypsy; length=200"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["name"] == "gypsy12"
        assert meta["class"] == "LTR"
        assert meta["family"] == "Gypsy"

    def test_parse_header_with_parent(self, flybase_config):
        adapter = FlyBaseAdapter(flybase_config)
        header = ">FBtr0070000 type=three_prime_UTR; loc=2L:100..200; parent=FBgn0031081;"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["type"] == "three_prime_UTR"
        assert meta["parent"] == "FBgn0031081"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_flybase.py::TestFlyBaseLocalization tests/test_adapters/test_flybase.py::TestFlyBaseTEMetadata tests/test_adapters/test_flybase.py::TestFlyBaseFastaMetadata -v`
Expected: FAIL — `NotImplementedError`

- [ ] **Step 3: Implement localization, TE metadata, and FASTA header parsing**

Replace remaining stubs in `flybase.py`:

```python
    def load_localization(self) -> dict[str, list[str]] | None:
        """Load FlyFISH subcellular localization data.

        Uses the symbol column (not FBcl clone IDs) and optionally
        maps to FBgn via the configured symbol map.
        """
        func = self.config.source.functional
        if not func or not func.localization:
            return None

        path = Path(func.localization)
        if not path.exists():
            return None

        localizations: dict[str, set[str]] = defaultdict(set)

        with open(path, encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f)
            header = None
            gene_col = None
            term_col = None

            for row in reader:
                if not row:
                    continue
                # Skip comment lines
                if row[0].startswith("#"):
                    continue

                if header is None:
                    header = [col.lower() for col in row]
                    for i, col in enumerate(header):
                        if col in ("gene", "symbol"):
                            gene_col = i
                        if col in ("term", "pattern", "localization"):
                            term_col = i
                    continue

                if gene_col is None or term_col is None:
                    continue
                if gene_col >= len(row) or term_col >= len(row):
                    continue

                symbol = row[gene_col].strip()
                term = row[term_col].strip()

                if not symbol or not term:
                    continue
                if term.lower() in ("", "na", "n/a", "none"):
                    continue

                localizations[symbol].add(term.lower())

        return {gene: sorted(pats) for gene, pats in localizations.items()}

    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE metadata from FlyBase-format FASTA.

        FlyBase TE headers use: >FBteXXXXXXX name=...; class=...; ...
        """
        path = Path(path)
        from fossil_finder.io.fasta import parse_fasta

        sequences = parse_fasta(path)
        te_info = {}

        with open(path) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                te_id = line.strip().split()[0][1:]
                meta = self.parse_fasta_metadata(line)
                name = meta.get("name", te_id)
                te_info[te_id] = {
                    "name": name,
                    "te_class": meta.get("class", "Unknown"),
                    "length": len(sequences.get(te_id, "")),
                }

        return te_info

    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract key=value metadata from FlyBase FASTA headers.

        Handles headers like:
            >FBte0000001 name=gypsy12; class=LTR; family=Gypsy; length=200
            >FBtr0070000 type=three_prime_UTR; loc=2L:100..200; parent=FBgn0031081;
        """
        meta = {}
        for match in re.finditer(r"(\w+)=([^;]+)", header):
            key = match.group(1)
            value = match.group(2).strip()
            meta[key] = value
        return meta
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_flybase.py -v`
Expected: 14 passed (5 gene ID + 3 expression + 5 GO/groups + 7 loc/TE/header — wait, recount)

Actually total: 3 gene IDs + 2 symbol map + 3 expression + 3 GO + 2 gene groups + 3 localization + 2 TE metadata + 2 FASTA headers = **20 tests** (but only 14 unique test methods spread across 8 test classes, since we're counting by adding up: 3+2+3+3+2+3+2+2 = 20 tests methods, but only 14 added so far because we started with 5 in task 3... let me recount properly)

Run: `pytest tests/test_adapters/test_flybase.py -v`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/flybase.py tests/test_adapters/test_flybase.py
git commit -m "feat: FlyBaseAdapter localization, TE metadata, and FASTA header parsing"
```

---

## Chunk 2: CustomAdapter + TE Classification

### Task 7: CustomAdapter (Minimal Implementation)

**Files:**
- Create: `src/fossil_finder/adapters/custom.py`
- Create: `tests/test_adapters/test_custom.py`

The `CustomAdapter` handles genomes with standard GFF3 + FASTA but no organism-specific database. It uses the config's `gene_id_prefix` for filtering and extracts what it can from standard file formats.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_adapters/test_custom.py
"""Tests for CustomAdapter (generic genome, no database)."""

import pytest

from fossil_finder.adapters.custom import CustomAdapter


class TestCustomGeneIDs:
    def test_load_gene_ids_with_custom_prefix(self, custom_config, tmp_path):
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("GENE001\nGENE002\nOTHER003\n")
        adapter = CustomAdapter(custom_config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["GENE001", "GENE002"]

    def test_load_gene_ids_no_prefix_loads_all(self, test_data_dir, tmp_path):
        """When no gene_id_prefix is configured, load all non-comment lines."""
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("# comment\nGENE001\nGENE002\n")
        adapter = CustomAdapter(config)
        ids = adapter.load_gene_ids(gene_file)
        assert ids == ["GENE001", "GENE002"]


class TestCustomOptionalMethods:
    def test_symbol_map_returns_empty(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_gene_id_symbol_map() == {}

    def test_expression_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_expression() is None

    def test_go_annotations_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_go_annotations() is None

    def test_localization_returns_none(self, custom_config):
        adapter = CustomAdapter(custom_config)
        assert adapter.load_localization() is None


class TestCustomTEMetadata:
    def test_load_te_metadata_from_plain_fasta(self, custom_config, mini_tes_fasta):
        adapter = CustomAdapter(custom_config)
        te_meta = adapter.load_te_metadata(mini_tes_fasta)
        assert len(te_meta) == 3
        # Plain headers: class and family from key=value pairs
        assert te_meta["TE_LTR1"]["te_class"] == "LTR"
        assert te_meta["TE_DNA1"]["te_class"] == "DNA"
        assert te_meta["TE_LINE1"]["te_class"] == "LINE"

    def test_parse_fasta_metadata_key_value(self, custom_config):
        adapter = CustomAdapter(custom_config)
        header = ">TE001 class=LTR; family=gypsy; length=200"
        meta = adapter.parse_fasta_metadata(header)
        assert meta["class"] == "LTR"
        assert meta["family"] == "gypsy"

    def test_parse_fasta_metadata_dfam_notation(self, custom_config):
        """Dfam-style headers: >Gypsy-2_DM#LTR/Gypsy"""
        adapter = CustomAdapter(custom_config)
        header = ">Gypsy-2_DM#LTR/Gypsy"
        meta = adapter.parse_fasta_metadata(header)
        assert meta.get("class") == "LTR"
        assert meta.get("family") == "Gypsy"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_adapters/test_custom.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.adapters.custom'`

- [ ] **Step 3: Write CustomAdapter**

```python
# src/fossil_finder/adapters/custom.py
"""Custom genome adapter for standard GFF3 + FASTA inputs.

Handles genomes with no organism-specific database. Extracts what it
can from standard file formats and FASTA header conventions.
"""

import re
from pathlib import Path

from fossil_finder.adapters.base import GenomeAdapter
from fossil_finder.config.schema import GenomeConfig


class CustomAdapter(GenomeAdapter):
    """Adapter for genomes with standard GFF3 + FASTA (no database)."""

    def __init__(self, config: GenomeConfig):
        super().__init__(config)
        self._gene_id_prefix = config.source.gene_id_prefix or ""

    def load_gene_ids(self, path: str | Path) -> list[str]:
        """Load gene IDs, filtering by configured prefix (if any)."""
        path = Path(path)
        prefix = self._gene_id_prefix
        genes = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "\t" in line:
                    # TSV: take first column
                    gene_id = line.split("\t")[0]
                else:
                    gene_id = line

                if not prefix or gene_id.startswith(prefix):
                    genes.append(gene_id)

        return genes

    def load_gene_id_symbol_map(self) -> dict[str, str]:
        """No symbol map available for custom genomes."""
        return {}

    def load_expression(self) -> None:
        """No expression data for custom genomes."""
        return None

    def load_go_annotations(self) -> None:
        """No GO annotations for custom genomes."""
        return None

    def load_gene_groups(self) -> None:
        """No gene groups for custom genomes."""
        return None

    def load_localization(self) -> None:
        """No localization data for custom genomes."""
        return None

    def load_te_metadata(self, path: str | Path) -> dict[str, dict]:
        """Load TE metadata from FASTA with key=value or Dfam-style headers."""
        path = Path(path)
        from fossil_finder.io.fasta import parse_fasta

        sequences = parse_fasta(path)
        te_info = {}

        with open(path) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                te_id = line.strip().split()[0][1:]
                meta = self.parse_fasta_metadata(line)
                te_info[te_id] = {
                    "name": meta.get("name", te_id),
                    "te_class": meta.get("class", "Unknown"),
                    "length": len(sequences.get(te_id, "")),
                }

        return te_info

    def parse_fasta_metadata(self, header: str) -> dict[str, str]:
        """Extract metadata from FASTA headers.

        Supports two conventions:
        1. Key=value pairs: >ID name=foo; class=LTR; family=gypsy
        2. Dfam notation: >Gypsy-2_DM#LTR/Gypsy
        """
        meta = {}

        # Try key=value pairs first
        for match in re.finditer(r"(\w+)=([^;]+)", header):
            meta[match.group(1)] = match.group(2).strip()

        # Try Dfam notation: ID#class/family
        if not meta:
            first_token = header.strip().split()[0].lstrip(">")
            if "#" in first_token:
                _, class_family = first_token.split("#", 1)
                parts = class_family.split("/", 1)
                meta["class"] = parts[0]
                if len(parts) > 1:
                    meta["family"] = parts[1]

        return meta
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_custom.py -v`
Expected: 9 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/adapters/custom.py tests/test_adapters/__init__.py \
       tests/test_adapters/test_custom.py
git commit -m "feat: CustomAdapter for standard GFF3 + FASTA genomes"
```

---

### Task 8: Factory Function Integration Tests

**Files:**
- Modify: `tests/test_adapters/test_base.py`

- [ ] **Step 1: Write factory integration tests**

Add to `tests/test_adapters/test_base.py`:

```python
from fossil_finder.adapters import get_adapter
from fossil_finder.adapters.flybase import FlyBaseAdapter
from fossil_finder.adapters.custom import CustomAdapter


class TestGetAdapterFactory:
    def test_flybase_adapter_from_config(self, flybase_config):
        adapter = get_adapter(flybase_config)
        assert isinstance(adapter, FlyBaseAdapter)

    def test_custom_adapter_from_config(self, custom_config):
        adapter = get_adapter(custom_config)
        assert isinstance(adapter, CustomAdapter)

    def test_unknown_adapter_raises(self, test_data_dir):
        from fossil_finder.config.schema import GenomeConfig
        # Build config manually with patched adapter to bypass Literal validation
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": str(test_data_dir / "mini_genome.fasta"),
                "annotation_gff": str(test_data_dir / "mini_annotation.gff3"),
            },
        )
        # Monkey-patch the adapter name to test error handling
        config.source.adapter = "nonexistent"
        with pytest.raises(ValueError, match="Unknown adapter"):
            get_adapter(config)
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `pytest tests/test_adapters/test_base.py -v`
Expected: 7 passed (4 ABC tests + 3 factory tests)

- [ ] **Step 3: Commit**

```bash
git add tests/test_adapters/test_base.py
git commit -m "test: factory function integration tests for adapter layer"
```

---

### Task 9: TE Taxonomy Module — Class Inference

**Files:**
- Create: `src/fossil_finder/te/__init__.py`
- Create: `src/fossil_finder/te/taxonomy.py`
- Create: `tests/test_te/__init__.py`
- Create: `tests/test_te/test_taxonomy.py`

**Ported from:** `scripts/utils/data_loaders.py:208-258` (`parse_te_class`) and `scripts/utils/te_domain_classifier.py:71-132` (`infer_te_class`, family name sets)

The taxonomy module consolidates TE class inference. It supports three classification strategies in priority order:
1. Dfam `#class/family` notation in the TE name (most reliable)
2. Configurable family name lists (defaults include common cross-species families)
3. Keyword fallback matching

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_te/test_taxonomy.py
"""Tests for TE taxonomy / class inference."""

import pytest

from fossil_finder.te.taxonomy import infer_te_class, DEFAULT_FAMILY_LISTS


class TestDfamNotation:
    def test_ltr_gypsy_dfam(self):
        assert infer_te_class("Gypsy-2_DM#LTR/Gypsy") == "LTR"

    def test_line_jockey_dfam(self):
        assert infer_te_class("Jockey-1_DM#LINE/Jockey") == "LINE"

    def test_dna_mariner_dfam(self):
        assert infer_te_class("Mariner-N1_DM#DNA/Tc1-Mariner") == "DNA"

    def test_sine_dfam(self):
        assert infer_te_class("SINE1A#SINE/tRNA") == "SINE"

    def test_helitron_dfam(self):
        assert infer_te_class("Helitron-1_DM#RC/Helitron") == "Helitron"


class TestFamilyNameMatching:
    def test_gypsy_without_dfam_tag(self):
        assert infer_te_class("gypsy12") == "LTR"

    def test_jockey_without_dfam_tag(self):
        assert infer_te_class("jockey3") == "LINE"

    def test_mariner_without_dfam_tag(self):
        assert infer_te_class("mariner2") == "DNA"

    def test_helitron_without_dfam_tag(self):
        assert infer_te_class("DINE-1") == "Helitron"

    def test_unknown_te(self):
        assert infer_te_class("completely_novel_element") == "Unknown"


class TestCustomFamilyLists:
    def test_override_with_custom_families(self):
        """Custom family lists should work for non-standard genomes."""
        custom = {
            "LTR": {"metavirus", "chromovirus"},
            "DNA": {"mutator", "cacta"},
        }
        assert infer_te_class("Metavirus-1", family_lists=custom) == "LTR"
        assert infer_te_class("CACTA_Os", family_lists=custom) == "DNA"
        # Default families should NOT match when custom is provided
        assert infer_te_class("gypsy12", family_lists=custom) == "Unknown"


class TestDefaultFamilyLists:
    def test_all_classes_present(self):
        assert "LTR" in DEFAULT_FAMILY_LISTS
        assert "LINE" in DEFAULT_FAMILY_LISTS
        assert "DNA" in DEFAULT_FAMILY_LISTS
        assert "Helitron" in DEFAULT_FAMILY_LISTS

    def test_common_families_included(self):
        """Cross-species families should be in defaults."""
        assert "gypsy" in DEFAULT_FAMILY_LISTS["LTR"]
        assert "copia" in DEFAULT_FAMILY_LISTS["LTR"]
        assert "jockey" in DEFAULT_FAMILY_LISTS["LINE"]
        assert "mariner" in DEFAULT_FAMILY_LISTS["DNA"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_te/test_taxonomy.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.te'`

- [ ] **Step 3: Write the taxonomy module**

```python
# src/fossil_finder/te/taxonomy.py
"""TE class inference from names and identifiers.

Supports three classification strategies in priority order:
1. Dfam notation: ID#class/family (most reliable, used by RepeatMasker/Dfam)
2. Configurable family name matching (defaults cover common cross-species families)
3. Keyword fallback (SINE, LTR in name)

Family lists can be overridden per-genome via config for organisms with
non-standard TE naming conventions.
"""

# Default family lists — covers common families found across many species.
# These are biased toward Drosophila (where the pipeline was developed) but
# include widely-distributed families found in vertebrates and plants too.
DEFAULT_FAMILY_LISTS: dict[str, set[str]] = {
    "LTR": {
        "gypsy", "copia", "bel", "pao", "mdg", "roo", "412", "297",
        "blood", "accord", "tirant", "springer", "opus", "diver",
        "quasimodo", "idefix", "invader", "gtwin", "tabor", "stalker",
        "max", "transpac", "micropia", "blastopia", "17.6", "1731",
    },
    "LINE": {
        "jockey", "doc", "i-element", "f-element", "g-element",
        "x-element", "het-a", "tart", "tahre", "r1", "r2", "cr1",
        "rt1", "baggins", "juan", "ivk", "bs", "fw",
    },
    "DNA": {
        "p-element", "hobo", "pogo", "bari", "s-element",
        "transib", "tc1", "mariner", "piggybac", "mite",
        "protop", "dnarep", "looper", "1360",
    },
    "Helitron": {"helitron", "dine"},
}


def infer_te_class(
    te_name: str,
    family_lists: dict[str, set[str]] | None = None,
) -> str:
    """Infer TE class from a name or identifier.

    Args:
        te_name: TE family name, ID, or Dfam-annotated name.
        family_lists: Optional custom family name lists. If provided,
            replaces DEFAULT_FAMILY_LISTS entirely (use for genomes
            with non-standard TE naming).

    Returns:
        One of: 'LTR', 'LINE', 'DNA', 'SINE', 'Helitron', 'Unknown'.
    """
    name_lower = te_name.lower()

    # Strategy 1: Dfam notation — ID#class/family
    if "#" in name_lower:
        class_part = name_lower.split("#")[1].split("/")[0]
        if "ltr" in class_part:
            return "LTR"
        if "line" in class_part:
            return "LINE"
        if "dna" in class_part:
            return "DNA"
        if "sine" in class_part:
            return "SINE"
        if "rc" in class_part or "helitron" in class_part:
            return "Helitron"

    # Strategy 2: Family name matching
    families = family_lists if family_lists is not None else DEFAULT_FAMILY_LISTS
    for te_class, names in families.items():
        for fam in names:
            if fam in name_lower:
                return te_class

    # Strategy 3: Keyword fallback (only when using defaults)
    if family_lists is None:
        if "sine" in name_lower or name_lower.startswith("ine-"):
            return "SINE"
        if "ltr" in name_lower:
            return "LTR"

    return "Unknown"
```

```python
# src/fossil_finder/te/__init__.py
"""TE (Transposable Element) classification and domain analysis."""

from .taxonomy import infer_te_class, DEFAULT_FAMILY_LISTS

__all__ = ["infer_te_class", "DEFAULT_FAMILY_LISTS"]
```

Create empty init:
```python
# tests/test_te/__init__.py
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_te/test_taxonomy.py -v`
Expected: 12 passed

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/te/__init__.py src/fossil_finder/te/taxonomy.py \
       tests/test_te/__init__.py tests/test_te/test_taxonomy.py
git commit -m "feat: TE taxonomy module with Dfam notation and configurable family lists"
```

---

### Task 10: TE Domain Classifier

**Files:**
- Create: `src/fossil_finder/te/classifier.py`
- Modify: `src/fossil_finder/te/__init__.py`
- Create: `tests/test_te/test_classifier.py`

**Ported from:** `scripts/utils/te_domain_classifier.py` (273 lines) — domain boundary maps, `classify_te_domain`, `get_relative_position`, helper functions.

The domain boundary maps (LTR: 5'LTR-gag-pol-3'LTR, LINE: 5'UTR-ORF1-ORF2, DNA: TIR-transposase-TIR, Helitron: 5'end-helicase-internal-hairpin) are based on general TE biology and are NOT organism-specific. They port directly.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_te/test_classifier.py
"""Tests for TE domain position classifier."""

import pytest

from fossil_finder.te.classifier import (
    classify_te_domain,
    get_relative_position,
    get_domain_description,
    is_coding_domain,
    is_regulatory_domain,
)


class TestRelativePosition:
    def test_start_of_element(self):
        start, end = get_relative_position(1, 50, 1000)
        assert abs(start - 0.001) < 0.01
        assert abs(end - 0.05) < 0.01

    def test_middle_of_element(self):
        start, end = get_relative_position(400, 600, 1000)
        assert abs(start - 0.4) < 0.01
        assert abs(end - 0.6) < 0.01

    def test_reversed_orientation(self):
        """sstart > send means minus strand; positions should still be ordered."""
        start, end = get_relative_position(600, 400, 1000)
        assert start < end


class TestClassifyTEDomain:
    def test_ltr_gag_region(self):
        """Hit at 20-30% of an LTR element should be gag."""
        result = classify_te_domain("gypsy1", 200, 300, 1000, te_class="LTR")
        assert result["te_class"] == "LTR"
        assert result["domain"] == "gag"

    def test_ltr_pol_with_subdomain(self):
        """Hit at 50% of an LTR element should be pol/rt."""
        result = classify_te_domain("copia1", 450, 550, 1000, te_class="LTR")
        assert result["domain"] == "pol"
        assert result["domain_detail"] == "rt"

    def test_line_orf2_rt(self):
        """Hit at 50% of a LINE should be orf2_rt."""
        result = classify_te_domain("jockey1", 450, 550, 1000, te_class="LINE")
        assert result["domain"] == "orf2_rt"

    def test_dna_transposase(self):
        """Hit at 50% of a DNA transposon should be transposase."""
        result = classify_te_domain("mariner1", 400, 600, 1000, te_class="DNA")
        assert result["domain"] == "transposase"

    def test_infers_class_from_name(self):
        """When te_class is not provided, infer from name."""
        result = classify_te_domain("gypsy12", 200, 300, 1000)
        assert result["te_class"] == "LTR"


class TestDomainHelpers:
    def test_coding_domain_gag(self):
        assert is_coding_domain("gag") is True

    def test_regulatory_domain_ltr(self):
        assert is_regulatory_domain("5_ltr") is True

    def test_domain_description_exists(self):
        desc = get_domain_description("rt")
        assert "Reverse transcriptase" in desc
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_te/test_classifier.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'fossil_finder.te.classifier'`

- [ ] **Step 3: Write the classifier module**

```python
# src/fossil_finder/te/classifier.py
"""TE domain position classifier.

Uses position-based heuristics to classify which region of a TE element
(gag, pol, env, LTR, transposase, etc.) a BLAST hit maps to.

Domain boundaries are approximate and based on typical TE structures:
- LTR retrotransposons: 5'LTR - gag - pol - (env) - 3'LTR
- LINE elements: 5'UTR - ORF1 - ORF2 (RT + endonuclease)
- DNA transposons: TIR - transposase - TIR
- Helitrons: 5'end - Rep/helicase - internal - 3'hairpin

Note: These are heuristic approximations. For precise domain mapping,
use Dfam HMMs or curated per-family annotations.
"""

from fossil_finder.te.taxonomy import infer_te_class

# Domain boundaries as fractions of total element length.
# Based on typical TE structures (not organism-specific).

LTR_DOMAINS = {
    "5_ltr": (0.0, 0.12),
    "gag": (0.12, 0.35),
    "pol": (0.35, 0.88),
    "3_ltr": (0.88, 1.0),
}

POL_SUBDOMAINS = {
    "protease": (0.35, 0.42),
    "rt": (0.42, 0.62),
    "rnaseH": (0.62, 0.72),
    "integrase": (0.72, 0.88),
}

LINE_DOMAINS = {
    "5_utr": (0.0, 0.08),
    "orf1": (0.08, 0.35),
    "orf2_rt": (0.35, 0.70),
    "orf2_en": (0.70, 0.92),
    "3_utr": (0.92, 1.0),
}

DNA_DOMAINS = {
    "5_tir": (0.0, 0.10),
    "transposase": (0.10, 0.90),
    "3_tir": (0.90, 1.0),
}

HELITRON_DOMAINS = {
    "5_end": (0.0, 0.05),
    "rep_helicase": (0.05, 0.60),
    "internal": (0.60, 0.95),
    "3_hairpin": (0.95, 1.0),
}

CLASS_TO_DOMAINS = {
    "LTR": LTR_DOMAINS,
    "LINE": LINE_DOMAINS,
    "DNA": DNA_DOMAINS,
    "Helitron": HELITRON_DOMAINS,
}

_CODING_DOMAINS = {
    "gag", "pol", "protease", "rt", "rnaseH", "integrase",
    "orf1", "orf2_rt", "orf2_en", "transposase", "rep_helicase",
}

_REGULATORY_DOMAINS = {
    "5_ltr", "3_ltr", "5_utr", "3_utr", "5_tir", "3_tir",
    "5_end", "3_hairpin", "internal",
}

_DOMAIN_DESCRIPTIONS = {
    "5_ltr": "5' Long Terminal Repeat (regulatory)",
    "3_ltr": "3' Long Terminal Repeat (regulatory)",
    "gag": "Gag capsid protein (structural)",
    "pol": "Pol polyprotein (RT, integrase)",
    "protease": "Protease domain",
    "rt": "Reverse transcriptase domain",
    "rnaseH": "RNase H domain",
    "integrase": "Integrase domain",
    "5_utr": "5' UTR (non-coding)",
    "3_utr": "3' UTR (non-coding)",
    "orf1": "ORF1 (RNA binding/chaperone)",
    "orf2_rt": "ORF2 reverse transcriptase",
    "orf2_en": "ORF2 endonuclease",
    "5_tir": "5' Terminal Inverted Repeat",
    "3_tir": "3' Terminal Inverted Repeat",
    "transposase": "Transposase (DNA cutting/joining)",
    "5_end": "5' end motif",
    "rep_helicase": "Rep/helicase domain",
    "internal": "Internal variable region",
    "3_hairpin": "3' hairpin structure",
    "unknown": "Unknown region",
}


def get_relative_position(sstart: int, send: int, te_length: int) -> tuple[float, float]:
    """Calculate relative position within TE as fraction (0-1).

    Handles both orientations (sstart < send and sstart > send).
    """
    pos1 = min(sstart, send) / te_length
    pos2 = max(sstart, send) / te_length
    return (
        max(0.0, min(1.0, pos1)),
        max(0.0, min(1.0, pos2)),
    )


def classify_te_domain(
    te_id: str,
    sstart: int,
    send: int,
    te_length: int,
    te_class: str | None = None,
) -> dict:
    """Classify which TE domain a BLAST hit maps to.

    Args:
        te_id: TE identifier (used to infer class if not provided).
        sstart: Subject alignment start position.
        send: Subject alignment end position.
        te_length: Total length of TE element.
        te_class: Optional pre-computed TE class.

    Returns:
        Dict with keys: te_class, domain, domain_detail, position,
        rel_start, rel_end.
    """
    if te_class is None:
        te_class = infer_te_class(te_id)

    rel_start, rel_end = get_relative_position(sstart, send, te_length)
    midpoint = (rel_start + rel_end) / 2

    domains = CLASS_TO_DOMAINS.get(te_class, {})
    primary_domain = "unknown"
    domain_detail = None

    for domain_name, (dom_start, dom_end) in domains.items():
        if dom_start <= midpoint <= dom_end:
            primary_domain = domain_name
            break

    if te_class == "LTR" and primary_domain == "pol":
        for sub_name, (sub_start, sub_end) in POL_SUBDOMAINS.items():
            if sub_start <= midpoint <= sub_end:
                domain_detail = sub_name
                break

    if midpoint < 0.2:
        position = "start"
    elif midpoint > 0.8:
        position = "end"
    else:
        position = "middle"

    return {
        "te_class": te_class,
        "domain": primary_domain,
        "domain_detail": domain_detail,
        "position": position,
        "rel_start": round(rel_start, 3),
        "rel_end": round(rel_end, 3),
    }


def get_domain_description(domain: str) -> str:
    """Get human-readable description of a TE domain."""
    return _DOMAIN_DESCRIPTIONS.get(domain, domain)


def is_coding_domain(domain: str) -> bool:
    """Check if domain is a coding/ORF region."""
    return domain in _CODING_DOMAINS


def is_regulatory_domain(domain: str) -> bool:
    """Check if domain is a regulatory/non-coding region."""
    return domain in _REGULATORY_DOMAINS
```

Update `src/fossil_finder/te/__init__.py`:

```python
"""TE (Transposable Element) classification and domain analysis."""

from .classifier import classify_te_domain, get_domain_description, is_coding_domain, is_regulatory_domain
from .taxonomy import DEFAULT_FAMILY_LISTS, infer_te_class

__all__ = [
    "classify_te_domain",
    "get_domain_description",
    "infer_te_class",
    "is_coding_domain",
    "is_regulatory_domain",
    "DEFAULT_FAMILY_LISTS",
]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_te/test_classifier.py -v`
Expected: 10 passed

- [ ] **Step 5: Run full test suite**

Run: `pytest tests/ -v`
Expected: All tests pass (Phase 1 tests: 49 + Phase 2 tests: ~45 ≈ 94 total)

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/te/classifier.py src/fossil_finder/te/__init__.py \
       tests/test_te/test_classifier.py
git commit -m "feat: TE domain classifier ported from legacy te_domain_classifier.py"
```

---

## Chunk 3: Integration + Final Validation

### Task 11: End-to-End Integration Tests

**Files:**
- Modify: `tests/test_adapters/test_flybase.py` (add integration class)
- Modify: `tests/test_adapters/test_base.py` (add cross-cutting tests)

- [ ] **Step 1: Write integration tests**

Add to `tests/test_adapters/test_flybase.py`:

```python
class TestFlyBaseIntegration:
    """End-to-end: config → adapter → all data loading methods."""

    def test_full_pipeline_load(self, flybase_config):
        """Load all data sources through the adapter interface."""
        from fossil_finder.adapters import get_adapter

        adapter = get_adapter(flybase_config)

        # Symbol map
        symbols = adapter.load_gene_id_symbol_map()
        assert len(symbols) > 0

        # Expression
        expr = adapter.load_expression()
        assert expr is not None
        assert len(expr) > 0

        # GO
        go = adapter.load_go_annotations()
        assert go is not None
        assert len(go) > 0

        # Gene groups
        groups = adapter.load_gene_groups()
        assert groups is not None
        assert len(groups) > 0

        # Localization
        loc = adapter.load_localization()
        assert loc is not None
        assert len(loc) > 0

    def test_te_metadata_feeds_classifier(self, flybase_config, mini_te_consensus_fb):
        """TE metadata from adapter feeds directly into classifier."""
        from fossil_finder.adapters import get_adapter
        from fossil_finder.te.classifier import classify_te_domain

        adapter = get_adapter(flybase_config)
        te_meta = adapter.load_te_metadata(mini_te_consensus_fb)

        for te_id, info in te_meta.items():
            result = classify_te_domain(
                te_id=te_id,
                sstart=1,
                send=info["length"] // 2,
                te_length=info["length"],
                te_class=info["te_class"],
            )
            assert result["te_class"] in ("LTR", "LINE", "DNA", "Helitron", "Unknown")
            assert result["domain"] != ""
```

Add to `tests/test_adapters/test_base.py`:

```python
class TestAdapterConfigIntegration:
    def test_both_adapters_share_same_interface(self, flybase_config, custom_config):
        """Both adapters return compatible types from all methods."""
        from fossil_finder.adapters import get_adapter

        fb = get_adapter(flybase_config)
        custom = get_adapter(custom_config)

        # Both return dict from symbol map (possibly empty)
        assert isinstance(fb.load_gene_id_symbol_map(), dict)
        assert isinstance(custom.load_gene_id_symbol_map(), dict)

        # Both return None or dict from optional methods
        for method_name in ["load_expression", "load_go_annotations",
                            "load_gene_groups", "load_localization"]:
            fb_result = getattr(fb, method_name)()
            custom_result = getattr(custom, method_name)()
            assert fb_result is None or isinstance(fb_result, dict)
            assert custom_result is None or isinstance(custom_result, dict)
```

- [ ] **Step 2: Run the full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: All tests pass

- [ ] **Step 3: Commit**

```bash
git add tests/test_adapters/test_flybase.py tests/test_adapters/test_base.py
git commit -m "test: adapter layer integration tests — config → adapter → TE classifier"
```

---

### Task 12: Update Package Exports + Documentation

**Files:**
- Modify: `src/fossil_finder/__init__.py`
- Modify: `docs/FILE_MAP.md` (if exists and relevant)

- [ ] **Step 1: Update package __init__.py**

```python
# src/fossil_finder/__init__.py
"""fossil_finder — genome-agnostic TE fossil mining pipeline."""

__version__ = "0.2.0"
```

- [ ] **Step 2: Verify everything works with a clean import**

Run:
```bash
python -c "
from fossil_finder.adapters import get_adapter, GenomeAdapter
from fossil_finder.adapters.flybase import FlyBaseAdapter
from fossil_finder.adapters.custom import CustomAdapter
from fossil_finder.te import infer_te_class, classify_te_domain
from fossil_finder.te.taxonomy import DEFAULT_FAMILY_LISTS
from fossil_finder.te.classifier import get_domain_description, is_coding_domain
print('All Phase 2 imports OK')
"
```
Expected: `All Phase 2 imports OK`

- [ ] **Step 3: Run full test suite one final time**

Run: `pytest tests/ -v --tb=short`
Expected: All tests pass

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/__init__.py
git commit -m "chore: bump version to 0.2.0 for Phase 2 adapter layer"
```

---

## Summary

### What Phase 2 delivers:
- **`GenomeAdapter` ABC** with 8 abstract methods defining the organism-specific interface
- **`FlyBaseAdapter`** porting all Dmel-specific parsing from `data_loaders.py` and `annotation_loaders.py`
- **`CustomAdapter`** for any genome with standard GFF3 + FASTA
- **`get_adapter(config)` factory** returning the right adapter based on YAML config
- **TE taxonomy module** with Dfam notation as primary, configurable family lists as fallback
- **TE domain classifier** ported from `te_domain_classifier.py` into the package
- **~45 new tests** across 6 test files

### What remains organism-specific after Phase 2:
- `_clean_tissue_name()` — FlyBase naming conventions (lives inside `FlyBaseAdapter`)
- FlyFISH localization loading — Dmel-only dataset (lives inside `FlyBaseAdapter`)
- FBgn/FBtr prefix handling — configured via `source.gene_id_prefix` in YAML

### What's now generic:
- Gene ID loading (prefix-driven)
- TE class inference (Dfam + configurable families)
- TE domain classification (biology-based, not organism-specific)
- GO annotations from GAF format (standard)
- Gene group loading from TSV (standard format)
- Expression matrix loading (standard TSV matrix)

### Stubs left for future phases:
- `EnsemblAdapter` — Phase 9+ (uses Ensembl BioMart / GTF attributes)
- `NCBIAdapter` — Phase 9+ (uses NCBI Gene / RefSeq)
- `build_comprehensive_symbol_map()` — FlyBase annotation ID file parsing; not needed by the adapter interface (annotation IDs are a FlyBase-specific secondary lookup). Can be added as a FlyBase utility if needed in Phase 7 reporting.
