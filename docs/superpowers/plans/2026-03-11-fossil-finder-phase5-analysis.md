# Fossil Finder Phase 5: Core Analysis Modules

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the analytical core that transforms raw BLAST hits into biological insights — per-gene aggregation, strand bias, TE family distribution, statistical enrichment testing, and RepeatMasker overlap classification.

**Architecture:** Five focused modules under `analysis/`, each with a single responsibility: (1) `aggregation` maps BLAST hits to genes and computes per-gene metrics (hit_count, hit_bp, density); (2) `strand` analyzes sense/antisense bias at gene, TE-family, and genome-wide levels; (3) `families` computes TE family/class distributions and fold-enrichment between gene sets; (4) `enrichment` provides formal statistical tests (Fisher's exact, Mann-Whitney U, Benjamini-Hochberg FDR); (5) `repeatmasker` parses RepeatMasker `.out` files and classifies BLAST hits as known-annotation vs. novel. All modules consume DataFrames from the existing `io/blast.py` and `blast/filter.py` layers.

**Tech Stack:** Python 3.11, pandas, scipy.stats, pytest

---

## File Structure (Phase 5)

```
src/fossil_finder/
├── analysis/
│   ├── __init__.py           # Package exports
│   ├── aggregation.py        # Per-gene hit aggregation (hit_count, hit_bp, density)
│   ├── strand.py             # Strand bias analysis (per-gene, per-TE, genome-wide)
│   ├── families.py           # TE family distribution + fold enrichment
│   ├── enrichment.py         # Fisher's exact, Mann-Whitney U, FDR correction
│   └── repeatmasker.py       # RepeatMasker .out parsing + overlap classification
tests/
├── test_analysis/
│   ├── __init__.py
│   ├── test_aggregation.py   # 11 tests
│   ├── test_strand.py        # 10 tests
│   ├── test_families.py      # 10 tests
│   ├── test_enrichment.py    # 12 tests
│   └── test_repeatmasker.py  # 10 tests
├── data/
│   └── mini_repeatmasker.out # NEW: synthetic RM fixture (8 lines)
```

**Key design decisions:**

1. **All modules operate on DataFrames** — no file I/O inside analysis modules. Parsing happens in `io/` and `analysis/repeatmasker.py` (the one exception, since RM `.out` format is analysis-specific). Analysis functions accept and return DataFrames/dicts.

2. **`aggregation.py` is the bridge** between raw BLAST results and all other analysis modules. It adds `gene_id` and computes density — these are then consumed by strand, families, and enrichment modules.

3. **`enrichment.py` is generic** — it accepts any two groups of values and runs the appropriate test. Not tied to TE analysis specifically.

4. **scipy is required** for Fisher's exact and Mann-Whitney U. It's already in `environment.yml`.

5. **No HTML/report generation** — that belongs in Phase 7 (reporting). This phase produces data structures and DataFrames only.

---

## Chunk 1: Per-Gene Aggregation

### Task 1: Per-gene hit aggregation

The foundation of all downstream analysis. Maps BLAST hits (keyed by query/transcript ID) to genes, and computes per-gene metrics: hit_count, hit_bp, density (hits per kb of query sequence).

**Context for the implementor:** BLAST results have `qseqid` (transcript/region ID) and `qlen` (query length). To aggregate by gene, we need a `query_to_gene` mapping (e.g., transcript → gene). The mapping comes from the adapter layer or from region metadata. This module takes that mapping as an argument — it does not load it from files.

**Files:**
- Create: `src/fossil_finder/analysis/__init__.py`
- Create: `src/fossil_finder/analysis/aggregation.py`
- Create: `tests/test_analysis/__init__.py`
- Create: `tests/test_analysis/test_aggregation.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_analysis/__init__.py` (empty) and `tests/test_analysis/test_aggregation.py`:

```python
"""Tests for per-gene hit aggregation."""

import pandas as pd
import pytest

from fossil_finder.analysis.aggregation import aggregate_by_gene, compute_density


@pytest.fixture
def blast_hits():
    """BLAST results with multiple transcripts mapping to same genes."""
    return pd.DataFrame({
        "qseqid": ["tr1", "tr1", "tr2", "tr2", "tr3"],
        "sseqid": ["TE1", "TE2", "TE1", "TE3", "TE1"],
        "pident": [80.0, 75.0, 90.0, 60.0, 85.0],
        "length": [100, 50, 200, 80, 150],
        "evalue": [1e-10, 1e-5, 1e-20, 0.01, 1e-12],
        "bitscore": [80.0, 40.0, 150.0, 30.0, 100.0],
        "qstart": [10, 200, 50, 300, 100],
        "qend": [110, 250, 250, 380, 250],
        "sstart": [1, 1, 1, 1, 1],
        "send": [100, 50, 200, 80, 150],
        "qlen": [500, 500, 800, 800, 600],
        "slen": [5000, 4000, 5000, 3000, 5000],
        "strand": ["plus", "minus", "plus", "plus", "plus"],
    })


@pytest.fixture
def query_to_gene():
    """Mapping from transcript/query ID to gene ID."""
    return {"tr1": "geneA", "tr2": "geneA", "tr3": "geneB"}


class TestAggregateByGene:
    def test_groups_by_gene(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert len(result) == 2  # geneA, geneB
        assert "geneA" in result.index
        assert "geneB" in result.index

    def test_hit_count(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "hit_count"] == 4  # tr1(2) + tr2(2)
        assert result.loc["geneB", "hit_count"] == 1

    def test_hit_bp(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        # geneA: 100 + 50 + 200 + 80 = 430
        assert result.loc["geneA", "hit_bp"] == 430
        # geneB: 150
        assert result.loc["geneB", "hit_bp"] == 150

    def test_sense_antisense_counts(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "sense_hits"] == 3
        assert result.loc["geneA", "antisense_hits"] == 1
        assert result.loc["geneB", "sense_hits"] == 1
        assert result.loc["geneB", "antisense_hits"] == 0

    def test_query_length_uses_max_per_gene(self, blast_hits, query_to_gene):
        """Gene query length = max qlen across its transcripts."""
        result = aggregate_by_gene(blast_hits, query_to_gene)
        # geneA has tr1(500) and tr2(800) → max = 800
        assert result.loc["geneA", "query_len"] == 800
        assert result.loc["geneB", "query_len"] == 600

    def test_unmapped_queries_skipped(self, blast_hits):
        partial_map = {"tr1": "geneA"}  # tr2, tr3 not mapped
        result = aggregate_by_gene(blast_hits, partial_map)
        assert len(result) == 1
        assert result.loc["geneA", "hit_count"] == 2

    def test_empty_dataframe(self, query_to_gene):
        df = pd.DataFrame(columns=["qseqid", "sseqid", "length", "strand", "qlen"])
        result = aggregate_by_gene(df, query_to_gene)
        assert len(result) == 0

    def test_te_families_tracked(self, blast_hits, query_to_gene):
        result = aggregate_by_gene(blast_hits, query_to_gene)
        assert result.loc["geneA", "n_te_families"] == 3  # TE1, TE2, TE3
        assert result.loc["geneB", "n_te_families"] == 1


class TestComputeDensity:
    def test_density_calculation(self):
        """Density = hit_bp / query_len * 1000 (per kb)."""
        gene_stats = pd.DataFrame({
            "hit_bp": [500, 100],
            "query_len": [1000, 500],
        }, index=["geneA", "geneB"])
        result = compute_density(gene_stats)
        assert result.loc["geneA", "density"] == pytest.approx(500.0)
        assert result.loc["geneB", "density"] == pytest.approx(200.0)

    def test_density_zero_length_safe(self):
        gene_stats = pd.DataFrame({
            "hit_bp": [100],
            "query_len": [0],
        }, index=["geneA"])
        result = compute_density(gene_stats)
        assert result.loc["geneA", "density"] == 0.0

    def test_density_added_as_column(self, blast_hits, query_to_gene):
        agg = aggregate_by_gene(blast_hits, query_to_gene)
        result = compute_density(agg)
        assert "density" in result.columns
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_analysis/test_aggregation.py -v`
Expected: FAIL — `fossil_finder.analysis.aggregation` does not exist

- [ ] **Step 3: Implement aggregation module**

Create `src/fossil_finder/analysis/__init__.py`:

```python
"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density

__all__ = ["aggregate_by_gene", "compute_density"]
```

Create `src/fossil_finder/analysis/aggregation.py`:

```python
"""Per-gene BLAST hit aggregation.

Maps BLAST hits (keyed by query/transcript) to genes and computes
per-gene summary metrics: hit_count, hit_bp, density, strand counts,
and TE family diversity.
"""

import pandas as pd


def aggregate_by_gene(
    df: pd.DataFrame,
    query_to_gene: dict[str, str],
) -> pd.DataFrame:
    """Aggregate BLAST hits by gene.

    Args:
        df: BLAST results DataFrame with at least columns:
            qseqid, length, strand, qlen, sseqid.
        query_to_gene: Mapping from query ID (qseqid) to gene ID.

    Returns:
        DataFrame indexed by gene_id with columns:
        hit_count, hit_bp, sense_hits, antisense_hits,
        query_len, n_te_families.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["hit_count", "hit_bp", "sense_hits", "antisense_hits",
                      "query_len", "n_te_families"],
        )

    # Map queries to genes, drop unmapped
    work = df.copy()
    work["gene_id"] = work["qseqid"].map(query_to_gene)
    work = work.dropna(subset=["gene_id"])

    if work.empty:
        return pd.DataFrame(
            columns=["hit_count", "hit_bp", "sense_hits", "antisense_hits",
                      "query_len", "n_te_families"],
        )

    grouped = work.groupby("gene_id")

    result = pd.DataFrame({
        "hit_count": grouped.size(),
        "hit_bp": grouped["length"].sum(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
        "query_len": grouped["qlen"].max(),
        "n_te_families": grouped["sseqid"].nunique(),
    })

    return result


def compute_density(gene_stats: pd.DataFrame) -> pd.DataFrame:
    """Add TE density column (hit_bp per kb of query sequence).

    Args:
        gene_stats: DataFrame with hit_bp and query_len columns.

    Returns:
        Same DataFrame with added 'density' column.
    """
    result = gene_stats.copy()
    mask = result["query_len"] > 0
    result["density"] = 0.0
    result.loc[mask, "density"] = (
        result.loc[mask, "hit_bp"] / result.loc[mask, "query_len"] * 1000
    )
    return result
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_analysis/test_aggregation.py -v`
Expected: ALL PASS

- [ ] **Step 5: Run full suite**

Run: `pytest tests/ --tb=short`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/analysis/__init__.py src/fossil_finder/analysis/aggregation.py \
  tests/test_analysis/__init__.py tests/test_analysis/test_aggregation.py
git commit -m "feat: per-gene BLAST hit aggregation (hit_count, hit_bp, density)"
```

---

## Chunk 2: Strand Bias + TE Family Distribution

### Task 2: Strand bias analysis

Computes sense/antisense strand bias at three levels: per-gene, per-TE-family, and genome-wide. Classifies bias strength.

**Files:**
- Create: `src/fossil_finder/analysis/strand.py`
- Create: `tests/test_analysis/test_strand.py`
- Modify: `src/fossil_finder/analysis/__init__.py` (add exports)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_analysis/test_strand.py`:

```python
"""Tests for strand bias analysis."""

import pandas as pd
import pytest

from fossil_finder.analysis.strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_te_strand_bias,
    compute_genome_strand_bias,
)


@pytest.fixture
def blast_hits_with_genes():
    """BLAST hits with gene_id column already populated."""
    return pd.DataFrame({
        "qseqid": ["tr1", "tr1", "tr1", "tr2", "tr2"],
        "sseqid": ["TE_gypsy", "TE_gypsy", "TE_copia", "TE_gypsy", "TE_jockey"],
        "gene_id": ["gA", "gA", "gA", "gB", "gB"],
        "strand": ["plus", "plus", "minus", "plus", "minus"],
        "length": [100, 80, 60, 120, 90],
    })


class TestClassifyBias:
    def test_strong_sense(self):
        assert classify_bias(0.85) == "strong_sense"

    def test_sense(self):
        assert classify_bias(0.60) == "sense"

    def test_balanced(self):
        assert classify_bias(0.50) == "balanced"

    def test_antisense(self):
        assert classify_bias(0.35) == "antisense"

    def test_strong_antisense(self):
        assert classify_bias(0.15) == "strong_antisense"


class TestComputeGeneStrandBias:
    def test_per_gene_counts(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "sense_hits"] == 2
        assert result.loc["gA", "antisense_hits"] == 1
        assert result.loc["gB", "sense_hits"] == 1
        assert result.loc["gB", "antisense_hits"] == 1

    def test_sense_pct(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "sense_pct"] == pytest.approx(2 / 3 * 100)
        assert result.loc["gB", "sense_pct"] == pytest.approx(50.0)

    def test_bias_classification(self, blast_hits_with_genes):
        result = compute_gene_strand_bias(blast_hits_with_genes)
        assert result.loc["gA", "bias"] == "sense"
        assert result.loc["gB", "bias"] == "balanced"

    def test_empty_input(self):
        df = pd.DataFrame(columns=["gene_id", "strand", "length"])
        result = compute_gene_strand_bias(df)
        assert len(result) == 0


class TestComputeTEStrandBias:
    def test_per_te_counts(self, blast_hits_with_genes):
        result = compute_te_strand_bias(blast_hits_with_genes)
        assert result.loc["TE_gypsy", "sense_hits"] == 3
        assert result.loc["TE_gypsy", "antisense_hits"] == 0

    def test_min_hits_filter(self, blast_hits_with_genes):
        result = compute_te_strand_bias(blast_hits_with_genes, min_hits=3)
        # Only TE_gypsy has 3+ hits
        assert "TE_gypsy" in result.index
        assert "TE_copia" not in result.index
        assert "TE_jockey" not in result.index


class TestComputeGenomeStrandBias:
    def test_genome_totals(self, blast_hits_with_genes):
        result = compute_genome_strand_bias(blast_hits_with_genes)
        assert result["total_hits"] == 5
        assert result["sense_hits"] == 3
        assert result["antisense_hits"] == 2
        assert result["sense_pct"] == pytest.approx(60.0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_analysis/test_strand.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement strand bias module**

Create `src/fossil_finder/analysis/strand.py`:

```python
"""Strand bias analysis for TE fossil hits.

Computes sense/antisense bias at gene, TE-family, and genome-wide levels.
Strand convention: "plus" = sense (same orientation as query), "minus" = antisense.
"""

import pandas as pd


def classify_bias(sense_fraction: float) -> str:
    """Classify strand bias from sense fraction.

    Args:
        sense_fraction: Fraction of sense hits [0.0, 1.0].

    Returns:
        One of: "strong_sense", "sense", "balanced",
        "antisense", "strong_antisense".
    """
    if sense_fraction >= 0.70:
        return "strong_sense"
    elif sense_fraction >= 0.55:
        return "sense"
    elif sense_fraction >= 0.45:
        return "balanced"
    elif sense_fraction >= 0.30:
        return "antisense"
    else:
        return "strong_antisense"


def compute_gene_strand_bias(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-gene strand bias.

    Args:
        df: BLAST results with gene_id and strand columns.

    Returns:
        DataFrame indexed by gene_id with columns:
        total_hits, sense_hits, antisense_hits, sense_pct, bias.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["total_hits", "sense_hits", "antisense_hits",
                      "sense_pct", "bias"],
        )

    grouped = df.groupby("gene_id")

    result = pd.DataFrame({
        "total_hits": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
    })

    result["sense_pct"] = result.apply(
        lambda r: (r["sense_hits"] / r["total_hits"] * 100)
        if r["total_hits"] > 0 else 0.0,
        axis=1,
    )
    result["bias"] = result["sense_pct"].apply(
        lambda pct: classify_bias(pct / 100)
    )

    return result


def compute_te_strand_bias(
    df: pd.DataFrame, min_hits: int = 1,
) -> pd.DataFrame:
    """Compute per-TE-family strand bias.

    Args:
        df: BLAST results with sseqid and strand columns.
        min_hits: Minimum total hits to include a TE family.

    Returns:
        DataFrame indexed by TE ID with strand bias columns.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["total_hits", "sense_hits", "antisense_hits",
                      "sense_pct", "bias"],
        )

    grouped = df.groupby("sseqid")

    result = pd.DataFrame({
        "total_hits": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
    })

    result = result[result["total_hits"] >= min_hits]

    result["sense_pct"] = result.apply(
        lambda r: (r["sense_hits"] / r["total_hits"] * 100)
        if r["total_hits"] > 0 else 0.0,
        axis=1,
    )
    result["bias"] = result["sense_pct"].apply(
        lambda pct: classify_bias(pct / 100)
    )

    return result


def compute_genome_strand_bias(df: pd.DataFrame) -> dict:
    """Compute genome-wide strand bias summary.

    Args:
        df: BLAST results with strand column.

    Returns:
        Dict with total_hits, sense_hits, antisense_hits, sense_pct, bias.
    """
    if df.empty:
        return {
            "total_hits": 0, "sense_hits": 0, "antisense_hits": 0,
            "sense_pct": 0.0, "bias": "balanced",
        }

    total = len(df)
    sense = (df["strand"] == "plus").sum()
    antisense = (df["strand"] == "minus").sum()
    pct = sense / total * 100 if total > 0 else 0.0

    return {
        "total_hits": total,
        "sense_hits": int(sense),
        "antisense_hits": int(antisense),
        "sense_pct": pct,
        "bias": classify_bias(pct / 100),
    }
```

- [ ] **Step 4: Update `__init__.py` with strand exports**

Update `src/fossil_finder/analysis/__init__.py`:

```python
"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density
from .strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)

__all__ = [
    "aggregate_by_gene",
    "classify_bias",
    "compute_density",
    "compute_gene_strand_bias",
    "compute_genome_strand_bias",
    "compute_te_strand_bias",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_analysis/test_strand.py -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/analysis/strand.py tests/test_analysis/test_strand.py \
  src/fossil_finder/analysis/__init__.py
git commit -m "feat: strand bias analysis (per-gene, per-TE, genome-wide)"
```

---

### Task 3: TE family distribution and fold enrichment

Per-family hit counts, class distribution, and fold-enrichment between two gene sets (e.g., germ plasm vs housekeeping).

**Files:**
- Create: `src/fossil_finder/analysis/families.py`
- Create: `tests/test_analysis/test_families.py`
- Modify: `src/fossil_finder/analysis/__init__.py` (add exports)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_analysis/test_families.py`:

```python
"""Tests for TE family distribution analysis."""

import pandas as pd
import pytest

from fossil_finder.analysis.families import (
    compute_family_stats,
    compute_class_distribution,
    compute_fold_enrichment,
)


@pytest.fixture
def blast_hits():
    return pd.DataFrame({
        "sseqid": ["TE_gypsy", "TE_gypsy", "TE_gypsy", "TE_copia", "TE_jockey"],
        "strand": ["plus", "plus", "minus", "plus", "minus"],
        "pident": [80.0, 75.0, 85.0, 70.0, 90.0],
        "evalue": [1e-10, 1e-8, 1e-12, 1e-5, 1e-15],
        "length": [100, 80, 120, 60, 150],
    })


@pytest.fixture
def te_metadata():
    return {
        "TE_gypsy": {"te_class": "LTR"},
        "TE_copia": {"te_class": "LTR"},
        "TE_jockey": {"te_class": "LINE"},
    }


class TestComputeFamilyStats:
    def test_hit_counts(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "hit_count"] == 3
        assert result.loc["TE_copia", "hit_count"] == 1
        assert result.loc["TE_jockey", "hit_count"] == 1

    def test_strand_counts(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "sense_hits"] == 2
        assert result.loc["TE_gypsy", "antisense_hits"] == 1

    def test_mean_metrics(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "mean_pident"] == pytest.approx(80.0)
        assert result.loc["TE_gypsy", "total_bp"] == 300  # 100+80+120

    def test_frequency(self, blast_hits):
        result = compute_family_stats(blast_hits)
        assert result.loc["TE_gypsy", "frequency"] == pytest.approx(3 / 5)

    def test_empty_input(self):
        df = pd.DataFrame(columns=["sseqid", "strand", "pident", "evalue", "length"])
        result = compute_family_stats(df)
        assert len(result) == 0


class TestComputeClassDistribution:
    def test_class_counts(self, blast_hits, te_metadata):
        result = compute_class_distribution(blast_hits, te_metadata)
        assert result["LTR"] == 4  # gypsy(3) + copia(1)
        assert result["LINE"] == 1

    def test_unknown_class(self, blast_hits):
        """TEs not in metadata get class 'Unknown'."""
        result = compute_class_distribution(blast_hits, {})
        assert result["Unknown"] == 5


class TestComputeFoldEnrichment:
    def test_enrichment_calculation(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1", "TE1", "TE2"],
            "strand": ["plus"] * 4,
            "pident": [80.0] * 4,
            "evalue": [1e-5] * 4,
            "length": [100] * 4,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE1", "TE2", "TE2", "TE2"],
            "strand": ["plus"] * 4,
            "pident": [80.0] * 4,
            "evalue": [1e-5] * 4,
            "length": [100] * 4,
        })
        result = compute_fold_enrichment(set_a, set_b)
        # TE1: freq_a=3/4=0.75, freq_b=1/4=0.25 → enrichment=3.0
        assert result.loc["TE1", "fold_enrichment"] == pytest.approx(3.0)
        # TE2: freq_a=1/4=0.25, freq_b=3/4=0.75 → enrichment=1/3
        assert result.loc["TE2", "fold_enrichment"] == pytest.approx(1 / 3)

    def test_log2_enrichment(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1"], "strand": ["plus"] * 2,
            "pident": [80.0] * 2, "evalue": [1e-5] * 2, "length": [100] * 2,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE1"], "strand": ["plus"],
            "pident": [80.0], "evalue": [1e-5], "length": [100],
        })
        result = compute_fold_enrichment(set_a, set_b)
        # freq_a=1.0, freq_b=1.0 → enrichment=1.0, log2=0.0
        # Both sets have TE1 at 100% frequency → fold=1.0
        assert result.loc["TE1", "log2_enrichment"] == pytest.approx(0.0)

    def test_family_unique_to_one_set(self):
        set_a = pd.DataFrame({
            "sseqid": ["TE1", "TE1"], "strand": ["plus"] * 2,
            "pident": [80.0] * 2, "evalue": [1e-5] * 2, "length": [100] * 2,
        })
        set_b = pd.DataFrame({
            "sseqid": ["TE2"], "strand": ["plus"],
            "pident": [80.0], "evalue": [1e-5], "length": [100],
        })
        result = compute_fold_enrichment(set_a, set_b)
        # TE1 only in set_a → fold_enrichment = inf, log2 = inf
        assert result.loc["TE1", "fold_enrichment"] == float("inf")
        assert result.loc["TE1", "log2_enrichment"] == float("inf")
        # TE2 only in set_b → fold_enrichment = 0.0, log2 = -inf
        assert result.loc["TE2", "fold_enrichment"] == 0.0
        assert result.loc["TE2", "log2_enrichment"] == float("-inf")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_analysis/test_families.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement families module**

Create `src/fossil_finder/analysis/families.py`:

```python
"""TE family distribution and fold-enrichment analysis.

Computes per-TE-family hit statistics and comparative enrichment
between gene sets (e.g., germ plasm vs housekeeping).
"""

import math

import pandas as pd


def compute_family_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-TE-family summary statistics.

    Args:
        df: BLAST results with sseqid, strand, pident, evalue, length columns.

    Returns:
        DataFrame indexed by TE ID with columns:
        hit_count, sense_hits, antisense_hits, mean_pident,
        mean_evalue, total_bp, frequency.
    """
    if df.empty:
        return pd.DataFrame(
            columns=["hit_count", "sense_hits", "antisense_hits",
                      "mean_pident", "mean_evalue", "total_bp", "frequency"],
        )

    grouped = df.groupby("sseqid")
    total = len(df)

    result = pd.DataFrame({
        "hit_count": grouped.size(),
        "sense_hits": grouped["strand"].apply(lambda s: (s == "plus").sum()),
        "antisense_hits": grouped["strand"].apply(lambda s: (s == "minus").sum()),
        "mean_pident": grouped["pident"].mean(),
        "mean_evalue": grouped["evalue"].mean(),
        "total_bp": grouped["length"].sum(),
        "frequency": grouped.size() / total,
    })

    return result


def compute_class_distribution(
    df: pd.DataFrame,
    te_metadata: dict[str, dict],
) -> dict[str, int]:
    """Count BLAST hits by TE class.

    Args:
        df: BLAST results with sseqid column.
        te_metadata: Dict mapping TE ID to metadata dict
            containing at least a 'te_class' key.

    Returns:
        Dict mapping TE class name to hit count.
    """
    counts: dict[str, int] = {}
    for te_id in df["sseqid"]:
        meta = te_metadata.get(te_id, {})
        te_class = meta.get("te_class", "Unknown")
        counts[te_class] = counts.get(te_class, 0) + 1
    return counts


def compute_fold_enrichment(
    set_a: pd.DataFrame,
    set_b: pd.DataFrame,
) -> pd.DataFrame:
    """Compute fold-enrichment of TE families between two hit sets.

    For each TE family, fold_enrichment = freq_in_a / freq_in_b.
    A value > 1 means enriched in set A; < 1 means enriched in set B.

    Args:
        set_a: BLAST results for gene set A.
        set_b: BLAST results for gene set B.

    Returns:
        DataFrame indexed by TE ID with columns:
        count_a, count_b, freq_a, freq_b, fold_enrichment, log2_enrichment.
    """
    stats_a = compute_family_stats(set_a) if not set_a.empty else pd.DataFrame()
    stats_b = compute_family_stats(set_b) if not set_b.empty else pd.DataFrame()

    all_families = set()
    if not stats_a.empty:
        all_families.update(stats_a.index)
    if not stats_b.empty:
        all_families.update(stats_b.index)

    rows = []
    for fam in sorted(all_families):
        count_a = int(stats_a.loc[fam, "hit_count"]) if fam in stats_a.index else 0
        count_b = int(stats_b.loc[fam, "hit_count"]) if fam in stats_b.index else 0
        freq_a = float(stats_a.loc[fam, "frequency"]) if fam in stats_a.index else 0.0
        freq_b = float(stats_b.loc[fam, "frequency"]) if fam in stats_b.index else 0.0

        if freq_b > 0:
            fold = freq_a / freq_b
        elif freq_a > 0:
            fold = float("inf")
        else:
            fold = 0.0

        if fold == float("inf"):
            log2 = float("inf")
        elif fold > 0:
            log2 = math.log2(fold)
        else:
            log2 = float("-inf")

        rows.append({
            "te_id": fam,
            "count_a": count_a,
            "count_b": count_b,
            "freq_a": freq_a,
            "freq_b": freq_b,
            "fold_enrichment": fold,
            "log2_enrichment": log2,
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.set_index("te_id")
    return result
```

- [ ] **Step 4: Update `__init__.py` exports**

```python
"""Core analysis modules for TE fossil mining."""

from .aggregation import aggregate_by_gene, compute_density
from .families import compute_class_distribution, compute_family_stats, compute_fold_enrichment
from .strand import (
    classify_bias,
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)

__all__ = [
    "aggregate_by_gene",
    "classify_bias",
    "compute_class_distribution",
    "compute_density",
    "compute_family_stats",
    "compute_fold_enrichment",
    "compute_gene_strand_bias",
    "compute_genome_strand_bias",
    "compute_te_strand_bias",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_analysis/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/analysis/families.py src/fossil_finder/analysis/strand.py \
  src/fossil_finder/analysis/__init__.py \
  tests/test_analysis/test_strand.py tests/test_analysis/test_families.py
git commit -m "feat: strand bias + TE family distribution analysis"
```

---

## Chunk 3: Statistical Enrichment

### Task 4: Enrichment testing (Fisher's exact, Mann-Whitney, FDR)

The only formal statistical tests in the pipeline. Generic framework: accepts gene sets and metric vectors, returns test results with multiple-testing correction.

**Files:**
- Create: `src/fossil_finder/analysis/enrichment.py`
- Create: `tests/test_analysis/test_enrichment.py`
- Modify: `src/fossil_finder/analysis/__init__.py` (add exports)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_analysis/test_enrichment.py`:

```python
"""Tests for statistical enrichment testing."""

import numpy as np
import pandas as pd
import pytest

from fossil_finder.analysis.enrichment import (
    fisher_exact_enrichment,
    mannwhitney_enrichment,
    correct_multiple_testing,
    test_gene_set_enrichment,
)


class TestFisherExactEnrichment:
    def test_enriched_set(self):
        """Set with more TE-positive genes than background → OR > 1."""
        result = fisher_exact_enrichment(
            set_positive=8, set_negative=2,
            bg_positive=20, bg_negative=80,
        )
        assert result["odds_ratio"] > 1.0
        assert result["p_value"] < 0.05

    def test_depleted_set(self):
        """Set with fewer TE-positive genes → OR < 1."""
        result = fisher_exact_enrichment(
            set_positive=1, set_negative=9,
            bg_positive=50, bg_negative=50,
        )
        assert result["odds_ratio"] < 1.0

    def test_no_difference(self):
        """Equal proportions → OR ≈ 1."""
        result = fisher_exact_enrichment(
            set_positive=5, set_negative=5,
            bg_positive=50, bg_negative=50,
        )
        assert result["odds_ratio"] == pytest.approx(1.0, abs=0.5)
        assert result["p_value"] > 0.05

    def test_returns_required_keys(self):
        result = fisher_exact_enrichment(5, 5, 50, 50)
        assert "odds_ratio" in result
        assert "p_value" in result


class TestMannWhitneyEnrichment:
    def test_higher_group(self):
        """Group with clearly higher values → p < 0.05."""
        group_a = [10.0, 12.0, 15.0, 11.0, 14.0]
        group_b = [1.0, 2.0, 3.0, 1.5, 2.5]
        result = mannwhitney_enrichment(group_a, group_b)
        assert result["p_value"] < 0.05

    def test_equal_groups(self):
        """Identical distributions → p > 0.05."""
        group = [5.0, 5.0, 5.0, 5.0, 5.0]
        result = mannwhitney_enrichment(group, group)
        assert result["p_value"] > 0.05

    def test_returns_required_keys(self):
        result = mannwhitney_enrichment([1.0, 2.0], [3.0, 4.0])
        assert "u_statistic" in result
        assert "p_value" in result

    def test_too_few_samples(self):
        """With < 2 samples, return NaN p-value."""
        result = mannwhitney_enrichment([1.0], [2.0, 3.0])
        assert np.isnan(result["p_value"])


class TestCorrectMultipleTesting:
    def test_bonferroni_correction(self):
        p_values = [0.01, 0.04, 0.06]
        q_values = correct_multiple_testing(p_values, method="bonferroni")
        assert len(q_values) == 3
        # Bonferroni: p * n_tests, capped at 1.0
        assert q_values[0] == pytest.approx(0.03)
        assert q_values[1] == pytest.approx(0.12)
        assert q_values[2] == pytest.approx(0.18)

    def test_fdr_correction(self):
        p_values = [0.01, 0.04, 0.06]
        q_values = correct_multiple_testing(p_values, method="fdr")
        assert len(q_values) == 3
        # BH: rank 1 → 0.01*3/1=0.03, rank 2 → 0.04*3/2=0.06, rank 3 → 0.06*3/3=0.06
        # After monotonicity sweep: [0.03, 0.06, 0.06]
        assert q_values[0] == pytest.approx(0.03)
        assert q_values[1] == pytest.approx(0.06)
        assert q_values[2] == pytest.approx(0.06)

    def test_empty_input(self):
        assert correct_multiple_testing([], method="fdr") == []

    def test_single_pvalue(self):
        q = correct_multiple_testing([0.05], method="bonferroni")
        assert q[0] == pytest.approx(0.05)


class TestGeneSetEnrichment:
    def test_full_enrichment_pipeline(self):
        """Test the combined pipeline with gene set vs background."""
        gene_densities = pd.Series(
            [5.0, 8.0, 12.0, 3.0, 15.0],
            index=["g1", "g2", "g3", "g4", "g5"],
        )
        gene_set = {"g1", "g2", "g3"}  # genes of interest
        te_positive = {"g1", "g2", "g3", "g5"}  # genes with TE hits

        result = test_gene_set_enrichment(
            gene_set=gene_set,
            te_positive_genes=te_positive,
            gene_densities=gene_densities,
        )
        assert "fisher" in result
        assert "mannwhitney" in result
        assert result["fisher"]["odds_ratio"] is not None
        assert result["mannwhitney"]["p_value"] is not None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_analysis/test_enrichment.py -v`
Expected: FAIL — module does not exist

- [ ] **Step 3: Implement enrichment module**

Create `src/fossil_finder/analysis/enrichment.py`:

```python
"""Statistical enrichment testing for TE fossil analysis.

Provides Fisher's exact test (TE presence enrichment), Mann-Whitney U
(density distribution comparison), and multiple testing correction.
"""

import math

import numpy as np
from scipy import stats


def fisher_exact_enrichment(
    set_positive: int,
    set_negative: int,
    bg_positive: int,
    bg_negative: int,
) -> dict:
    """Fisher's exact test for TE presence enrichment.

    Tests whether a gene set has more TE-positive genes than expected
    by the background rate.

    Args:
        set_positive: Genes in set WITH TE hits.
        set_negative: Genes in set WITHOUT TE hits.
        bg_positive: Background genes WITH TE hits (excluding set).
        bg_negative: Background genes WITHOUT TE hits (excluding set).

    Returns:
        Dict with odds_ratio, p_value.
    """
    table = [[set_positive, set_negative], [bg_positive, bg_negative]]
    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")
    return {"odds_ratio": odds_ratio, "p_value": p_value}


def mannwhitney_enrichment(
    group_a: list[float],
    group_b: list[float],
) -> dict:
    """Mann-Whitney U test for density distribution comparison.

    Args:
        group_a: Density values for gene set.
        group_b: Density values for background.

    Returns:
        Dict with u_statistic, p_value.
    """
    if len(group_a) < 2 or len(group_b) < 2:
        return {"u_statistic": float("nan"), "p_value": float("nan")}

    u_stat, p_value = stats.mannwhitneyu(
        group_a, group_b, alternative="two-sided",
    )
    return {"u_statistic": u_stat, "p_value": p_value}


def correct_multiple_testing(
    p_values: list[float],
    method: str = "fdr",
) -> list[float]:
    """Apply multiple testing correction.

    Args:
        p_values: List of raw p-values.
        method: "fdr" (Benjamini-Hochberg) or "bonferroni".

    Returns:
        List of corrected p-values (q-values).
    """
    if not p_values:
        return []

    n = len(p_values)

    if method == "bonferroni":
        return [min(p * n, 1.0) for p in p_values]

    # Benjamini-Hochberg FDR
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [0.0] * n

    for rank, (orig_idx, p) in enumerate(indexed, 1):
        q = min(p * n / rank, 1.0)
        q_values[orig_idx] = q

    # Enforce monotonicity (q-values must be non-decreasing in sorted order)
    sorted_indices = [i for i, _ in indexed]
    for k in range(n - 2, -1, -1):
        idx = sorted_indices[k]
        next_idx = sorted_indices[k + 1]
        if q_values[idx] > q_values[next_idx]:
            q_values[idx] = q_values[next_idx]

    return q_values


def test_gene_set_enrichment(
    gene_set: set[str],
    te_positive_genes: set[str],
    gene_densities: "pd.Series",
) -> dict:
    """Run full enrichment analysis for a gene set.

    Combines Fisher's exact (TE presence) and Mann-Whitney U (density)
    into a single result.

    Args:
        gene_set: Gene IDs in the set of interest.
        te_positive_genes: All gene IDs with at least one TE hit.
        gene_densities: Series mapping gene_id → density value.

    Returns:
        Dict with 'fisher' and 'mannwhitney' sub-dicts.
    """
    all_genes = set(gene_densities.index)

    # Fisher's exact: TE presence
    set_pos = len(gene_set & te_positive_genes)
    set_neg = len(gene_set - te_positive_genes)
    bg_genes = all_genes - gene_set
    bg_pos = len(bg_genes & te_positive_genes)
    bg_neg = len(bg_genes - te_positive_genes)

    fisher_result = fisher_exact_enrichment(set_pos, set_neg, bg_pos, bg_neg)

    # Mann-Whitney U: density distribution
    set_densities = [
        gene_densities[g] for g in gene_set if g in gene_densities.index
    ]
    bg_densities = [
        gene_densities[g] for g in bg_genes if g in gene_densities.index
    ]

    mw_result = mannwhitney_enrichment(set_densities, bg_densities)

    return {"fisher": fisher_result, "mannwhitney": mw_result}
```

- [ ] **Step 4: Update `__init__.py` exports**

Add to `src/fossil_finder/analysis/__init__.py`:

```python
from .enrichment import (
    correct_multiple_testing,
    fisher_exact_enrichment,
    mannwhitney_enrichment,
    test_gene_set_enrichment,
)
```

And add to `__all__`:

```python
    "correct_multiple_testing",
    "fisher_exact_enrichment",
    "mannwhitney_enrichment",
    "test_gene_set_enrichment",
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_analysis/test_enrichment.py -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/analysis/enrichment.py tests/test_analysis/test_enrichment.py \
  src/fossil_finder/analysis/__init__.py
git commit -m "feat: statistical enrichment (Fisher, Mann-Whitney, FDR)"
```

---

## Chunk 4: RepeatMasker Overlap + Version Bump

### Task 5: RepeatMasker parsing and overlap classification

Parses RepeatMasker `.out` files into a DataFrame, detects overlaps with query regions, and classifies BLAST hits as "known" (in RepeatMasker annotation) vs "novel" (not previously annotated).

**Files:**
- Create: `src/fossil_finder/analysis/repeatmasker.py`
- Create: `tests/test_analysis/test_repeatmasker.py`
- Create: `tests/data/mini_repeatmasker.out`
- Modify: `src/fossil_finder/analysis/__init__.py` (add exports)
- Modify: `tests/conftest.py` (add fixture)

- [ ] **Step 1: Create test fixture**

Create `tests/data/mini_repeatmasker.out`:

```
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID
  450  12.3  1.2  0.5  chr1          50       200  (312)  +  TE_gypsy1      LTR/Gypsy              1    150  (4850)   1
  320   8.1  0.0  0.0  chr1         280       350  (162)  C  TE_jockey1     LINE/Jockey          (3850)  150     1     2
  180  15.5  2.0  1.0  chr1         400       450  (62)   +  TE_mariner1    DNA/TcMar-Mariner      1     50  (100)    3
  500   5.2  0.5  0.3  chr2          80       300  (212)  +  TE_gypsy1      LTR/Gypsy              1    220  (4780)   4
  250  18.0  1.5  0.8  chr2         350       420  (92)   C  TE_copia1      LTR/Copia            (4400)  100     1     5
```

- [ ] **Step 2: Write the failing tests**

Add to `tests/conftest.py`:

```python
@pytest.fixture
def mini_repeatmasker(test_data_dir) -> Path:
    """Path to synthetic RepeatMasker .out file."""
    return test_data_dir / "mini_repeatmasker.out"
```

Create `tests/test_analysis/test_repeatmasker.py`:

```python
"""Tests for RepeatMasker parsing and overlap detection."""

import pandas as pd
import pytest

from fossil_finder.analysis.repeatmasker import (
    parse_repeatmasker_out,
    find_overlaps,
    classify_hits,
)


class TestParseRepeatMaskerOut:
    def test_loads_all_records(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        assert len(df) == 5

    def test_columns_present(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        for col in ["chrom", "start", "end", "strand", "repeat_name", "repeat_class"]:
            assert col in df.columns

    def test_coordinate_values(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        first = df.iloc[0]
        assert first["chrom"] == "chr1"
        assert first["start"] == 50
        assert first["end"] == 200
        assert first["repeat_name"] == "TE_gypsy1"
        assert first["repeat_class"] == "LTR/Gypsy"

    def test_strand_parsing(self, mini_repeatmasker):
        df = parse_repeatmasker_out(mini_repeatmasker)
        assert df.iloc[0]["strand"] == "+"
        assert df.iloc[1]["strand"] == "-"  # C → minus

    def test_nonexistent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_repeatmasker_out(tmp_path / "nope.out")


class TestFindOverlaps:
    def test_overlap_detected(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR/Gypsy"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr1"],
            "start": [150], "end": [300],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 1
        assert overlaps.iloc[0]["region_id"] == "utr1"
        assert overlaps.iloc[0]["overlap_bp"] == 51  # 150-200 inclusive

    def test_no_overlap(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr1"],
            "start": [300], "end": [400],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 0

    def test_different_chromosomes_no_overlap(self):
        rm_regions = pd.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "repeat_name": ["TE1"], "repeat_class": ["LTR"],
            "strand": ["+"], "divergence": [10.0],
        })
        query_regions = pd.DataFrame({
            "region_id": ["utr1"], "chrom": ["chr2"],
            "start": [100], "end": [200],
        })
        overlaps = find_overlaps(rm_regions, query_regions)
        assert len(overlaps) == 0


class TestClassifyHits:
    def test_known_vs_novel(self):
        blast_hits = pd.DataFrame({
            "qseqid": ["utr1", "utr1", "utr2"],
            "qstart": [50, 250, 10],
            "qend": [150, 350, 60],
            "sseqid": ["TE1", "TE2", "TE3"],
        })
        rm_overlaps = pd.DataFrame({
            "region_id": ["utr1"],
            "rm_start_in_query": [31],  # 1-based relative to query start
            "rm_end_in_query": [181],   # 1-based relative to query start
            "repeat_name": ["TE_gypsy1"],
            "repeat_class": ["LTR/Gypsy"],
        })
        known, novel = classify_hits(blast_hits, rm_overlaps)
        # Hit 1 (qstart=50, qend=150) overlaps RM region [31,181] → known
        assert len(known) == 1
        assert known.iloc[0]["qseqid"] == "utr1"
        assert known.iloc[0]["qstart"] == 50
        # Hit 2 (qstart=250, qend=350) doesn't overlap → novel
        # Hit 3 (utr2) has no RM data → novel
        assert len(novel) == 2
```

- [ ] **Step 3: Implement RepeatMasker module**

Create `src/fossil_finder/analysis/repeatmasker.py`:

```python
"""RepeatMasker .out file parsing and overlap classification.

Parses RepeatMasker output, detects overlaps with query regions,
and classifies BLAST hits as known (in RM annotation) vs novel.
"""

from pathlib import Path

import pandas as pd


def parse_repeatmasker_out(path: str | Path) -> pd.DataFrame:
    """Parse RepeatMasker .out file into a DataFrame.

    Args:
        path: Path to RepeatMasker .out file.

    Returns:
        DataFrame with columns: score, divergence, chrom, start, end,
        strand, repeat_name, repeat_class.

    Raises:
        FileNotFoundError: If file does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"RepeatMasker file not found: {path}")

    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("SW") or line.startswith("score"):
                continue

            parts = line.split()
            if len(parts) < 15:
                continue

            try:
                score = int(parts[0])
                divergence = float(parts[1])
                chrom = parts[4]
                start = int(parts[5])
                end = int(parts[6])
                strand_raw = parts[8]
                strand = "-" if strand_raw == "C" else "+"
                repeat_name = parts[9]
                repeat_class = parts[10]
            except (ValueError, IndexError):
                continue

            records.append({
                "score": score,
                "divergence": divergence,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "repeat_name": repeat_name,
                "repeat_class": repeat_class,
            })

    return pd.DataFrame(records)


def find_overlaps(
    rm_regions: pd.DataFrame,
    query_regions: pd.DataFrame,
) -> pd.DataFrame:
    """Find overlaps between RepeatMasker regions and query regions.

    Uses simple interval intersection. Both inputs use 1-based inclusive
    coordinates (GFF3/RM convention).

    Args:
        rm_regions: RepeatMasker regions with chrom, start, end,
            repeat_name, repeat_class, strand, divergence.
        query_regions: Query regions with region_id, chrom, start, end.

    Returns:
        DataFrame of overlaps with region_id, overlap_bp, and RM metadata.
    """
    overlaps = []

    for _, rm in rm_regions.iterrows():
        for _, qr in query_regions.iterrows():
            if rm["chrom"] != qr["chrom"]:
                continue

            overlap_start = max(rm["start"], qr["start"])
            overlap_end = min(rm["end"], qr["end"])

            if overlap_start <= overlap_end:
                overlap_bp = overlap_end - overlap_start + 1
                # Compute RM interval in query-relative 1-based coordinates
                # to match BLAST qstart/qend (1-based inclusive)
                rm_start_rel = max(1, rm["start"] - qr["start"] + 1)
                rm_end_rel = min(qr["end"], rm["end"]) - qr["start"] + 1
                overlaps.append({
                    "region_id": qr["region_id"],
                    "chrom": rm["chrom"],
                    "overlap_start": overlap_start,
                    "overlap_end": overlap_end,
                    "overlap_bp": overlap_bp,
                    "repeat_name": rm["repeat_name"],
                    "repeat_class": rm["repeat_class"],
                    "divergence": rm.get("divergence", None),
                    "rm_start_in_query": rm_start_rel,
                    "rm_end_in_query": rm_end_rel,
                })

    return pd.DataFrame(overlaps)


def classify_hits(
    blast_hits: pd.DataFrame,
    rm_overlaps: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classify BLAST hits as known (RepeatMasker) or novel.

    A hit is "known" if its query range [qstart, qend] overlaps with
    any RepeatMasker region mapped to the same query.

    Args:
        blast_hits: BLAST results with qseqid, qstart, qend columns.
        rm_overlaps: RM overlap data with region_id, rm_start_in_query,
            rm_end_in_query, repeat_name, repeat_class.

    Returns:
        Tuple of (known_hits, novel_hits) DataFrames.
    """
    # Build lookup: region_id → list of RM intervals in query space
    rm_by_region: dict[str, list[dict]] = {}
    for _, row in rm_overlaps.iterrows():
        rid = row["region_id"]
        if rid not in rm_by_region:
            rm_by_region[rid] = []
        rm_by_region[rid].append({
            "start": row["rm_start_in_query"],
            "end": row["rm_end_in_query"],
            "repeat_name": row["repeat_name"],
            "repeat_class": row["repeat_class"],
        })

    known_rows = []
    novel_rows = []

    for _, hit in blast_hits.iterrows():
        qid = hit["qseqid"]
        regions = rm_by_region.get(qid, [])

        is_known = False
        for rm in regions:
            if hit["qstart"] <= rm["end"] and hit["qend"] >= rm["start"]:
                known_row = hit.to_dict()
                known_row["rm_repeat_name"] = rm["repeat_name"]
                known_row["rm_repeat_class"] = rm["repeat_class"]
                known_rows.append(known_row)
                is_known = True
                break

        if not is_known:
            novel_rows.append(hit.to_dict())

    known = pd.DataFrame(known_rows) if known_rows else pd.DataFrame()
    novel = pd.DataFrame(novel_rows) if novel_rows else pd.DataFrame()

    return known, novel
```

- [ ] **Step 4: Update conftest and __init__.py**

Add to `tests/conftest.py`:

```python
@pytest.fixture
def mini_repeatmasker(test_data_dir) -> Path:
    """Path to synthetic RepeatMasker .out file."""
    return test_data_dir / "mini_repeatmasker.out"
```

Update `src/fossil_finder/analysis/__init__.py` — add:

```python
from .repeatmasker import classify_hits, find_overlaps, parse_repeatmasker_out
```

And add to `__all__`:

```python
    "classify_hits",
    "find_overlaps",
    "parse_repeatmasker_out",
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_analysis/ -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/analysis/repeatmasker.py tests/test_analysis/test_repeatmasker.py \
  tests/data/mini_repeatmasker.out tests/conftest.py src/fossil_finder/analysis/__init__.py
git commit -m "feat: RepeatMasker parsing + hit classification (known vs novel)"
```

---

### Task 6: Version bump + final integration check

**Files:**
- Modify: `src/fossil_finder/__init__.py`

- [ ] **Step 1: Bump version**

Update `src/fossil_finder/__init__.py`:

```python
"""Fossil Finder: Multi-genome TE fossil mining and regulatory analysis framework."""

__version__ = "0.5.0"
```

- [ ] **Step 2: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

- [ ] **Step 3: Verify imports work cleanly**

Run: `python -c "from fossil_finder.analysis import aggregate_by_gene, classify_bias, fisher_exact_enrichment, parse_repeatmasker_out; print('Phase 5 imports OK')"`
Expected: `Phase 5 imports OK`

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/__init__.py
git commit -m "chore: bump version to 0.5.0 for Phase 5"
```

---

## Summary

| Task | Module | Tests | What it replaces |
|------|--------|-------|------------------|
| 1 | `analysis/aggregation.py` | 11 | Per-gene aggregation from `analyze_genome_wide_te.py` |
| 2 | `analysis/strand.py` | 10 | Strand bias from `analyze_genome_wide_te.py`, `rank_gene_sets.py` |
| 3 | `analysis/families.py` | 10 | TE family analysis from `analyze_te_families.py` |
| 4 | `analysis/enrichment.py` | 12 | Statistical tests from `analyze_functional_enrichment.py` |
| 5 | `analysis/repeatmasker.py` | 10 | RM overlap from `analyze_repeatmasker_overlap.py` |
| 6 | Version bump | 0 | — |

**Total new tests:** ~53
**Estimated total after Phase 5:** ~248 tests
