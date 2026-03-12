# Fossil Finder Phase 6: Pipeline Orchestrator

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the pipeline orchestrator that ties config → regions → BLAST → filter → dedup → analysis into a single, config-driven `run()` function with checkpointing and structured results.

**Architecture:** Two focused modules under `pipeline/`: (1) `steps.py` contains individual pipeline steps as pure functions (extract, blast, filter, analyze) that accept explicit inputs and return structured outputs — no global state; (2) `runner.py` provides `PipelineRunner` which reads config, wires steps together, manages output directories, and provides a `run()` entry point. Each step produces intermediate files that can be inspected or resumed from. The runner returns a `PipelineResult` dataclass containing all analysis outputs.

**Tech Stack:** Python 3.11, pandas, dataclasses, pytest

---

## File Structure (Phase 6)

```
src/fossil_finder/
├── pipeline/
│   ├── __init__.py          # Package exports
│   ├── steps.py             # Individual pipeline steps (pure functions)
│   └── runner.py            # PipelineRunner orchestrator
tests/
├── test_pipeline/
│   ├── __init__.py
│   ├── test_steps.py        # 14 tests
│   └── test_runner.py       # 10 tests
```

**Key design decisions:**

1. **Steps are pure functions** — each step takes explicit inputs (DataFrames, paths, config) and returns outputs. No side effects beyond file writes. This makes them independently testable.

2. **`PipelineRunner` is the only stateful object** — it holds config, output paths, and orchestrates step sequencing. Steps don't know about each other.

3. **`PipelineResult` is a dataclass** — not a dict. Typed fields for each output (gene_stats, strand_bias, family_stats, enrichment_results, etc.). This makes downstream code self-documenting.

4. **Checkpointing via output directory** — each step writes intermediate files (FASTA, TSV, JSON). If a file already exists and `force=False`, the step is skipped. This enables resume-from-failure.

5. **No BLAST database creation** — `makeblastdb` is an external prerequisite. The runner validates the database exists but does not create it. This avoids subprocess complexity for a one-time setup step.

6. **No reporting** — that's Phase 7. The runner returns structured data; rendering is a separate concern.

---

## Chunk 1: Pipeline Steps

### Task 1: Individual pipeline steps

Each step is a pure function that performs one stage of the pipeline. Steps are composable — the runner calls them in sequence, passing outputs as inputs to the next step.

**Files:**
- Create: `src/fossil_finder/pipeline/__init__.py`
- Create: `src/fossil_finder/pipeline/steps.py`
- Create: `tests/test_pipeline/__init__.py`
- Create: `tests/test_pipeline/test_steps.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_pipeline/__init__.py` (empty) and `tests/test_pipeline/test_steps.py`:

```python
"""Tests for individual pipeline steps."""

import json
from pathlib import Path

import pandas as pd
import pytest

from fossil_finder.pipeline.steps import (
    step_extract_regions,
    step_load_and_filter,
    step_deduplicate,
    step_aggregate,
    step_strand_analysis,
    step_family_analysis,
    step_enrichment_analysis,
    step_repeatmasker_overlap,
)


@pytest.fixture
def blast_df():
    """Minimal BLAST results DataFrame."""
    return pd.DataFrame({
        "qseqid": ["utr1", "utr1", "utr2", "utr2", "utr3"],
        "sseqid": ["TE_gypsy1", "TE_copia1", "TE_gypsy1", "TE_jockey1", "TE_gypsy1"],
        "pident": [80.0, 65.0, 85.0, 70.0, 90.0],
        "length": [120, 85, 200, 60, 150],
        "mismatch": [20, 25, 30, 15, 12],
        "gapopen": [3, 4, 5, 2, 1],
        "qstart": [50, 200, 10, 300, 100],
        "qend": [170, 284, 210, 360, 250],
        "sstart": [1, 300, 50, 100, 1],
        "send": [120, 216, 250, 160, 150],
        "evalue": [1.5e-10, 3.2e-5, 2.1e-15, 0.005, 8.3e-20],
        "bitscore": [85.2, 42.1, 120.5, 30.8, 155.0],
        "qlen": [500, 500, 800, 800, 600],
        "slen": [5000, 4500, 5000, 3000, 4000],
        "qseq": ["A"] * 5,
        "sseq": ["A"] * 5,
        "strand": ["plus", "minus", "plus", "plus", "plus"],
    })


@pytest.fixture
def query_to_gene():
    return {"utr1": "gene1", "utr2": "gene2", "utr3": "gene3"}


class TestStepExtractRegions:
    def test_extracts_to_fasta_and_metadata(self, custom_config, tmp_path):
        fasta_path = tmp_path / "regions.fa"
        meta_path = tmp_path / "regions.tsv"
        regions = step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
        )
        assert fasta_path.exists()
        assert meta_path.exists()
        assert len(regions) > 0
        assert "sequence" in regions[0]

    def test_returns_region_dicts(self, custom_config, tmp_path):
        regions = step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=tmp_path / "r.fa",
            metadata_out=tmp_path / "r.tsv",
        )
        for r in regions:
            assert "region_id" in r
            assert "chrom" in r
            assert "start" in r

    def test_skip_if_exists(self, custom_config, tmp_path):
        fasta_path = tmp_path / "regions.fa"
        meta_path = tmp_path / "regions.tsv"
        # First run
        step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
        )
        mtime = fasta_path.stat().st_mtime
        # Second run with force=False should skip
        step_extract_regions(
            config=custom_config,
            feature_type="three_prime_UTR",
            fasta_out=fasta_path,
            metadata_out=meta_path,
            force=False,
        )
        assert fasta_path.stat().st_mtime == mtime


class TestStepLoadAndFilter:
    def test_loads_from_path(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results)
        assert len(df) == 5
        assert "strand" in df.columns

    def test_applies_evalue_filter(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results, max_evalue=1e-10)
        assert all(df["evalue"] <= 1e-10)

    def test_applies_pident_filter(self, mini_blast_results):
        df = step_load_and_filter(mini_blast_results, min_pident=80.0)
        assert all(df["pident"] >= 80.0)


class TestStepDeduplicate:
    def test_removes_duplicates(self):
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 1],
            "te_end": [100, 100],
            "source_transcript": ["tr1", "tr2"],
        })
        result, stats = step_deduplicate(df)
        assert len(result) == 1
        assert stats["duplicates_removed"] == 1

    def test_no_duplicates_passthrough(self, blast_df):
        # blast_df lacks default dedup key columns, so specify BLAST-level keys
        result, stats = step_deduplicate(
            blast_df,
            key_columns=["qseqid", "sseqid", "qstart", "qend"],
        )
        assert stats["duplicates_removed"] == 0


class TestStepAggregate:
    def test_produces_gene_stats(self, blast_df, query_to_gene):
        result = step_aggregate(blast_df, query_to_gene)
        assert "gene1" in result.index
        assert "density" in result.columns
        assert result.loc["gene1", "hit_count"] == 2

    def test_empty_input(self):
        df = pd.DataFrame(columns=["qseqid", "sseqid", "length", "strand", "qlen"])
        result = step_aggregate(df, {})
        assert len(result) == 0


class TestStepStrandAnalysis:
    def test_returns_three_levels(self, blast_df):
        blast_df["gene_id"] = blast_df["qseqid"].map(
            {"utr1": "g1", "utr2": "g2", "utr3": "g3"}
        )
        result = step_strand_analysis(blast_df)
        assert "gene" in result
        assert "te_family" in result
        assert "genome" in result


class TestStepFamilyAnalysis:
    def test_returns_family_stats(self, blast_df):
        result = step_family_analysis(blast_df)
        assert "family_stats" in result
        assert "TE_gypsy1" in result["family_stats"].index

    def test_class_distribution_with_metadata(self, blast_df):
        te_meta = {"TE_gypsy1": {"te_class": "LTR"}, "TE_copia1": {"te_class": "LTR"}}
        result = step_family_analysis(blast_df, te_metadata=te_meta)
        assert "class_distribution" in result
        assert "LTR" in result["class_distribution"]


class TestStepEnrichmentAnalysis:
    def test_enrichment_for_gene_set(self):
        """Use a 4-gene setup so Mann-Whitney gets >= 2 values per group."""
        gene_stats = pd.DataFrame({
            "hit_count": [5, 3, 1, 2],
            "hit_bp": [500, 300, 100, 200],
            "query_len": [1000, 800, 600, 700],
            "density": [500.0, 375.0, 166.7, 285.7],
        }, index=["gene1", "gene2", "gene3", "gene4"])
        te_positive = {"gene1", "gene2", "gene3", "gene4"}
        result = step_enrichment_analysis(
            gene_set={"gene1", "gene2"},
            te_positive_genes=te_positive,
            gene_densities=gene_stats["density"],
        )
        assert "fisher" in result
        assert "mannwhitney" in result
        # Both groups have 2 values, so Mann-Whitney should produce a real p-value
        assert result["mannwhitney"]["p_value"] is not None


class TestStepRepeatMaskerOverlap:
    def test_classifies_hits(self, blast_df, mini_repeatmasker):
        # Need qseqid-based region mapping for RM overlap
        query_regions = pd.DataFrame({
            "region_id": ["utr1", "utr2", "utr3"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [1, 200, 50],
            "end": [500, 600, 400],
        })
        result = step_repeatmasker_overlap(
            blast_hits=blast_df,
            repeatmasker_path=mini_repeatmasker,
            query_regions=query_regions,
        )
        assert "known" in result
        assert "novel" in result
        assert "rm_stats" in result
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_pipeline/test_steps.py -v`
Expected: FAIL — `fossil_finder.pipeline.steps` does not exist

- [ ] **Step 3: Implement pipeline steps**

Create `src/fossil_finder/pipeline/__init__.py`:

```python
"""Pipeline orchestration for fossil_finder."""

from .steps import (
    step_aggregate,
    step_deduplicate,
    step_enrichment_analysis,
    step_extract_regions,
    step_family_analysis,
    step_load_and_filter,
    step_repeatmasker_overlap,
    step_strand_analysis,
)

__all__ = [
    "step_aggregate",
    "step_deduplicate",
    "step_enrichment_analysis",
    "step_extract_regions",
    "step_family_analysis",
    "step_load_and_filter",
    "step_repeatmasker_overlap",
    "step_strand_analysis",
]
```

Create `src/fossil_finder/pipeline/steps.py`:

```python
"""Individual pipeline steps for TE fossil analysis.

Each step is a pure function: explicit inputs in, structured outputs out.
Steps do not depend on each other — the runner wires them together.
"""

from pathlib import Path

import pandas as pd

from fossil_finder.analysis.aggregation import aggregate_by_gene, compute_density
from fossil_finder.analysis.enrichment import test_gene_set_enrichment
from fossil_finder.analysis.families import (
    compute_class_distribution,
    compute_family_stats,
)
from fossil_finder.analysis.repeatmasker import (
    classify_hits,
    find_overlaps,
    parse_repeatmasker_out,
)
from fossil_finder.analysis.strand import (
    compute_gene_strand_bias,
    compute_genome_strand_bias,
    compute_te_strand_bias,
)
from fossil_finder.blast.dedup import HitDeduplicator
from fossil_finder.blast.filter import apply_filters
from fossil_finder.config.schema import GenomeConfig
from fossil_finder.io.blast import load_blast_results
from fossil_finder.regions.extractor import RegionExtractor


def step_extract_regions(
    config: GenomeConfig,
    feature_type: str,
    fasta_out: str | Path,
    metadata_out: str | Path,
    gene_ids: list[str] | None = None,
    deduplicate: bool = False,
    force: bool = True,
    base_dir: Path | None = None,
) -> list[dict]:
    """Extract genomic regions and write FASTA + metadata.

    Args:
        config: Validated GenomeConfig.
        feature_type: GFF3 feature type to extract.
        fasta_out: Path to write query FASTA.
        metadata_out: Path to write region metadata TSV.
        gene_ids: Optional gene ID filter.
        deduplicate: Collapse identical sequences.
        force: If False, skip extraction when output files exist.
        base_dir: Base directory for resolving config paths.

    Returns:
        List of region dicts (same as RegionExtractor.extract_features).
    """
    fasta_out = Path(fasta_out)
    metadata_out = Path(metadata_out)

    if not force and fasta_out.exists() and metadata_out.exists():
        return []  # Skip — files already exist

    extractor = RegionExtractor(config, base_dir=base_dir)
    regions = extractor.extract_features(
        feature_type, gene_ids=gene_ids, deduplicate=deduplicate,
    )

    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    metadata_out.parent.mkdir(parents=True, exist_ok=True)

    extractor.write_fasta(regions, fasta_out)
    extractor.write_metadata(regions, metadata_out)

    return regions


def step_load_and_filter(
    blast_results: str | Path,
    max_evalue: float | None = None,
    min_pident: float | None = None,
    min_length: int | None = None,
) -> pd.DataFrame:
    """Load BLAST results and apply quality filters.

    Args:
        blast_results: Path to BLAST TSV output.
        max_evalue: Maximum e-value threshold.
        min_pident: Minimum percent identity threshold.
        min_length: Minimum alignment length threshold.

    Returns:
        Filtered BLAST results DataFrame.
    """
    df = load_blast_results(blast_results)
    return apply_filters(
        df,
        max_evalue=max_evalue,
        min_pident=min_pident,
        min_length=min_length,
    )


def step_deduplicate(
    df: pd.DataFrame,
    key_columns: list[str] | None = None,
    normalize_te_coords: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Deduplicate BLAST hits by coordinate key.

    Args:
        df: BLAST results (or annotated hits) DataFrame.
        key_columns: Columns forming the dedup key.
        normalize_te_coords: Normalize TE start/end before keying.

    Returns:
        Tuple of (deduplicated DataFrame, stats dict).
    """
    dedup = HitDeduplicator(
        key_columns=key_columns,
        normalize_te_coords=normalize_te_coords,
    )
    result = dedup.deduplicate(df)
    return result, dedup.stats


def step_aggregate(
    df: pd.DataFrame,
    query_to_gene: dict[str, str],
) -> pd.DataFrame:
    """Aggregate BLAST hits by gene and compute density.

    Args:
        df: BLAST results DataFrame.
        query_to_gene: Mapping from query ID to gene ID.

    Returns:
        Per-gene stats DataFrame with density column.
    """
    gene_stats = aggregate_by_gene(df, query_to_gene)
    return compute_density(gene_stats)


def step_strand_analysis(df: pd.DataFrame) -> dict:
    """Run strand bias analysis at all three levels.

    Args:
        df: BLAST results with gene_id and strand columns.

    Returns:
        Dict with 'gene', 'te_family', and 'genome' keys.
    """
    return {
        "gene": compute_gene_strand_bias(df),
        "te_family": compute_te_strand_bias(df),
        "genome": compute_genome_strand_bias(df),
    }


def step_family_analysis(
    df: pd.DataFrame,
    te_metadata: dict[str, dict] | None = None,
) -> dict:
    """Compute TE family statistics and optional class distribution.

    Args:
        df: BLAST results DataFrame.
        te_metadata: Optional TE ID → metadata mapping.

    Returns:
        Dict with 'family_stats' DataFrame and optional 'class_distribution'.
    """
    result = {"family_stats": compute_family_stats(df)}
    if te_metadata is not None:
        result["class_distribution"] = compute_class_distribution(df, te_metadata)
    return result


def step_enrichment_analysis(
    gene_set: set[str],
    te_positive_genes: set[str],
    gene_densities: "pd.Series",
) -> dict:
    """Run enrichment analysis for a gene set.

    Args:
        gene_set: Gene IDs in the set of interest.
        te_positive_genes: All gene IDs with at least one TE hit.
        gene_densities: Series mapping gene_id → density value.

    Returns:
        Dict with 'fisher' and 'mannwhitney' sub-dicts.
    """
    return test_gene_set_enrichment(gene_set, te_positive_genes, gene_densities)


def step_repeatmasker_overlap(
    blast_hits: pd.DataFrame,
    repeatmasker_path: str | Path,
    query_regions: pd.DataFrame,
) -> dict:
    """Classify BLAST hits as known (RepeatMasker) or novel.

    Args:
        blast_hits: BLAST results DataFrame.
        repeatmasker_path: Path to RepeatMasker .out file.
        query_regions: DataFrame with region_id, chrom, start, end.

    Returns:
        Dict with 'known' DataFrame, 'novel' DataFrame,
        and 'rm_stats' summary dict.
    """
    rm_regions = parse_repeatmasker_out(repeatmasker_path)
    overlaps = find_overlaps(rm_regions, query_regions)

    if overlaps.empty:
        return {
            "known": pd.DataFrame(),
            "novel": blast_hits.copy(),
            "rm_stats": {
                "total_rm_regions": len(rm_regions),
                "overlapping_regions": 0,
                "known_hits": 0,
                "novel_hits": len(blast_hits),
            },
        }

    known, novel = classify_hits(blast_hits, overlaps)

    return {
        "known": known,
        "novel": novel,
        "rm_stats": {
            "total_rm_regions": len(rm_regions),
            "overlapping_regions": len(overlaps),
            "known_hits": len(known),
            "novel_hits": len(novel),
        },
    }
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_pipeline/test_steps.py -v`
Expected: ALL PASS

- [ ] **Step 5: Run full suite**

Run: `pytest tests/ --tb=short -q`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/pipeline/__init__.py src/fossil_finder/pipeline/steps.py \
  tests/test_pipeline/__init__.py tests/test_pipeline/test_steps.py
git commit -m "feat: pipeline steps — composable functions for each analysis stage"
```

---

## Chunk 2: Pipeline Runner

### Task 2: PipelineResult dataclass and PipelineRunner

The runner wires steps together, manages output paths, and returns a typed result.

**Files:**
- Create: `src/fossil_finder/pipeline/runner.py`
- Create: `tests/test_pipeline/test_runner.py`
- Modify: `src/fossil_finder/pipeline/__init__.py` (add exports)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_pipeline/test_runner.py`:

```python
"""Tests for PipelineRunner orchestration."""

import json
from pathlib import Path

import pandas as pd
import pytest

from fossil_finder.pipeline.runner import PipelineResult, PipelineRunner


@pytest.fixture
def pipeline_config(test_data_dir):
    """Config with blast results pre-computed (no BLAST+ needed)."""
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
            "te_consensus": str(test_data_dir / "mini_tes.fasta"),
        },
        blast={"word_size": 7, "dust": False},
    )


class TestPipelineResult:
    def test_default_fields(self):
        result = PipelineResult()
        assert result.gene_stats is None
        assert result.strand_bias is None
        assert result.family_stats is None
        assert result.enrichment == {}
        assert result.rm_overlap is None
        assert result.blast_hits is None

    def test_summary_dict(self):
        result = PipelineResult()
        result.gene_stats = pd.DataFrame(
            {"hit_count": [5, 3]}, index=["g1", "g2"]
        )
        summary = result.summary()
        assert "n_genes_analyzed" in summary
        assert summary["n_genes_analyzed"] == 2


class TestPipelineRunner:
    def test_init_creates_output_dir(self, pipeline_config, tmp_path):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        assert (tmp_path / "output").exists()

    def test_extract_step(self, pipeline_config, tmp_path):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        regions = runner.extract(feature_type="three_prime_UTR")
        assert len(regions) > 0
        assert (tmp_path / "output" / "regions.fa").exists()
        assert (tmp_path / "output" / "regions.tsv").exists()

    def test_analyze_with_precomputed_blast(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        """Full analysis using pre-computed BLAST results (no BLAST+ needed)."""
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        assert isinstance(result, PipelineResult)
        assert result.gene_stats is not None
        assert len(result.gene_stats) > 0
        assert result.strand_bias is not None
        assert result.family_stats is not None

    def test_analyze_empty_blast(self, pipeline_config, tmp_path):
        """Empty BLAST results should produce empty but valid result."""
        empty_blast = tmp_path / "empty.tsv"
        empty_blast.write_text("")
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=empty_blast,
            query_to_gene={},
        )
        assert isinstance(result, PipelineResult)
        assert result.gene_stats is not None
        assert len(result.gene_stats) == 0

    def test_analyze_with_gene_sets(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        gene_sets = {
            "test_set": {"gene1", "gene2"},
        }
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
            gene_sets=gene_sets,
        )
        assert "test_set" in result.enrichment

    def test_analyze_with_repeatmasker(
        self, pipeline_config, mini_blast_results, mini_repeatmasker, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        query_regions = pd.DataFrame({
            "region_id": ["gene1_utr", "gene2_utr", "gene3_utr"],
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [281, 350, 381],
            "end": [300, 369, 400],
        })
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
            repeatmasker_path=mini_repeatmasker,
            query_regions=query_regions,
        )
        assert result.rm_overlap is not None
        assert "known" in result.rm_overlap
        assert "novel" in result.rm_overlap

    def test_saves_summary_json(
        self, pipeline_config, mini_blast_results, tmp_path
    ):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        result = runner.analyze(
            blast_results=mini_blast_results,
            query_to_gene={"gene1_utr": "gene1", "gene2_utr": "gene2",
                           "gene3_utr": "gene3"},
        )
        summary_path = tmp_path / "output" / "summary.json"
        assert summary_path.exists()
        summary = json.loads(summary_path.read_text())
        assert "n_genes_analyzed" in summary

    def test_force_rerun(self, pipeline_config, mini_blast_results, tmp_path):
        runner = PipelineRunner(
            config=pipeline_config,
            output_dir=tmp_path / "output",
        )
        q2g = {"gene1_utr": "gene1", "gene2_utr": "gene2",
               "gene3_utr": "gene3"}
        result1 = runner.analyze(blast_results=mini_blast_results,
                                 query_to_gene=q2g)
        result2 = runner.analyze(blast_results=mini_blast_results,
                                 query_to_gene=q2g, force=True)
        assert isinstance(result2, PipelineResult)
        assert result2.gene_stats is not None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_pipeline/test_runner.py -v`
Expected: FAIL — `fossil_finder.pipeline.runner` does not exist

- [ ] **Step 3: Implement PipelineResult and PipelineRunner**

Create `src/fossil_finder/pipeline/runner.py`:

```python
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
        force: bool = True,
    ) -> PipelineResult:
        """Run full analysis on BLAST results.

        Args:
            blast_results: Path to BLAST TSV output.
            query_to_gene: Mapping from query ID to gene ID.
            max_evalue: E-value filter threshold.
            min_pident: Percent identity filter threshold.
            min_length: Alignment length filter threshold.
            gene_sets: Named gene sets for enrichment testing.
            te_metadata: TE ID → metadata mapping for class distribution.
            repeatmasker_path: Path to RepeatMasker .out file.
            query_regions: DataFrame for RM overlap (region_id, chrom, start, end).
            force: Re-run even if outputs exist.

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
```

- [ ] **Step 4: Update `__init__.py` exports**

Update `src/fossil_finder/pipeline/__init__.py`:

```python
"""Pipeline orchestration for fossil_finder."""

from .runner import PipelineResult, PipelineRunner
from .steps import (
    step_aggregate,
    step_deduplicate,
    step_enrichment_analysis,
    step_extract_regions,
    step_family_analysis,
    step_load_and_filter,
    step_repeatmasker_overlap,
    step_strand_analysis,
)

__all__ = [
    "PipelineResult",
    "PipelineRunner",
    "step_aggregate",
    "step_deduplicate",
    "step_enrichment_analysis",
    "step_extract_regions",
    "step_family_analysis",
    "step_load_and_filter",
    "step_repeatmasker_overlap",
    "step_strand_analysis",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_pipeline/ -v`
Expected: ALL PASS

- [ ] **Step 6: Run full suite**

Run: `pytest tests/ --tb=short -q`
Expected: ALL PASS

- [ ] **Step 7: Commit**

```bash
git add src/fossil_finder/pipeline/runner.py tests/test_pipeline/test_runner.py \
  src/fossil_finder/pipeline/__init__.py
git commit -m "feat: PipelineRunner — config-driven orchestration with typed results"
```

---

## Chunk 3: Version Bump + Integration Verification

### Task 3: Version bump and integration check

**Files:**
- Modify: `src/fossil_finder/__init__.py`

- [ ] **Step 1: Bump version**

Update `src/fossil_finder/__init__.py`:

```python
"""Fossil Finder: Multi-genome TE fossil mining and regulatory analysis framework."""

__version__ = "0.6.0"
```

- [ ] **Step 2: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

- [ ] **Step 3: Verify imports work cleanly**

Run: `python -c "from fossil_finder.pipeline import PipelineRunner, PipelineResult, step_extract_regions, step_aggregate; print('Phase 6 imports OK')"`
Expected: `Phase 6 imports OK`

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/__init__.py
git commit -m "chore: bump version to 0.6.0 for Phase 6"
```

---

## Summary

| Task | Module | Tests | What it provides |
|------|--------|-------|------------------|
| 1 | `pipeline/steps.py` | 14 | Composable step functions for each pipeline stage |
| 2 | `pipeline/runner.py` | 10 | `PipelineRunner` orchestrator + `PipelineResult` dataclass |
| 3 | Version bump | 0 | Version 0.6.0 |

**Total new tests:** ~24
**Estimated total after Phase 6:** ~274 tests

**What Phase 6 does NOT include (deferred to later phases):**
- `makeblastdb` — external prerequisite, not orchestrated
- HTML/Markdown reporting — Phase 7
- CLI entry point — Phase 8
- Shuffled-sequence controls — could be a future step function
