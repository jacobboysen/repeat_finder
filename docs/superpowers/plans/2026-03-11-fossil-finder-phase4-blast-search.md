# Fossil Finder Phase 4: BLAST Search & Hit Deduplication

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a config-driven BLAST runner and coordinate-based hit deduplication pipeline, replacing the legacy `scripts/blast_runner.py` and `scripts/deduplicate_*_te_hits.py` with generic, multi-genome modules.

**Architecture:** Three focused modules: (1) `BlastRunner` builds and executes BLAST commands from `GenomeConfig.blast` (BlastSpec), returning parsed results; (2) `HitFilter` applies quality thresholds (e-value, identity, length) as a composable pipeline; (3) `HitDeduplicator` removes duplicate hits using a configurable coordinate-based key, replacing the legacy 7-tuple `(gene_id, chrom, genomic_start, genomic_end, te_id, te_start, te_end)` dedup. All three modules work with the existing `io/blast.py` parsing layer — no duplication of I/O code.

**Tech Stack:** Python 3.11, subprocess (BLAST+), pandas, Pydantic v2 (config), pytest

---

## File Structure (Phase 4)

```
src/fossil_finder/
├── blast/
│   ├── __init__.py          # Package exports: BlastRunner, HitFilter, HitDeduplicator
│   ├── runner.py            # BlastRunner — builds command, executes, returns results
│   ├── filter.py            # HitFilter — composable quality filtering pipeline
│   └── dedup.py             # HitDeduplicator — coordinate-based hit deduplication
├── config/
│   └── schema.py            # MODIFY: add max_hsps, program fields to BlastSpec
tests/
├── test_blast/
│   ├── __init__.py
│   ├── test_runner.py       # 12 tests: command building, execution, error handling
│   ├── test_filter.py       # 10 tests: threshold filtering, chaining, edge cases
│   └── test_dedup.py        # 12 tests: dedup key, stats, coordinate normalization
├── data/
│   └── mini_blast_results.tsv  # EXISTS — 5 hits, reused as-is
```

**Key design decisions:**

1. **`blast/` subpackage** (not `execution/`): mirrors the `te/`, `regions/`, `adapters/` naming convention — named for the domain, not the action.

2. **`BlastRunner` does NOT shell out in tests**: The actual `blastn` subprocess call is isolated behind a `_run_subprocess()` method so tests can verify command construction without requiring BLAST+ installed. One integration test with `shutil.which("blastn")` guard validates end-to-end.

3. **`HitFilter` is stateless**: Pure functions that accept a DataFrame and return a filtered DataFrame. No classes needed — just `apply_filters(df, ...)` with keyword args.

4. **`HitDeduplicator` is generic**: The dedup key is configurable (default: the legacy 7-tuple), but the user can override it for different analysis types. Stats tracking built in.

5. **No RepeatMasker or domain annotation here**: Those belong in Phase 5 (analysis modules). This phase is search + dedup only.

---

## Chunk 1: BlastRunner + Config Updates

### Task 1: Extend BlastSpec with missing fields

The existing `BlastSpec` in `config/schema.py` is missing `max_hsps` and `program` (both used by the legacy runner). Add them.

**Files:**
- Modify: `src/fossil_finder/config/schema.py:57-73`
- Test: `tests/test_config/test_schema.py` (add cases to existing file)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_config/test_schema.py`:

```python
class TestBlastSpecExtended:
    def test_max_hsps_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.max_hsps == 10

    def test_program_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.program == "blastn"

    def test_program_override(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec(program="tblastx")
        assert spec.program == "tblastx"

    def test_soft_masking_default(self):
        from fossil_finder.config.schema import BlastSpec

        spec = BlastSpec()
        assert spec.soft_masking is True

    def test_blast_spec_in_config(self, mini_genome_config):
        from fossil_finder.config.schema import load_genome_config

        config = load_genome_config(mini_genome_config)
        assert config.blast.program == "blastn"
        assert config.blast.max_hsps == 10
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config/test_schema.py::TestBlastSpecExtended -v`
Expected: FAIL — `BlastSpec` has no `max_hsps` or `program` attribute

- [ ] **Step 3: Implement — add fields to BlastSpec**

In `src/fossil_finder/config/schema.py`, update the `BlastSpec` class:

```python
class BlastSpec(BaseModel):
    """BLAST search parameters.

    Defaults based on DUST_FILTERING_ANALYSIS.md (2026-01-22):
    dust=no with stringent e-value captures 52% more high-quality hits
    and higher real/shuffled enrichment ratios than dust=yes.
    """

    program: str = "blastn"
    word_size: int = 7
    gapopen: int = 2
    gapextend: int = 1
    penalty: int = -1
    reward: int = 1
    dust: bool = False
    soft_masking: bool = True
    evalue: float = 0.001
    max_target_seqs: int = 1000
    max_hsps: int = 10
    num_threads: int = 4
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_config/test_schema.py -v`
Expected: ALL PASS (existing tests unaffected by new defaults)

- [ ] **Step 5: Commit**

```bash
git add src/fossil_finder/config/schema.py tests/test_config/test_schema.py
git commit -m "feat(config): add program, max_hsps, soft_masking to BlastSpec"
```

---

### Task 2: BlastRunner — command construction

The core of the runner: build a BLAST command-line from `BlastSpec` + query/database paths.

**Files:**
- Create: `src/fossil_finder/blast/__init__.py`
- Create: `src/fossil_finder/blast/runner.py`
- Create: `tests/test_blast/__init__.py`
- Create: `tests/test_blast/test_runner.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_blast/__init__.py` (empty) and `tests/test_blast/test_runner.py`:

```python
"""Tests for BlastRunner."""

import pytest

from fossil_finder.blast.runner import BlastRunner
from fossil_finder.config.schema import BlastSpec


class TestBuildCommand:
    """Test BLAST command construction from BlastSpec."""

    def test_default_command(self, tmp_path):
        query = tmp_path / "query.fasta"
        query.write_text(">seq1\nATCG\n")
        db = tmp_path / "db"
        output = tmp_path / "results.tsv"

        runner = BlastRunner(BlastSpec())
        cmd = runner.build_command(query, db, output)

        assert cmd[0] == "blastn"
        assert "-query" in cmd
        assert str(query) in cmd
        assert "-db" in cmd
        assert str(db) in cmd
        assert "-out" in cmd
        assert str(output) in cmd

    def test_evalue_in_command(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        runner = BlastRunner(BlastSpec(evalue=10.0))
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")
        idx = cmd.index("-evalue")
        assert cmd[idx + 1] == "10.0"

    def test_dust_no_in_command(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        runner = BlastRunner(BlastSpec(dust=False))
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")
        idx = cmd.index("-dust")
        assert cmd[idx + 1] == "no"

    def test_dust_yes_in_command(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        runner = BlastRunner(BlastSpec(dust=True))
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")
        idx = cmd.index("-dust")
        assert cmd[idx + 1] == "yes"

    def test_outfmt_includes_17_columns(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        runner = BlastRunner(BlastSpec())
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")
        idx = cmd.index("-outfmt")
        outfmt_str = cmd[idx + 1]
        assert outfmt_str.startswith("6 ")
        fields = outfmt_str.split()[1:]  # skip the "6"
        assert "qseqid" in fields
        assert "sseqid" in fields
        assert "qseq" in fields
        assert "sseq" in fields
        # strand is computed post-hoc, not in outfmt
        assert "strand" not in fields

    def test_custom_program(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        runner = BlastRunner(BlastSpec(program="tblastx"))
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")
        assert cmd[0] == "tblastx"

    def test_word_size_and_gap_params(self, tmp_path):
        query = tmp_path / "q.fasta"
        query.write_text(">s\nA\n")

        spec = BlastSpec(word_size=11, gapopen=5, gapextend=2)
        runner = BlastRunner(spec)
        cmd = runner.build_command(query, tmp_path / "db", tmp_path / "out.tsv")

        idx_ws = cmd.index("-word_size")
        assert cmd[idx_ws + 1] == "11"
        idx_go = cmd.index("-gapopen")
        assert cmd[idx_go + 1] == "5"
        idx_ge = cmd.index("-gapextend")
        assert cmd[idx_ge + 1] == "2"

    def test_query_not_found_raises(self, tmp_path):
        runner = BlastRunner(BlastSpec())
        with pytest.raises(FileNotFoundError, match="Query file not found"):
            runner.build_command(
                tmp_path / "nonexistent.fasta",
                tmp_path / "db",
                tmp_path / "out.tsv",
            )


class TestRunBlast:
    """Test BLAST execution (mocked subprocess)."""

    def test_run_returns_result_path(self, tmp_path, monkeypatch):
        """Verify run() calls subprocess and returns output path."""
        query = tmp_path / "query.fasta"
        query.write_text(">seq1\nATCGATCG\n")
        db = tmp_path / "db"
        output = tmp_path / "results.tsv"

        # Write fake BLAST output so load_blast_results works
        output.write_text(
            "seq1\tTE1\t80.0\t100\t15\t2\t1\t100\t1\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG\n"
        )

        import subprocess

        calls = []

        def mock_run(cmd, **kwargs):
            calls.append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        monkeypatch.setattr("subprocess.run", mock_run)

        runner = BlastRunner(BlastSpec())
        result_path = runner.run(query, db, output)
        assert result_path == output
        assert len(calls) == 1
        assert calls[0][0] == "blastn"

    def test_run_raises_on_blast_failure(self, tmp_path, monkeypatch):
        query = tmp_path / "query.fasta"
        query.write_text(">seq1\nATCG\n")

        import subprocess

        def mock_run(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd, stderr="BLAST error")

        monkeypatch.setattr("subprocess.run", mock_run)

        runner = BlastRunner(BlastSpec())
        with pytest.raises(subprocess.CalledProcessError):
            runner.run(query, tmp_path / "db", tmp_path / "out.tsv")

    def test_run_raises_on_missing_blast(self, tmp_path, monkeypatch):
        query = tmp_path / "query.fasta"
        query.write_text(">seq1\nATCG\n")

        def mock_run(cmd, **kwargs):
            raise FileNotFoundError("blastn not found")

        monkeypatch.setattr("subprocess.run", mock_run)

        runner = BlastRunner(BlastSpec())
        with pytest.raises(FileNotFoundError):
            runner.run(query, tmp_path / "db", tmp_path / "out.tsv")

    def test_run_and_load_parses_results(self, tmp_path, monkeypatch):
        """run_and_load() should return a parsed DataFrame."""
        query = tmp_path / "query.fasta"
        query.write_text(">seq1\nATCGATCG\n")
        output = tmp_path / "results.tsv"

        # Fake BLAST output: 16 columns (no strand), strand computed post-hoc
        output.write_text(
            "seq1\tTE1\t80.0\t100\t15\t2\t1\t100\t1\t100\t1e-5\t50.0\t500\t3000\tATCG\tATCG\n"
        )

        import subprocess

        def mock_run(cmd, **kwargs):
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        monkeypatch.setattr("subprocess.run", mock_run)

        runner = BlastRunner(BlastSpec())
        df = runner.run_and_load(query, tmp_path / "db", output)
        assert len(df) == 1
        assert "strand" in df.columns
        assert df.iloc[0]["strand"] == "plus"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_blast/test_runner.py -v`
Expected: FAIL — `fossil_finder.blast.runner` does not exist

- [ ] **Step 3: Implement BlastRunner**

Create `src/fossil_finder/blast/__init__.py`:

```python
"""BLAST search execution and hit processing."""

from .runner import BlastRunner

__all__ = ["BlastRunner"]
```

Create `src/fossil_finder/blast/runner.py`:

```python
"""Config-driven BLAST command execution.

Builds BLAST command-lines from BlastSpec, executes them via subprocess,
and returns parsed results using the existing io/blast.py layer.
"""

import logging
import subprocess
from pathlib import Path

from fossil_finder.config.schema import BlastSpec
from fossil_finder.io.blast import BLAST_COLUMNS_NO_STRAND, load_blast_results

import pandas as pd

logger = logging.getLogger(__name__)

# Derive outfmt fields from the single source of truth in io/blast.py
_OUTFMT_FIELDS = " ".join(BLAST_COLUMNS_NO_STRAND)


class BlastRunner:
    """Execute BLAST searches using BlastSpec configuration.

    Wraps subprocess BLAST execution with:
    - Command construction from BlastSpec parameters
    - Consistent 16-column outfmt (strand computed post-hoc)
    - Result parsing via fossil_finder.io.blast
    """

    def __init__(self, spec: BlastSpec):
        self.spec = spec

    def build_command(
        self,
        query: str | Path,
        database: str | Path,
        output: str | Path,
    ) -> list[str]:
        """Build BLAST command-line from spec + paths.

        Args:
            query: Path to query FASTA file (must exist).
            database: Path to BLAST database (without extension).
            output: Path for output TSV file.

        Returns:
            Command as list of strings for subprocess.run().

        Raises:
            FileNotFoundError: If query file does not exist.
        """
        query = Path(query)
        if not query.exists():
            raise FileNotFoundError(f"Query file not found: {query}")

        s = self.spec
        cmd = [
            s.program,
            "-query", str(query),
            "-db", str(database),
            "-out", str(output),
            "-evalue", str(s.evalue),
            "-word_size", str(s.word_size),
            "-gapopen", str(s.gapopen),
            "-gapextend", str(s.gapextend),
            "-penalty", str(s.penalty),
            "-reward", str(s.reward),
            "-dust", "yes" if s.dust else "no",
            "-soft_masking", str(s.soft_masking).lower(),
            "-max_target_seqs", str(s.max_target_seqs),
            "-max_hsps", str(s.max_hsps),
            "-num_threads", str(s.num_threads),
            "-outfmt", f"6 {_OUTFMT_FIELDS}",
        ]

        return cmd

    def run(
        self,
        query: str | Path,
        database: str | Path,
        output: str | Path,
    ) -> Path:
        """Execute BLAST and return output file path.

        Args:
            query: Path to query FASTA.
            database: Path to BLAST database.
            output: Path for results TSV.

        Returns:
            Path to the results file.

        Raises:
            FileNotFoundError: If query or BLAST binary not found.
            subprocess.CalledProcessError: If BLAST exits non-zero.
        """
        output = Path(output)
        cmd = self.build_command(query, database, output)

        logger.info("Running BLAST: %s", " ".join(cmd[:3]))
        logger.debug("Full command: %s", " ".join(cmd))

        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
        )

        if result.stderr:
            logger.warning("BLAST stderr: %s", result.stderr.strip())

        return output

    def run_and_load(
        self,
        query: str | Path,
        database: str | Path,
        output: str | Path,
    ) -> pd.DataFrame:
        """Execute BLAST, parse results, and return DataFrame.

        Convenience method combining run() + load_blast_results().
        """
        result_path = self.run(query, database, output)
        return load_blast_results(result_path)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_blast/test_runner.py -v`
Expected: ALL PASS

- [ ] **Step 5: Run full suite to verify no regressions**

Run: `pytest tests/ -v --tb=short`
Expected: ALL 155+ tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/blast/__init__.py src/fossil_finder/blast/runner.py \
  tests/test_blast/__init__.py tests/test_blast/test_runner.py \
  src/fossil_finder/config/schema.py tests/test_config/test_schema.py
git commit -m "feat: BlastRunner — config-driven BLAST command execution"
```

---

### Task 3: BlastRunner integration test (optional, requires BLAST+)

A single integration test that actually runs `blastn` against the mini test data, gated by `shutil.which("blastn")`.

**Files:**
- Modify: `tests/test_blast/test_runner.py` (append)

- [ ] **Step 1: Write the integration test**

Append to `tests/test_blast/test_runner.py`:

```python
import shutil
import subprocess


@pytest.mark.skipif(
    shutil.which("blastn") is None,
    reason="BLAST+ not installed",
)
class TestBlastRunnerIntegration:
    """End-to-end test requiring BLAST+ installed."""

    def test_blast_against_mini_tes(self, tmp_path, test_data_dir):
        """Run blastn with mini genome queries against mini TE database."""
        from fossil_finder.io.fasta import write_fasta

        # Create a tiny query from the test genome
        query = tmp_path / "query.fasta"
        write_fasta({"test_query": "ATCGATCGATCGATCGATCG"}, query)

        # Build a BLAST database from mini TEs
        te_fasta = test_data_dir / "mini_tes.fasta"
        db_path = tmp_path / "test_db"
        subprocess.run(
            ["makeblastdb", "-in", str(te_fasta), "-dbtype", "nucl",
             "-out", str(db_path)],
            check=True, capture_output=True,
        )

        # Run BLAST
        output = tmp_path / "results.tsv"
        runner = BlastRunner(BlastSpec(evalue=10, word_size=7))
        df = runner.run_and_load(query, db_path, output)

        # Should find hits (the mini genome shares sequence with mini TEs)
        assert len(df) > 0
        assert "qseqid" in df.columns
        assert "strand" in df.columns
```

- [ ] **Step 2: Run test (only passes if BLAST+ installed)**

Run: `pytest tests/test_blast/test_runner.py::TestBlastRunnerIntegration -v`
Expected: PASS if BLAST+ installed, SKIP otherwise

- [ ] **Step 3: Commit**

```bash
git add tests/test_blast/test_runner.py
git commit -m "test: add BLAST+ integration test (skipped without blastn)"
```

---

## Chunk 2: HitFilter

### Task 4: HitFilter — composable quality filtering

Pure functions that filter a BLAST results DataFrame by quality thresholds. This replaces ad-hoc filtering scattered across legacy analysis scripts.

**Note:** `io/blast.py`'s `load_blast_results()` already accepts `min_length`, `min_pident`, and `max_evalue` parameters for load-time filtering. Those remain for convenience. This filter module is the preferred approach for analysis pipelines where filters are applied post-load, chained, or varied across experiments.

**Files:**
- Create: `src/fossil_finder/blast/filter.py`
- Create: `tests/test_blast/test_filter.py`
- Modify: `src/fossil_finder/blast/__init__.py` (add export)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_blast/test_filter.py`:

```python
"""Tests for BLAST hit filtering."""

import pandas as pd
import pytest

from fossil_finder.blast.filter import apply_filters, filter_by_evalue, filter_by_pident, filter_by_length


@pytest.fixture
def sample_hits():
    """DataFrame with varied quality hits for filter testing."""
    return pd.DataFrame({
        "qseqid": ["q1", "q2", "q3", "q4", "q5"],
        "sseqid": ["s1", "s1", "s2", "s2", "s3"],
        "pident": [95.0, 65.0, 80.0, 45.0, 90.0],
        "length": [200, 50, 100, 30, 150],
        "evalue": [1e-20, 0.5, 1e-5, 10.0, 1e-10],
        "bitscore": [150.0, 20.0, 80.0, 10.0, 120.0],
        "strand": ["plus", "minus", "plus", "plus", "minus"],
    })


class TestFilterByEvalue:
    def test_strict_evalue(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=1e-5)
        assert len(result) == 3  # q1, q3, q5
        assert all(result["evalue"] <= 1e-5)

    def test_permissive_evalue(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=10.0)
        assert len(result) == 5  # all pass

    def test_no_hits_pass(self, sample_hits):
        result = filter_by_evalue(sample_hits, max_evalue=1e-30)
        assert len(result) == 0


class TestFilterByPident:
    def test_high_identity(self, sample_hits):
        result = filter_by_pident(sample_hits, min_pident=80.0)
        assert len(result) == 3  # q1, q3, q5

    def test_low_threshold_passes_all(self, sample_hits):
        result = filter_by_pident(sample_hits, min_pident=0.0)
        assert len(result) == 5


class TestFilterByLength:
    def test_min_length(self, sample_hits):
        result = filter_by_length(sample_hits, min_length=100)
        assert len(result) == 3  # q1, q3, q5

    def test_length_zero_passes_all(self, sample_hits):
        result = filter_by_length(sample_hits, min_length=0)
        assert len(result) == 5


class TestApplyFilters:
    def test_combined_filters(self, sample_hits):
        result = apply_filters(
            sample_hits, max_evalue=1e-5, min_pident=80.0, min_length=100,
        )
        assert len(result) == 3  # q1, q3, q5

    def test_no_filters_returns_copy(self, sample_hits):
        result = apply_filters(sample_hits)
        assert len(result) == 5
        assert result is not sample_hits  # returns copy, not same object

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["qseqid", "evalue", "pident", "length"])
        result = apply_filters(df, max_evalue=1e-5)
        assert len(result) == 0

    def test_filters_are_composable(self, sample_hits):
        """Applying filters sequentially equals applying all at once."""
        step1 = filter_by_evalue(sample_hits, max_evalue=0.5)
        step2 = filter_by_pident(step1, min_pident=80.0)

        combined = apply_filters(sample_hits, max_evalue=0.5, min_pident=80.0)
        assert list(combined["qseqid"]) == list(step2["qseqid"])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_blast/test_filter.py -v`
Expected: FAIL — `fossil_finder.blast.filter` does not exist

- [ ] **Step 3: Implement HitFilter**

Create `src/fossil_finder/blast/filter.py`:

```python
"""Composable BLAST hit quality filtering.

Pure functions that accept a DataFrame and return a filtered copy.
Replaces ad-hoc filtering scattered across legacy analysis scripts.
"""

import pandas as pd


def filter_by_evalue(df: pd.DataFrame, max_evalue: float) -> pd.DataFrame:
    """Keep hits with e-value <= threshold."""
    return df[df["evalue"] <= max_evalue].copy()


def filter_by_pident(df: pd.DataFrame, min_pident: float) -> pd.DataFrame:
    """Keep hits with percent identity >= threshold."""
    return df[df["pident"] >= min_pident].copy()


def filter_by_length(df: pd.DataFrame, min_length: int) -> pd.DataFrame:
    """Keep hits with alignment length >= threshold."""
    return df[df["length"] >= min_length].copy()


def apply_filters(
    df: pd.DataFrame,
    max_evalue: float | None = None,
    min_pident: float | None = None,
    min_length: int | None = None,
) -> pd.DataFrame:
    """Apply multiple quality filters in sequence.

    Args:
        df: BLAST results DataFrame.
        max_evalue: Maximum e-value (inclusive).
        min_pident: Minimum percent identity (inclusive).
        min_length: Minimum alignment length (inclusive).

    Returns:
        Filtered copy of the DataFrame.
    """
    result = df.copy()

    if max_evalue is not None:
        result = result[result["evalue"] <= max_evalue]
    if min_pident is not None:
        result = result[result["pident"] >= min_pident]
    if min_length is not None:
        result = result[result["length"] >= min_length]

    return result
```

- [ ] **Step 4: Update `blast/__init__.py` exports**

```python
"""BLAST search execution and hit processing."""

from .filter import apply_filters, filter_by_evalue, filter_by_length, filter_by_pident
from .runner import BlastRunner

__all__ = [
    "BlastRunner",
    "apply_filters",
    "filter_by_evalue",
    "filter_by_length",
    "filter_by_pident",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_blast/test_filter.py -v`
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/fossil_finder/blast/filter.py src/fossil_finder/blast/__init__.py \
  tests/test_blast/test_filter.py
git commit -m "feat: HitFilter — composable BLAST quality filtering"
```

---

## Chunk 3: HitDeduplicator

### Task 5: HitDeduplicator — coordinate-based dedup

Replaces the legacy `deduplicate_*_te_hits.py` scripts. Uses a configurable key tuple to identify duplicate hits. Default key matches the legacy 7-tuple: `(gene_id, chrom, genomic_start, genomic_end, te_id, te_start, te_end)`.

**Context for the implementor:** The legacy dedup scripts work on _annotated_ hits that already have genomic coordinates and gene IDs. The raw BLAST output has _query-relative_ coordinates. Coordinate conversion (query → genomic) belongs in Phase 5 analysis modules. This deduplicator works on already-annotated DataFrames with the right columns present.

**Files:**
- Create: `src/fossil_finder/blast/dedup.py`
- Create: `tests/test_blast/test_dedup.py`
- Modify: `src/fossil_finder/blast/__init__.py` (add export)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_blast/test_dedup.py`:

```python
"""Tests for BLAST hit deduplication."""

import pandas as pd
import pytest

from fossil_finder.blast.dedup import HitDeduplicator


DEFAULT_KEY = ["gene_id", "chrom", "genomic_start", "genomic_end",
               "te_id", "te_start", "te_end"]


@pytest.fixture
def annotated_hits():
    """DataFrame simulating annotated BLAST hits with genomic coordinates."""
    return pd.DataFrame({
        "gene_id": ["g1", "g1", "g1", "g2", "g2"],
        "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "genomic_start": [100, 100, 200, 500, 500],
        "genomic_end": [200, 200, 300, 600, 600],
        "te_id": ["TE_gypsy1", "TE_gypsy1", "TE_gypsy1", "TE_jockey1", "TE_jockey1"],
        "te_start": [1, 1, 50, 1, 1],
        "te_end": [100, 100, 150, 100, 100],
        "pident": [80.0, 80.0, 75.0, 90.0, 90.0],
        "evalue": [1e-10, 1e-10, 1e-5, 1e-15, 1e-15],
        "source_transcript": ["tr1", "tr2", "tr1", "tr3", "tr4"],
    })


class TestHitDeduplicator:
    def test_removes_exact_duplicates(self, annotated_hits):
        dedup = HitDeduplicator()
        result = dedup.deduplicate(annotated_hits)
        # Rows 0 and 1 are duplicates (same key, from different transcripts)
        # Rows 3 and 4 are duplicates
        assert len(result) == 3  # rows 0, 2, 3 (or 1, 2, 4)

    def test_keeps_different_te_regions(self, annotated_hits):
        dedup = HitDeduplicator()
        result = dedup.deduplicate(annotated_hits)
        # g1 has two distinct TE regions: (1,100) and (50,150)
        g1_hits = result[result["gene_id"] == "g1"]
        assert len(g1_hits) == 2

    def test_stats_tracking(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        assert dedup.stats["total_input"] == 5
        assert dedup.stats["unique_hits"] == 3
        assert dedup.stats["duplicates_removed"] == 2

    def test_custom_key_columns(self, annotated_hits):
        """Using a smaller key should produce more dedup."""
        dedup = HitDeduplicator(key_columns=["gene_id", "chrom", "te_id"])
        result = dedup.deduplicate(annotated_hits)
        # With this key: g1+chr1+TE_gypsy1 = 1 unique, g2+chr2+TE_jockey1 = 1 unique
        assert len(result) == 2

    def test_empty_dataframe(self):
        dedup = HitDeduplicator()
        df = pd.DataFrame(columns=DEFAULT_KEY + ["pident"])
        result = dedup.deduplicate(df)
        assert len(result) == 0
        assert dedup.stats["total_input"] == 0

    def test_no_duplicates(self):
        df = pd.DataFrame({
            "gene_id": ["g1", "g2"],
            "chrom": ["chr1", "chr2"],
            "genomic_start": [100, 200],
            "genomic_end": [200, 300],
            "te_id": ["TE1", "TE2"],
            "te_start": [1, 1],
            "te_end": [50, 50],
        })
        dedup = HitDeduplicator()
        result = dedup.deduplicate(df)
        assert len(result) == 2
        assert dedup.stats["duplicates_removed"] == 0


class TestTECoordinateNormalization:
    def test_normalizes_te_coords(self):
        """TE start/end should be normalized so start < end before keying."""
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 100],   # forward
            "te_end": [100, 1],     # reverse (same region, different strand)
            "pident": [80.0, 80.0],
        })
        dedup = HitDeduplicator(normalize_te_coords=True)
        result = dedup.deduplicate(df)
        # After normalization, both rows have te_start=1, te_end=100 → duplicate
        assert len(result) == 1

    def test_without_normalization_keeps_both(self):
        """Without normalization, (1,100) and (100,1) are different keys."""
        df = pd.DataFrame({
            "gene_id": ["g1", "g1"],
            "chrom": ["chr1", "chr1"],
            "genomic_start": [100, 100],
            "genomic_end": [200, 200],
            "te_id": ["TE1", "TE1"],
            "te_start": [1, 100],
            "te_end": [100, 1],
            "pident": [80.0, 80.0],
        })
        dedup = HitDeduplicator(normalize_te_coords=False)
        result = dedup.deduplicate(df)
        assert len(result) == 2


class TestDuplicationReport:
    def test_per_gene_stats(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        gene_stats = dedup.per_gene_stats(annotated_hits)
        assert "g1" in gene_stats
        assert gene_stats["g1"]["raw_hits"] == 3
        assert gene_stats["g1"]["unique_hits"] == 2
        assert gene_stats["g1"]["duplicates_removed"] == 1

    def test_per_gene_stats_empty(self):
        dedup = HitDeduplicator()
        df = pd.DataFrame(columns=DEFAULT_KEY)
        dedup.deduplicate(df)
        gene_stats = dedup.per_gene_stats(df)
        assert gene_stats == {}

    def test_duplication_rate(self, annotated_hits):
        dedup = HitDeduplicator()
        dedup.deduplicate(annotated_hits)
        assert dedup.stats["duplication_rate"] == pytest.approx(2 / 5)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_blast/test_dedup.py -v`
Expected: FAIL — `fossil_finder.blast.dedup` does not exist

- [ ] **Step 3: Implement HitDeduplicator**

Create `src/fossil_finder/blast/dedup.py`:

```python
"""Coordinate-based BLAST hit deduplication.

Removes duplicate hits that arise when multiple transcripts of the same
gene produce identical BLAST alignments against the same TE region.
Replaces the legacy scripts/deduplicate_*_te_hits.py.
"""

import pandas as pd


DEFAULT_KEY_COLUMNS = [
    "gene_id", "chrom", "genomic_start", "genomic_end",
    "te_id", "te_start", "te_end",
]


class HitDeduplicator:
    """Remove duplicate BLAST hits by coordinate key.

    Args:
        key_columns: Columns forming the dedup key. Defaults to the
            legacy 7-tuple (gene_id, chrom, genomic_start, genomic_end,
            te_id, te_start, te_end).
        normalize_te_coords: If True, normalize te_start/te_end so
            start < end before computing the key (handles forward/reverse
            hits to the same TE region).
    """

    def __init__(
        self,
        key_columns: list[str] | None = None,
        normalize_te_coords: bool = True,
    ):
        self.key_columns = key_columns or DEFAULT_KEY_COLUMNS
        self.normalize_te_coords = normalize_te_coords
        self.stats: dict = {
            "total_input": 0,
            "unique_hits": 0,
            "duplicates_removed": 0,
            "duplication_rate": 0.0,
        }

    def deduplicate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove duplicate hits from annotated BLAST results.

        Args:
            df: DataFrame with columns matching key_columns.

        Returns:
            Deduplicated DataFrame (first occurrence kept).
        """
        if df.empty:
            self.stats = {
                "total_input": 0,
                "unique_hits": 0,
                "duplicates_removed": 0,
                "duplication_rate": 0.0,
            }
            return df.copy()

        work = df.copy()

        if self.normalize_te_coords and "te_start" in work.columns and "te_end" in work.columns:
            te_min = work[["te_start", "te_end"]].min(axis=1)
            te_max = work[["te_start", "te_end"]].max(axis=1)
            work["te_start"] = te_min
            work["te_end"] = te_max

        result = work.drop_duplicates(subset=self.key_columns, keep="first")

        total = len(df)
        unique = len(result)
        self.stats = {
            "total_input": total,
            "unique_hits": unique,
            "duplicates_removed": total - unique,
            "duplication_rate": (total - unique) / total if total > 0 else 0.0,
        }

        return result.reset_index(drop=True)

    def per_gene_stats(
        self, df: pd.DataFrame, gene_col: str = "gene_id",
    ) -> dict[str, dict]:
        """Compute per-gene duplication statistics.

        Args:
            df: The original (pre-dedup) DataFrame.
            gene_col: Column containing gene identifiers.

        Returns:
            Dict mapping gene_id -> {raw_hits, unique_hits, duplicates_removed}.
        """
        if df.empty or gene_col not in df.columns:
            return {}

        result = {}
        for gene_id, group in df.groupby(gene_col):
            # Use a temporary deduplicator to avoid corrupting self.stats
            temp = HitDeduplicator(
                key_columns=self.key_columns,
                normalize_te_coords=self.normalize_te_coords,
            )
            deduped = temp.deduplicate(group)
            raw = len(group)
            unique = len(deduped)
            result[gene_id] = {
                "raw_hits": raw,
                "unique_hits": unique,
                "duplicates_removed": raw - unique,
            }

        return result
```

- [ ] **Step 4: Update `blast/__init__.py` exports**

```python
"""BLAST search execution and hit processing."""

from .dedup import HitDeduplicator
from .filter import apply_filters, filter_by_evalue, filter_by_length, filter_by_pident
from .runner import BlastRunner

__all__ = [
    "BlastRunner",
    "HitDeduplicator",
    "apply_filters",
    "filter_by_evalue",
    "filter_by_length",
    "filter_by_pident",
]
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pytest tests/test_blast/test_dedup.py -v`
Expected: ALL PASS

- [ ] **Step 6: Run full suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL tests PASS (~190 total: 155 existing + ~35 new)

- [ ] **Step 7: Commit**

```bash
git add src/fossil_finder/blast/dedup.py src/fossil_finder/blast/__init__.py \
  tests/test_blast/test_dedup.py
git commit -m "feat: HitDeduplicator — coordinate-based BLAST hit dedup"
```

---

### Task 6: Version bump + final integration check

**Files:**
- Modify: `src/fossil_finder/__init__.py`

- [ ] **Step 1: Bump version**

Update `src/fossil_finder/__init__.py`:

```python
"""Fossil Finder: Multi-genome TE fossil mining and regulatory analysis framework."""

__version__ = "0.4.0"
```

- [ ] **Step 2: Run full test suite**

Run: `pytest tests/ -v --tb=short`
Expected: ALL PASS

- [ ] **Step 3: Verify imports work cleanly**

Run: `python -c "from fossil_finder.blast import BlastRunner, HitDeduplicator, apply_filters; print('Phase 4 imports OK')"`
Expected: `Phase 4 imports OK`

- [ ] **Step 4: Commit**

```bash
git add src/fossil_finder/__init__.py
git commit -m "chore: bump version to 0.4.0 for Phase 4"
```

---

## Summary

| Task | Module | Tests | What it replaces |
|------|--------|-------|------------------|
| 1 | `config/schema.py` (extend BlastSpec) | 5 | Legacy `blast_runner.py` defaults |
| 2 | `blast/runner.py` (BlastRunner) | 12 | `scripts/blast_runner.py` (375 lines) |
| 3 | Integration test (optional) | 1 | Manual BLAST verification |
| 4 | `blast/filter.py` (HitFilter) | 10 | Ad-hoc filtering in ~6 analysis scripts |
| 5 | `blast/dedup.py` (HitDeduplicator) | 11 | `scripts/deduplicate_*_te_hits.py` (~800 lines combined) |
| 6 | Version bump | 0 | — |

**Total new tests:** ~39
**Estimated total after Phase 4:** ~194 tests
