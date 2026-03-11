"""Tests for BlastRunner."""

import shutil
import subprocess

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

    def test_outfmt_includes_16_columns(self, tmp_path):
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
