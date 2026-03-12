"""Tests for RepeatMasker runner — command construction, output discovery, validation."""

from pathlib import Path
from unittest.mock import patch

import pytest

from fossil_finder.config.schema import RepeatMaskerSpec
from fossil_finder.repeatmasker.runner import RepeatMaskerRunner


@pytest.fixture
def default_spec():
    return RepeatMaskerSpec()


@pytest.fixture
def runner(default_spec):
    return RepeatMaskerRunner(default_spec)


@pytest.fixture
def fake_fasta(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1\nACGTACGT\n")
    return fasta


class TestBuildCommand:
    def test_default_command(self, runner, fake_fasta, tmp_path):
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert cmd[0] == "RepeatMasker"
        assert "-dir" in cmd
        assert str(tmp_path / "out") in cmd
        assert "-e" in cmd
        assert "rmblast" in cmd
        assert "-pa" in cmd
        assert "4" in cmd
        assert "-species" in cmd
        assert "drosophila" in cmd
        assert str(fake_fasta) == cmd[-1]

    def test_gff_flag_default_on(self, runner, fake_fasta, tmp_path):
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-gff" in cmd

    def test_gff_flag_off(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(gff=False)
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-gff" not in cmd

    def test_custom_lib_overrides_species(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(lib="/path/to/custom.lib")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-lib" in cmd
        assert "/path/to/custom.lib" in cmd
        assert "-species" not in cmd

    def test_slow_sensitivity(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(sensitivity="slow")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-s" in cmd

    def test_quick_sensitivity(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(sensitivity="quick")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-q" in cmd

    def test_rush_sensitivity(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(sensitivity="rush")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-qq" in cmd

    def test_default_sensitivity_no_flag(self, runner, fake_fasta, tmp_path):
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-s" not in cmd
        assert "-q" not in cmd
        assert "-qq" not in cmd

    def test_no_is_flag(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(no_is=True)
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-no_is" in cmd

    def test_xsmall_flag(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(xsmall=True)
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert "-xsmall" in cmd

    def test_custom_parallel(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(parallel=8)
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        idx = cmd.index("-pa")
        assert cmd[idx + 1] == "8"

    def test_custom_species(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(species="human")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        idx = cmd.index("-species")
        assert cmd[idx + 1] == "human"

    def test_custom_engine(self, fake_fasta, tmp_path):
        spec = RepeatMaskerSpec(engine="crossmatch")
        runner = RepeatMaskerRunner(spec)
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        idx = cmd.index("-e")
        assert cmd[idx + 1] == "crossmatch"

    def test_nonexistent_fasta_raises(self, runner, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            runner.build_command(tmp_path / "nonexistent.fa", tmp_path / "out")

    def test_input_fasta_is_last_arg(self, runner, fake_fasta, tmp_path):
        cmd = runner.build_command(fake_fasta, tmp_path / "out")
        assert cmd[-1] == str(fake_fasta)


class TestFindOutput:
    def test_finds_existing_out_file(self, runner, tmp_path):
        fasta = tmp_path / "genome.fasta"
        fasta.touch()
        out_dir = tmp_path / "rm_out"
        out_dir.mkdir()
        (out_dir / "genome.fasta.out").write_text("dummy")
        result = runner.find_output(fasta, out_dir)
        assert result == out_dir / "genome.fasta.out"

    def test_raises_if_out_missing(self, runner, tmp_path):
        fasta = tmp_path / "genome.fasta"
        fasta.touch()
        out_dir = tmp_path / "rm_out"
        out_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="output not found"):
            runner.find_output(fasta, out_dir)


class TestIsAvailable:
    def test_available_when_on_path(self, runner):
        with patch("fossil_finder.repeatmasker.runner.shutil.which", return_value="/usr/local/bin/RepeatMasker"):
            assert runner.is_available() is True

    def test_not_available_when_missing(self, runner):
        with patch("fossil_finder.repeatmasker.runner.shutil.which", return_value=None):
            assert runner.is_available() is False


class TestRun:
    def test_raises_if_not_installed(self, runner, fake_fasta, tmp_path):
        with patch.object(runner, "is_available", return_value=False):
            with pytest.raises(RuntimeError, match="not found on PATH"):
                runner.run(fake_fasta, tmp_path / "out")

    def test_creates_output_dir(self, runner, fake_fasta, tmp_path):
        out_dir = tmp_path / "deep" / "nested" / "rm_out"
        with patch.object(runner, "is_available", return_value=True), \
             patch("fossil_finder.repeatmasker.runner.subprocess.run") as mock_run, \
             patch.object(runner, "find_output", return_value=out_dir / "genome.fasta.out"):
            mock_run.return_value.stderr = ""
            runner.run(fake_fasta, out_dir)
            assert out_dir.exists()

    def test_calls_subprocess_with_check(self, runner, fake_fasta, tmp_path):
        out_dir = tmp_path / "rm_out"
        with patch.object(runner, "is_available", return_value=True), \
             patch("fossil_finder.repeatmasker.runner.subprocess.run") as mock_run, \
             patch.object(runner, "find_output", return_value=out_dir / "genome.fasta.out"):
            mock_run.return_value.stderr = ""
            runner.run(fake_fasta, out_dir)
            mock_run.assert_called_once()
            call_kwargs = mock_run.call_args
            assert call_kwargs.kwargs["check"] is True
            assert call_kwargs.kwargs["capture_output"] is True

    def test_returns_out_path(self, runner, fake_fasta, tmp_path):
        out_dir = tmp_path / "rm_out"
        expected = out_dir / "genome.fasta.out"
        with patch.object(runner, "is_available", return_value=True), \
             patch("fossil_finder.repeatmasker.runner.subprocess.run") as mock_run, \
             patch.object(runner, "find_output", return_value=expected):
            mock_run.return_value.stderr = ""
            result = runner.run(fake_fasta, out_dir)
            assert result == expected


class TestSpecDefaults:
    def test_default_species(self):
        spec = RepeatMaskerSpec()
        assert spec.species == "drosophila"

    def test_default_engine(self):
        spec = RepeatMaskerSpec()
        assert spec.engine == "rmblast"

    def test_default_parallel(self):
        spec = RepeatMaskerSpec()
        assert spec.parallel == 4

    def test_default_sensitivity(self):
        spec = RepeatMaskerSpec()
        assert spec.sensitivity == "default"

    def test_default_lib_none(self):
        spec = RepeatMaskerSpec()
        assert spec.lib is None

    def test_default_gff_true(self):
        spec = RepeatMaskerSpec()
        assert spec.gff is True

    def test_invalid_sensitivity_raises(self):
        from pydantic import ValidationError
        with pytest.raises(ValidationError):
            RepeatMaskerSpec(sensitivity="invalid")

    def test_spec_in_genome_config(self):
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": "/path/genome.fa",
                "annotation_gff": "/path/ann.gff",
            },
            repeatmasker={"species": "human", "parallel": 8},
        )
        assert config.repeatmasker.species == "human"
        assert config.repeatmasker.parallel == 8

    def test_default_spec_in_genome_config(self):
        from fossil_finder.config.schema import GenomeConfig
        config = GenomeConfig(
            genome={"species": "Test", "assembly": "v1"},
            source={
                "adapter": "custom",
                "genome_fasta": "/path/genome.fa",
                "annotation_gff": "/path/ann.gff",
            },
        )
        assert config.repeatmasker.species == "drosophila"
        assert config.repeatmasker.engine == "rmblast"
