"""Tests for the Snakemake workflow and bridge script."""

import ast
import shutil
import subprocess
from pathlib import Path

import pytest

WORKFLOWS_DIR = Path(__file__).resolve().parents[2] / "workflows"
SNAKEFILE = WORKFLOWS_DIR / "Snakefile"
BRIDGE_SCRIPT = WORKFLOWS_DIR / "scripts" / "run_analysis.py"


def snakemake_available() -> bool:
    return shutil.which("snakemake") is not None


def blast_available() -> bool:
    return shutil.which("blastn") is not None


def _load_bridge_module():
    """Import the bridge script as a module (main() won't run on import)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "run_analysis", str(BRIDGE_SCRIPT),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestBridgeScript:
    def test_bridge_script_exists(self):
        assert BRIDGE_SCRIPT.exists()

    def test_bridge_script_valid_python(self):
        source = BRIDGE_SCRIPT.read_text()
        ast.parse(source)  # Raises SyntaxError if invalid

    def test_build_query_to_gene_from_metadata(self, test_data_dir, tmp_path):
        """query_to_gene resolves mRNA parent_id -> gene_id via GFF3."""
        mod = _load_bridge_module()

        # Write a mini regions.tsv with mRNA IDs as parent_id (GFF3 reality)
        regions_tsv = tmp_path / "regions.tsv"
        regions_tsv.write_text(
            "region_id\tchrom\tstart\tend\tstrand\tparent_id\tfeature_type\tlength\n"
            "utr3_001\tchr1\t281\t300\t+\tmRNA001\tthree_prime_UTR\t20\n"
            "utr3_002\tchr1\t350\t369\t-\tmRNA002\tthree_prime_UTR\t20\n"
        )

        annotation_gff = test_data_dir / "mini_annotation.gff3"
        q2g = mod.build_query_to_gene(regions_tsv, annotation_gff)
        # Should resolve mRNA001 -> gene001, mRNA002 -> gene002
        assert q2g["utr3_001"] == "gene001"
        assert q2g["utr3_002"] == "gene002"

    def test_build_query_to_gene_unknown_mrna_falls_back(self, test_data_dir, tmp_path):
        """Unknown mRNA IDs fall back to parent_id as-is."""
        mod = _load_bridge_module()

        regions_tsv = tmp_path / "regions.tsv"
        regions_tsv.write_text(
            "region_id\tchrom\tstart\tend\tstrand\tparent_id\tfeature_type\tlength\n"
            "utr_x\tchr1\t100\t200\t+\tunknown_mrna\tthree_prime_UTR\t100\n"
        )

        annotation_gff = test_data_dir / "mini_annotation.gff3"
        q2g = mod.build_query_to_gene(regions_tsv, annotation_gff)
        assert q2g["utr_x"] == "unknown_mrna"

    def test_load_gene_sets_from_config(self, tmp_path):
        """Gene set files are loaded into dict[str, set[str]]."""
        mod = _load_bridge_module()

        # Write gene list files
        gs_dir = tmp_path / "gene_lists"
        gs_dir.mkdir()
        (gs_dir / "set_a.txt").write_text("gene1\ngene2\ngene3\n")
        (gs_dir / "set_b.txt").write_text("gene2\ngene4\n")

        # Mock GeneSetSpec objects
        from fossil_finder.config.schema import GeneSetSpec
        gene_set_specs = {
            "set_a": GeneSetSpec(genes=str(gs_dir / "set_a.txt")),
            "set_b": GeneSetSpec(genes=str(gs_dir / "set_b.txt")),
        }

        result = mod.load_gene_sets(gene_set_specs, base_dir=None)
        assert result["set_a"] == {"gene1", "gene2", "gene3"}
        assert result["set_b"] == {"gene2", "gene4"}


@pytest.mark.skipif(not snakemake_available(), reason="snakemake not installed")
class TestSnakemakeDryRun:
    def test_snakefile_parses(self, test_data_dir, tmp_path):
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
        for rule_name in ["extract_regions", "make_blast_db", "run_blast",
                          "run_repeatmasker", "analyze"]:
            assert any(rule_name in r for r in rules), f"Missing rule: {rule_name}"

    def test_config_validation_fails_on_missing(self, tmp_path):
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
        assert result.returncode != 0


@pytest.mark.skipif(
    not (snakemake_available() and blast_available()),
    reason="snakemake + BLAST+ required",
)
class TestSnakemakeIntegration:
    def test_full_pipeline_mini_genome(self, test_data_dir, tmp_path):
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

        out = tmp_path / "output" / "test_v1"
        assert (out / "regions.fa").exists()
        assert (out / "regions.tsv").exists()
        assert (out / "blast_results.tsv").exists()
        assert (out / "summary.json").exists()
