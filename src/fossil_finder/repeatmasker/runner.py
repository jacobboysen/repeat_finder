"""Config-driven RepeatMasker command execution.

Builds RepeatMasker command-lines from RepeatMaskerSpec, executes them
via subprocess, and returns the path to the .out file for downstream
overlap analysis.
"""

import logging
import shutil
import subprocess
from pathlib import Path

from fossil_finder.config.schema import RepeatMaskerSpec

logger = logging.getLogger(__name__)


class RepeatMaskerRunner:
    """Execute RepeatMasker using RepeatMaskerSpec configuration.

    Wraps subprocess RepeatMasker execution with:
    - Command construction from RepeatMaskerSpec parameters
    - Output file discovery (RepeatMasker names outputs after the input)
    - Validation of inputs and binary availability

    Usage:
        runner = RepeatMaskerRunner(spec)
        out_path = runner.run(input_fasta, output_dir)
        # out_path is the .out file, ready for parse_repeatmasker_out()
    """

    def __init__(self, spec: RepeatMaskerSpec):
        self.spec = spec

    def build_command(
        self,
        input_fasta: str | Path,
        output_dir: str | Path,
    ) -> list[str]:
        """Build RepeatMasker command-line from spec.

        Args:
            input_fasta: Path to input FASTA file (must exist).
            output_dir: Directory for RepeatMasker output files.

        Returns:
            Command as list of strings for subprocess.run().

        Raises:
            FileNotFoundError: If input FASTA does not exist.
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)

        if not input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

        s = self.spec
        cmd = [
            "RepeatMasker",
            "-dir", str(output_dir),
            "-e", s.engine,
            "-pa", str(s.parallel),
        ]

        # Library source: custom lib overrides species
        if s.lib:
            cmd.extend(["-lib", s.lib])
        else:
            cmd.extend(["-species", s.species])

        # Sensitivity flags
        if s.sensitivity == "slow":
            cmd.append("-s")
        elif s.sensitivity == "quick":
            cmd.append("-q")
        elif s.sensitivity == "rush":
            cmd.append("-qq")
        # "default" = no flag

        if s.gff:
            cmd.append("-gff")

        if s.no_is:
            cmd.append("-no_is")

        if s.xsmall:
            cmd.append("-xsmall")

        # Input file last
        cmd.append(str(input_fasta))

        return cmd

    def find_output(
        self,
        input_fasta: str | Path,
        output_dir: str | Path,
    ) -> Path:
        """Locate the .out file RepeatMasker produces.

        RepeatMasker names output files as ``<input_basename>.out``
        inside the output directory.

        Args:
            input_fasta: Path to the input FASTA used in the run.
            output_dir: Directory where RM wrote its outputs.

        Returns:
            Path to the .out file.

        Raises:
            FileNotFoundError: If the expected .out file doesn't exist.
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        expected = output_dir / f"{input_fasta.name}.out"

        if not expected.exists():
            raise FileNotFoundError(
                f"RepeatMasker output not found: {expected}. "
                f"Check that RepeatMasker completed successfully."
            )

        return expected

    def is_available(self) -> bool:
        """Check if RepeatMasker is installed and on PATH."""
        return shutil.which("RepeatMasker") is not None

    def run(
        self,
        input_fasta: str | Path,
        output_dir: str | Path,
    ) -> Path:
        """Execute RepeatMasker and return path to .out file.

        Args:
            input_fasta: Path to input FASTA.
            output_dir: Directory for output files.

        Returns:
            Path to the RepeatMasker .out file.

        Raises:
            FileNotFoundError: If input FASTA or RepeatMasker binary not found.
            subprocess.CalledProcessError: If RepeatMasker exits non-zero.
            RuntimeError: If RepeatMasker is not installed.
        """
        if not self.is_available():
            raise RuntimeError(
                "RepeatMasker not found on PATH. "
                "Install RepeatMasker or add it to your PATH."
            )

        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = self.build_command(input_fasta, output_dir)

        logger.info("Running RepeatMasker on %s", input_fasta.name)
        logger.debug("Full command: %s", " ".join(cmd))

        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
        )

        if result.stderr:
            logger.warning("RepeatMasker stderr: %s", result.stderr.strip())

        return self.find_output(input_fasta, output_dir)
