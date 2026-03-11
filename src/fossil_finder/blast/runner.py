"""Config-driven BLAST command execution.

Builds BLAST command-lines from BlastSpec, executes them via subprocess,
and returns parsed results using the existing io/blast.py layer.
"""

import logging
import subprocess
from pathlib import Path

import pandas as pd

from fossil_finder.config.schema import BlastSpec
from fossil_finder.io.blast import BLAST_COLUMNS_NO_STRAND, load_blast_results

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
