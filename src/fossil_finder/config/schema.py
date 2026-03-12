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
    chromosomes: list[str] | None = None
    ucsc_prefix: str | None = None


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
    """Data source configuration."""

    adapter: Literal["flybase", "ensembl", "ncbi", "custom"]
    genome_fasta: str
    annotation_gff: str
    annotation_gff_full: str | None = None  # unfiltered GFF for large genomes
    gene_id_prefix: str | None = None
    transcript_id_prefix: str | None = None
    te_instances: str | None = None
    te_consensus: str | None = None
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


class RepeatMaskerSpec(BaseModel):
    """RepeatMasker execution parameters.

    Controls how RepeatMasker is invoked on input sequences.
    The species flag selects the built-in RM library for that organism.
    """

    species: str = "drosophila"
    engine: str = "rmblast"
    parallel: int = 4
    sensitivity: Literal["default", "slow", "quick", "rush"] = "default"
    lib: str | None = None  # custom RM library (overrides species)
    gff: bool = True  # also produce GFF output
    no_is: bool = False  # skip bacterial insertion element check
    xsmall: bool = False  # output soft-masked lowercase instead of N-masked


class ScoringSpec(BaseModel):
    """ML/LM scoring configuration."""

    model: str = "InstaDeepAI/NTv3_100M_post"
    species_key: str | None = None
    window_size: int = 32768


class GeneSetSpec(BaseModel):
    """A named gene set for comparative analysis."""

    description: str = ""
    genes: str
    tier: int = 1


class GenomeConfig(BaseModel):
    """Top-level genome configuration."""

    genome: GenomeSpec
    source: SourceSpec
    blast: BlastSpec = BlastSpec()
    repeatmasker: RepeatMaskerSpec = RepeatMaskerSpec()
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
    """Load and validate a genome config from YAML file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path) as f:
        raw = yaml.safe_load(f)

    return GenomeConfig(**raw)
