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
