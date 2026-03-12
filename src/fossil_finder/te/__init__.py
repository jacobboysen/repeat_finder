"""TE (Transposable Element) classification and domain analysis."""

from .classifier import (
    classify_te_domain,
    get_domain_description,
    is_coding_domain,
    is_regulatory_domain,
)
from .taxonomy import (
    DEFAULT_FAMILY_LISTS,
    FLYBASE_OVERRIDES,
    infer_te_class,
    parse_consensus_fasta,
    parse_instance_fasta,
    strip_instance_suffix,
)

__all__ = [
    "classify_te_domain",
    "get_domain_description",
    "infer_te_class",
    "is_coding_domain",
    "is_regulatory_domain",
    "DEFAULT_FAMILY_LISTS",
    "FLYBASE_OVERRIDES",
    "parse_consensus_fasta",
    "parse_instance_fasta",
    "strip_instance_suffix",
]
