"""TE (Transposable Element) classification and domain analysis."""

from .classifier import (
    classify_te_domain,
    get_domain_description,
    is_coding_domain,
    is_regulatory_domain,
)
from .taxonomy import DEFAULT_FAMILY_LISTS, infer_te_class

__all__ = [
    "classify_te_domain",
    "get_domain_description",
    "infer_te_class",
    "is_coding_domain",
    "is_regulatory_domain",
    "DEFAULT_FAMILY_LISTS",
]
