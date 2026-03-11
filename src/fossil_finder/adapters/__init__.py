"""Genome adapter layer -- pluggable organism-specific data loading."""

from .base import GenomeAdapter

__all__ = ["GenomeAdapter", "get_adapter"]


def get_adapter(config) -> GenomeAdapter:
    """Factory: return the appropriate adapter for a genome config.

    Args:
        config: A GenomeConfig instance (config.source.adapter selects the adapter).

    Returns:
        Concrete GenomeAdapter subclass instance.

    Raises:
        ValueError: If the adapter name is not recognized.
    """
    from .custom import CustomAdapter
    from .flybase import FlyBaseAdapter

    _registry: dict[str, type[GenomeAdapter]] = {
        "flybase": FlyBaseAdapter,
        "custom": CustomAdapter,
        # "ensembl": EnsemblAdapter,  # Phase 9+
        # "ncbi": NCBIAdapter,        # Phase 9+
    }

    adapter_name = config.source.adapter
    cls = _registry.get(adapter_name)
    if cls is None:
        raise ValueError(
            f"Unknown adapter '{adapter_name}'. "
            f"Available: {sorted(_registry.keys())}"
        )
    return cls(config)
