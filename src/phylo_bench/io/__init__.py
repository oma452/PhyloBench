"""Input and output adapters."""

from phylo_bench.io.atomic import atomic_write_bytes, atomic_write_text
from phylo_bench.io.newick import parse_newick

__all__ = ["atomic_write_bytes", "atomic_write_text", "parse_newick"]
