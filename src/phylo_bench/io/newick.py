"""Biopython-backed Newick parsing."""

from __future__ import annotations

from io import StringIO
import re
from math import isfinite

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


def parse_newick(newick: str) -> RootedTree:
    """Parse one Newick tree into the PhyloBench tree model.

    Parameters
    ----------
    newick:
        A single Newick tree, including its terminating semicolon.

    Raises
    ------
    ValueError
        If the input is empty, contains multiple trees, or fails validation.
    """

    stripped = newick.strip()
    if not stripped:
        raise ValueError("Newick input must not be empty")
    if not stripped.endswith(";"):
        raise ValueError("Newick input must end with exactly one semicolon")
    if ";" in stripped[:-1]:
        raise ValueError("Newick input must contain exactly one terminal semicolon")
    _validate_branch_lengths(stripped[:-1])

    try:
        trees = list(Phylo.parse(StringIO(newick), "newick"))
    except Exception as exc:
        raise ValueError(f"invalid Newick input: {exc}") from exc

    if len(trees) != 1:
        raise ValueError("Newick input must contain exactly one tree")

    next_id = 0
    root_clade = trees[0].root
    root = TreeNode(next_id, root_clade.name, root_clade.branch_length)
    next_id += 1
    stack = [(root_clade, root)]
    while stack:
        clade, parent = stack.pop()
        converted_children = []
        for child_clade in clade.clades:
            child = TreeNode(next_id, child_clade.name, child_clade.branch_length)
            next_id += 1
            parent.add_child(child)
            converted_children.append((child_clade, child))
        stack.extend(reversed(converted_children))

    tree = RootedTree(root)
    validate_tree(tree)
    return tree


_BRANCH_LENGTH = re.compile(r":([^,();\[\]\s]*)")


def _validate_branch_lengths(newick: str) -> None:
    """Reject malformed, negative, or non-finite Newick branch lengths."""

    for token in _BRANCH_LENGTH.findall(newick):
        if not token:
            raise ValueError("malformed branch length: empty value")
        try:
            value = float(token)
        except ValueError as exc:
            raise ValueError(f"malformed branch length: {token!r}") from exc
        if not isfinite(value) or value < 0:
            raise ValueError(f"invalid branch length: {token!r}")
