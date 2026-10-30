"""Iterative tree traversals."""

from __future__ import annotations

from collections.abc import Iterator

from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


def preorder(tree: RootedTree) -> Iterator[TreeNode]:
    """Yield nodes in root-first depth-first order."""

    validate_tree(tree, require_named_leaves=False)
    stack = [tree.root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(reversed(node.children))


def postorder(tree: RootedTree) -> Iterator[TreeNode]:
    """Yield nodes after all descendants, without recursion."""

    validate_tree(tree, require_named_leaves=False)
    stack: list[tuple[TreeNode, bool]] = [(tree.root, False)]
    while stack:
        node, expanded = stack.pop()
        if expanded:
            yield node
            continue
        stack.append((node, True))
        stack.extend((child, False) for child in reversed(node.children))
