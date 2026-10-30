"""Lowest Common Ancestor preprocessing and queries."""

from __future__ import annotations

from math import ceil, log2

from phylo_bench.algorithms.traversal import preorder
from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


class LCAIndex:
    """Binary-lifting index for LCA queries on a rooted tree.

    Preprocessing takes O(n log n) time and space. Each query takes O(log n)
    time, where n is the number of nodes.
    """

    def __init__(self, tree: RootedTree) -> None:
        validate_tree(tree)
        tree.freeze()
        self._tree = tree
        nodes = list(preorder(tree))
        self._depth = {node.node_id: 0 for node in nodes}
        self._nodes = {node.node_id: node for node in nodes}
        for node in nodes:
            if node.parent is not None:
                self._depth[node.node_id] = self._depth[node.parent.node_id] + 1

        levels = max(1, ceil(log2(max(1, len(nodes)))))
        self._up: list[dict[int, int]] = []
        self._up.append({
            node.node_id: node.parent.node_id if node.parent else node.node_id
            for node in nodes
        })
        for level in range(1, levels):
            previous = self._up[-1]
            self._up.append({node_id: previous[previous[node_id]] for node_id in previous})

    def lca(self, first: TreeNode, second: TreeNode) -> TreeNode:
        """Return the lowest common ancestor of two indexed nodes."""

        self._ensure_indexed(first)
        self._ensure_indexed(second)
        a, b = first.node_id, second.node_id

        if self._depth[a] < self._depth[b]:
            a, b = b, a
        a = self._raise(a, self._depth[a] - self._depth[b])
        if a == b:
            return self._nodes[a]

        for level in range(len(self._up) - 1, -1, -1):
            if self._up[level][a] != self._up[level][b]:
                a = self._up[level][a]
                b = self._up[level][b]
        return self._nodes[self._up[0][a]]

    def _raise(self, node_id: int, distance: int) -> int:
        level = 0
        while distance:
            if distance & 1:
                node_id = self._up[level][node_id]
            distance >>= 1
            level += 1
        return node_id

    def _ensure_indexed(self, node: TreeNode) -> None:
        if node.node_id not in self._nodes or self._nodes[node.node_id] is not node:
            raise ValueError("both nodes must belong to the indexed tree")
