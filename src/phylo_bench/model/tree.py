"""Typed rooted-tree data structures.

The model deliberately contains no reconciliation-specific behavior.  Node
identity is represented by a stable integer ID, while labels remain optional
metadata supplied by an input format such as Newick.
"""

from __future__ import annotations

from typing import Iterator


class TreeNode:
    """A node in a rooted tree."""

    __slots__ = (
        "_node_id", "_name", "_branch_length", "_children", "_parent", "_owner"
    )

    def __init__(
        self,
        node_id: int,
        name: str | None = None,
        branch_length: float | None = None,
        parent: "TreeNode | None" = None,
    ) -> None:
        self._node_id = node_id
        self._name = name
        self._branch_length = branch_length
        self._children: list[TreeNode] = []
        self._parent = parent
        self._owner: RootedTree | None = None

    @property
    def node_id(self) -> int:
        return self._node_id

    @property
    def name(self) -> str | None:
        return self._name

    @property
    def branch_length(self) -> float | None:
        return self._branch_length

    @property
    def children(self) -> tuple["TreeNode", ...]:
        return tuple(self._children)

    @property
    def parent(self) -> "TreeNode | None":
        return self._parent

    @property
    def is_leaf(self) -> bool:
        """Return whether this node has no children."""

        return not self.children

    def add_child(self, child: "TreeNode") -> None:
        """Attach ``child`` and maintain the parent relationship."""

        owner = getattr(self, "_owner", None)
        if owner is not None and owner._frozen:
            raise RuntimeError("tree topology is frozen")
        if child is self or self._is_descendant_of(child):
            raise ValueError("attaching the child would create a cycle")
        if child in self._children:
            raise ValueError("child is already attached to this parent")
        if child._parent is not None and child._parent is not self:
            raise ValueError("a node cannot have multiple parents")
        if child._owner is not None and child._owner is not owner:
            raise ValueError("a node cannot belong to multiple trees")
        child._parent = self
        self._children.append(child)

    def _is_descendant_of(self, candidate: "TreeNode") -> bool:
        current: TreeNode | None = self
        while current is not None:
            if current is candidate:
                return True
            current = current.parent
        return False


class RootedTree:
    """A rooted tree with a stable node index.

    Polytomies are allowed by default; callers requiring a binary tree should
    pass ``require_binary=True`` to :func:`validate_tree`. A named singleton
    leaf is valid.
    """

    __slots__ = ("_root", "_frozen")

    def __init__(self, root: TreeNode) -> None:
        self._root = root
        self._frozen = False
        self._register_owner()

    @property
    def root(self) -> TreeNode:
        return self._root

    def freeze(self) -> None:
        """Prevent topology mutations after an algorithm indexes the tree."""

        self._frozen = True
        for node in self.nodes():
            node._owner = self

    def _register_owner(self) -> None:
        seen: set[int] = set()
        for node in self.nodes():
            if id(node) in seen:
                raise ValueError("tree contains a cycle or shared node")
            seen.add(id(node))
            if node._owner is not None and node._owner is not self:
                raise ValueError("a node cannot belong to multiple trees")
            node._owner = self

    def nodes(self) -> Iterator[TreeNode]:
        """Yield all nodes in deterministic pre-order."""

        stack = [self.root]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def node_by_id(self, node_id: int) -> TreeNode:
        """Return the node with ``node_id`` or raise ``KeyError``."""

        for node in self.nodes():
            if node.node_id == node_id:
                return node
        raise KeyError(node_id)
