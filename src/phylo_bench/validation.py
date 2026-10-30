"""Validation for rooted tree inputs."""

from __future__ import annotations

from math import isfinite

from phylo_bench.model.tree import RootedTree, TreeNode


def validate_tree(
    tree: RootedTree,
    *,
    require_named_leaves: bool = True,
    require_binary: bool = False,
) -> None:
    """Validate structural invariants of a rooted tree.

    Internal node names may be absent or repeated because algorithms use stable
    node IDs. Leaf names are required by default for biological taxon mapping.
    Polytomies are allowed unless ``require_binary`` is true. A named singleton
    leaf is valid.
    """

    if not isinstance(tree.root, TreeNode):
        raise ValueError("tree root must be a TreeNode")
    if tree.root.parent is not None:
        raise ValueError("the root must not have a parent")

    seen_identity: set[int] = set()
    seen_ids: set[int] = set()
    seen_leaf_names: set[str] = set()
    stack: list[tuple[TreeNode, TreeNode | None]] = [(tree.root, None)]

    while stack:
        node, expected_parent = stack.pop()
        identity = id(node)
        if identity in seen_identity:
            raise ValueError("tree contains a cycle or shared node")
        seen_identity.add(identity)

        if isinstance(node.node_id, bool) or not isinstance(node.node_id, int):
            raise ValueError(f"node {node.node_id!r} has an invalid ID")
        if node.node_id < 0:
            raise ValueError(f"node {node.node_id} has a negative ID")
        if node.node_id in seen_ids:
            raise ValueError(f"duplicate node ID: {node.node_id}")
        seen_ids.add(node.node_id)

        if node.parent is not expected_parent:
            raise ValueError(f"incorrect parent pointer for node {node.node_id}")
        if node.branch_length is not None and (
            not isinstance(node.branch_length, (int, float))
            or isinstance(node.branch_length, bool)
            or not isfinite(float(node.branch_length))
            or node.branch_length < 0
        ):
            raise ValueError(f"node {node.node_id} has an invalid branch length")

        if node.is_leaf:
            if require_named_leaves and not node.name:
                raise ValueError(f"leaf {node.node_id} has no name")
            if node.name is not None:
                if node.name in seen_leaf_names:
                    raise ValueError(f"duplicate leaf name: {node.name}")
                seen_leaf_names.add(node.name)
        else:
            if require_binary and len(node.children) != 2:
                raise ValueError(
                    f"internal node {node.node_id} is not binary"
                )
            for child in reversed(node.children):
                stack.append((child, node))
