"""Classical LCA-based gene-tree/species-tree reconciliation.

The implementation assumes rooted binary trees, complete taxon sampling, a
fixed gene-leaf-to-species mapping, and only vertical speciation, duplication,
and loss events. The family-origin convention places the ancestral family at
the species-tree root. Losses are enumerated as branch-level events, so two
independent lineages lost on the same species branch produce two records.

For ``g`` gene-tree nodes and ``s`` species-tree nodes, reconciliation takes
O(g log s + g h) time and O(g + s log s + l) space, where ``h`` is the
species-tree height and ``l`` is the number of reported loss events. The
``g h`` term comes from enumerating each species-tree transition path.
"""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite

from phylo_bench.algorithms.lca import LCAIndex
from phylo_bench.algorithms.traversal import postorder
from phylo_bench.model.reconciliation import (
    EventType,
    LossEvent,
    ReconciliationEvent,
    ReconciliationResult,
)
from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


def map_gene_leaves(
    gene_tree: RootedTree,
    species_tree: RootedTree,
    leaf_to_species: Mapping[str, str] | None = None,
) -> dict[int, TreeNode]:
    """Map each gene leaf to a species-tree leaf.

    When ``leaf_to_species`` is omitted, gene leaf names must equal species
    leaf names. Otherwise it maps each gene leaf name to a species leaf name.
    """

    validate_tree(gene_tree, require_binary=True)
    validate_tree(species_tree, require_binary=True)
    species_by_name = {
        node.name: node
        for node in species_tree.nodes()
        if node.is_leaf and node.name is not None
    }
    mapping = leaf_to_species or {}
    result: dict[int, TreeNode] = {}
    for gene_leaf in (node for node in gene_tree.nodes() if node.is_leaf):
        if gene_leaf.name is None:
            raise ValueError(f"gene leaf {gene_leaf.node_id} has no name")
        if leaf_to_species is not None and gene_leaf.name not in leaf_to_species:
            raise ValueError(f"gene leaf {gene_leaf.name!r} has no species mapping")
        species_name = mapping.get(gene_leaf.name, gene_leaf.name)
        species_leaf = species_by_name.get(species_name)
        if species_leaf is None:
            raise ValueError(
                f"gene leaf {gene_leaf.name!r} maps to unknown species "
                f"{species_name!r}"
            )
        result[gene_leaf.node_id] = species_leaf
    return result


def reconcile(
    gene_tree: RootedTree,
    species_tree: RootedTree,
    leaf_to_species: Mapping[str, str] | None = None,
    *,
    duplication_cost: float = 1.0,
    loss_cost: float = 1.0,
) -> ReconciliationResult:
    """Reconcile a rooted binary gene tree against a rooted binary species tree.

    Internal gene nodes map to the species-tree LCA of their child mappings.
    A node is a duplication when either child maps to the same species node as
    its parent; otherwise it is a speciation. Losses use the approved
    transition convention: duplication transitions include the starting
    species split, while speciation transitions skip that already represented
    split. The family-origin transition starts at the species-tree root.
    """

    _validate_cost("duplication_cost", duplication_cost)
    _validate_cost("loss_cost", loss_cost)
    leaf_mapping = map_gene_leaves(gene_tree, species_tree, leaf_to_species)
    species_lca = LCAIndex(species_tree)
    node_mapping: dict[int, TreeNode] = dict(leaf_mapping)
    events: dict[int, ReconciliationEvent] = {
        node_id: ReconciliationEvent(
            gene_node_id=node_id,
            event_type=EventType.LEAF,
            species_node_id=species.node_id,
        )
        for node_id, species in leaf_mapping.items()
    }

    for gene_node in postorder(gene_tree):
        if gene_node.is_leaf:
            continue
        left, right = gene_node.children
        left_species = node_mapping[left.node_id]
        right_species = node_mapping[right.node_id]
        mapped = species_lca.lca(left_species, right_species)
        node_mapping[gene_node.node_id] = mapped
        event_type = (
            EventType.DUPLICATION
            if mapped is left_species or mapped is right_species
            else EventType.SPECIATION
        )
        events[gene_node.node_id] = ReconciliationEvent(
            gene_node_id=gene_node.node_id,
            event_type=event_type,
            species_node_id=mapped.node_id,
        )

    losses: list[LossEvent] = []
    root_mapping = node_mapping[gene_tree.root.node_id]
    losses.extend(
        _transition_losses(
            start=species_tree.root,
            target=root_mapping,
            include_start=True,
            gene_parent_id=None,
            gene_child_id=gene_tree.root.node_id,
            transition="origin",
        )
    )
    for parent in postorder(gene_tree):
        if parent.is_leaf:
            continue
        parent_event = events[parent.node_id].event_type
        include_start = parent_event is EventType.DUPLICATION
        for child in parent.children:
            losses.extend(
                _transition_losses(
                    start=node_mapping[parent.node_id],
                    target=node_mapping[child.node_id],
                    include_start=include_start,
                    gene_parent_id=parent.node_id,
                    gene_child_id=child.node_id,
                    transition=parent_event.value,
                )
            )

    duplication_count = sum(
        event.event_type is EventType.DUPLICATION for event in events.values()
    )
    loss_count = len(losses)
    total_cost = duplication_cost * duplication_count + loss_cost * loss_count
    return ReconciliationResult(
        leaf_mapping=leaf_mapping,
        node_mapping=node_mapping,
        events=events,
        losses=tuple(losses),
        duplication_count=duplication_count,
        loss_count=loss_count,
        duplication_cost=float(duplication_cost),
        loss_cost=float(loss_cost),
        total_cost=float(total_cost),
    )


def _transition_losses(
    *,
    start: TreeNode,
    target: TreeNode,
    include_start: bool,
    gene_parent_id: int | None,
    gene_child_id: int,
    transition: str,
) -> list[LossEvent]:
    path = _species_path(start, target)
    losses: list[LossEvent] = []
    first_split = True
    for current, selected in zip(path, path[1:]):
        if include_start or not first_split:
            for sibling in current.children:
                if sibling is not selected:
                    losses.append(
                        LossEvent(
                            species_parent_id=current.node_id,
                            species_child_id=sibling.node_id,
                            gene_parent_id=gene_parent_id,
                            gene_child_id=gene_child_id,
                            transition=transition,
                        )
                    )
        first_split = False
    return losses


def _species_path(start: TreeNode, target: TreeNode) -> list[TreeNode]:
    reverse_path: list[TreeNode] = []
    current: TreeNode | None = target
    while current is not None and current is not start:
        reverse_path.append(current)
        current = current.parent
    if current is None:
        raise ValueError("target mapping is not a descendant of the start mapping")
    reverse_path.append(start)
    return list(reversed(reverse_path))


def _validate_cost(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite non-negative number")
    if not isfinite(float(value)) or value < 0:
        raise ValueError(f"{name} must be a finite non-negative number")
