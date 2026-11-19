"""Exhaustive small-tree checks against an independent reconciliation oracle."""

from __future__ import annotations

from collections import Counter
from itertools import product

from phylo_bench.algorithms import reconcile
from phylo_bench.io import parse_newick
from phylo_bench.model import EventType


def _topologies(labels: tuple[str, ...]) -> tuple[str, ...]:
    if len(labels) == 2:
        return (f"({labels[0]},{labels[1]})",)
    first = labels[0]
    results: set[str] = set()
    remaining = labels[1:]
    for mask in range(1, 1 << len(remaining)):
        left = [first]
        right = []
        for index, label in enumerate(remaining):
            (left if mask & (1 << index) else right).append(label)
        if not right:
            continue
        left_trees = _topologies(tuple(left)) if len(left) > 1 else (left[0],)
        right_trees = _topologies(tuple(right)) if len(right) > 1 else (right[0],)
        for left_tree in left_trees:
            for right_tree in right_trees:
                children = tuple(sorted((left_tree, right_tree)))
                results.add(f"({children[0]},{children[1]})")
    return tuple(sorted(results))


def _postorder(node):
    for child in node.children:
        yield from _postorder(child)
    yield node


def _ancestor_lca(first, second):
    ancestors = set()
    current = first
    while current is not None:
        ancestors.add(current)
        current = current.parent
    current = second
    while current not in ancestors:
        current = current.parent
    return current


def _path(start, target):
    reverse = []
    current = target
    while current is not start:
        reverse.append(current)
        current = current.parent
    return [start, *reversed(reverse)]


def _oracle(gene_tree, species_tree, leaf_to_species):
    species_by_name = {
        node.name: node
        for node in species_tree.nodes()
        if node.is_leaf
    }
    mapping = {}
    events = {}
    for node in gene_tree.nodes():
        if node.is_leaf:
            mapping[node.node_id] = species_by_name[leaf_to_species[node.name]]
            events[node.node_id] = EventType.LEAF

    for node in _postorder(gene_tree.root):
        if node.is_leaf:
            continue
        children = node.children
        mapped = _ancestor_lca(
            mapping[children[0].node_id],
            mapping[children[1].node_id],
        )
        mapping[node.node_id] = mapped
        events[node.node_id] = (
            EventType.DUPLICATION
            if mapped in {
                mapping[children[0].node_id],
                mapping[children[1].node_id],
            }
            else EventType.SPECIATION
        )

    losses = []

    def transition(start, target, include_start):
        path = _path(start, target)
        for index, (current, selected) in enumerate(zip(path, path[1:])):
            if include_start or index:
                for sibling in current.children:
                    if sibling is not selected:
                        losses.append((current.node_id, sibling.node_id))

    root_mapping = mapping[gene_tree.root.node_id]
    transition(species_tree.root, root_mapping, True)
    for parent in _postorder(gene_tree.root):
        if parent.is_leaf:
            continue
        include_start = events[parent.node_id] is EventType.DUPLICATION
        for child in parent.children:
            transition(
                mapping[parent.node_id],
                mapping[child.node_id],
                include_start,
            )
    return mapping, events, Counter(losses)


def test_small_binary_reconciliations_match_independent_oracle():
    for species_count in (2, 3):
        species_names = tuple(chr(ord("A") + index) for index in range(species_count))
        species_topologies = _topologies(species_names)
        gene_names = tuple(f"G{index}" for index in range(3))
        for species_topology in species_topologies:
            species = parse_newick(f"{species_topology}S;")
            for gene_topology in _topologies(gene_names):
                gene = parse_newick(f"{gene_topology}G;")
                for assignment in product(species_names, repeat=len(gene_names)):
                    leaf_mapping = dict(zip(gene_names, assignment))
                    result = reconcile(gene, species, leaf_mapping)
                    expected_mapping, expected_events, expected_losses = _oracle(
                        gene,
                        species,
                        leaf_mapping,
                    )
                    assert {
                        node_id: node.node_id
                        for node_id, node in result.node_mapping.items()
                    } == {
                        node_id: node.node_id
                        for node_id, node in expected_mapping.items()
                    }
                    assert {
                        node_id: event.event_type
                        for node_id, event in result.events.items()
                    } == expected_events
                    assert Counter(
                        (event.species_parent_id, event.species_child_id)
                        for event in result.losses
                    ) == expected_losses
