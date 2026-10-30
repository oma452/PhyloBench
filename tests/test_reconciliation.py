from collections import Counter

import pytest

from phylo_bench.algorithms import map_gene_leaves, reconcile
from phylo_bench.model import EventType
from phylo_bench.io import parse_newick


def run_case(species_newick: str, gene_newick: str):
    return reconcile(parse_newick(gene_newick), parse_newick(species_newick))


def loss_branches(result, species_tree):
    nodes = {node.node_id: node for node in species_tree.nodes()}
    return Counter(
        (nodes[event.species_parent_id].name, nodes[event.species_child_id].name)
        for event in result.losses
    )


def event_for_name(result, gene_tree, name):
    node = next(node for node in gene_tree.nodes() if node.name == name)
    return result.events[node.node_id].event_type


def test_single_gene_in_a():
    species = parse_newick("(A,B)R;")
    gene = parse_newick("A1;")
    result = reconcile(species_tree=species, gene_tree=gene,
                       leaf_to_species={"A1": "A"})

    assert result.node_mapping[gene.root.node_id].name == "A"
    assert result.duplication_count == 0
    assert result.loss_count == 1
    assert result.total_cost == 1
    assert loss_branches(result, species) == Counter({("R", "B"): 1})


def test_one_copy_in_both_a_and_b():
    species = parse_newick("(A,B)R;")
    gene = parse_newick("(A1,B1)s;")
    result = reconcile(gene, species, {"A1": "A", "B1": "B"})

    assert event_for_name(result, gene, "s") is EventType.SPECIATION
    assert result.duplication_count == 0
    assert result.loss_count == 0
    assert result.total_cost == 0


def test_duplication_in_a():
    species = parse_newick("(A,B)R;")
    gene = parse_newick("(A1,A2)d;")
    result = reconcile(gene, species, {"A1": "A", "A2": "A"})

    assert event_for_name(result, gene, "d") is EventType.DUPLICATION
    assert result.duplication_count == 1
    assert result.loss_count == 1
    assert result.total_cost == 2
    assert loss_branches(result, species) == Counter({("R", "B"): 1})


def test_duplication_before_speciation():
    species = parse_newick("((A,C)X,B)R;")
    gene = parse_newick("((A1,C1)s,A2)d;")
    result = reconcile(gene, species, {"A1": "A", "C1": "C", "A2": "A"})

    assert event_for_name(result, gene, "s") is EventType.SPECIATION
    assert event_for_name(result, gene, "d") is EventType.DUPLICATION
    assert result.duplication_count == 1
    assert result.loss_count == 2
    assert result.total_cost == 3
    assert loss_branches(result, species) == Counter({
        ("R", "B"): 1,
        ("X", "C"): 1,
    })


def test_duplication_followed_by_loss():
    species = parse_newick("((A,C)X,B)R;")
    gene = parse_newick("(A1,(A2,C1)s)d;")
    result = reconcile(gene, species, {"A1": "A", "A2": "A", "C1": "C"})

    assert result.duplication_count == 1
    assert result.loss_count == 2
    assert result.total_cost == 3
    assert loss_branches(result, species) == Counter({
        ("R", "B"): 1,
        ("X", "C"): 1,
    })


def test_independent_copies_lost_on_same_branch_remain_separate():
    species = parse_newick("((A,C)X,B)R;")
    gene = parse_newick("(((A1,C1)s,A2)d2,A3)d1;")
    result = reconcile(
        gene,
        species,
        {"A1": "A", "C1": "C", "A2": "A", "A3": "A"},
    )

    assert result.duplication_count == 2
    assert result.loss_count == 3
    assert result.total_cost == 5
    assert loss_branches(result, species) == Counter({
        ("R", "B"): 1,
        ("X", "C"): 2,
    })


def test_explicit_gene_to_species_mapping():
    gene = parse_newick("(gene_a,gene_b)g;")
    species = parse_newick("(A,B)R;")
    mapping = map_gene_leaves(
        gene, species, {"gene_a": "A", "gene_b": "B"}
    )
    assert {node.name for node in mapping.values()} == {"A", "B"}


def test_weighted_cost_keeps_counts_separate():
    result = reconcile(
        parse_newick("(A1,A2)d;"),
        parse_newick("(A,B)R;"),
        {"A1": "A", "A2": "A"},
        duplication_cost=2.5,
        loss_cost=4,
    )
    assert result.duplication_count == 1
    assert result.loss_count == 1
    assert result.total_cost == 6.5


def test_missing_or_unknown_species_mapping_is_rejected():
    gene = parse_newick("(gene_a,gene_b)g;")
    species = parse_newick("(A,B)R;")
    with pytest.raises(ValueError, match="no species mapping"):
        map_gene_leaves(gene, species, {"gene_a": "A"})
    with pytest.raises(ValueError, match="unknown species"):
        map_gene_leaves(gene, species, {"gene_a": "A", "gene_b": "C"})
