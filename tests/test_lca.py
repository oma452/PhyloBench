import pytest

from phylo_bench.algorithms import LCAIndex
from phylo_bench.io import parse_newick


def test_lca_hand_verifiable_tree():
    tree = parse_newick("((A,B)X,(C,D)Y)R;")
    index = LCAIndex(tree)
    nodes = {node.name: node for node in tree.nodes() if node.name}

    assert index.lca(nodes["A"], nodes["B"]).name == "X"
    assert index.lca(nodes["A"], nodes["C"]).name == "R"
    assert index.lca(nodes["C"], nodes["D"]).name == "Y"
    assert index.lca(nodes["R"], nodes["A"]).name == "R"


def test_lca_rejects_node_from_another_tree():
    tree = parse_newick("(A,B)R;")
    other = parse_newick("(A,B)R;")
    index = LCAIndex(tree)
    with pytest.raises(ValueError, match="indexed tree"):
        index.lca(tree.root, other.root)


def test_lca_identical_root_and_ancestor_descendant_nodes():
    tree = parse_newick("(((A,B)X,C)Y,D)R;")
    index = LCAIndex(tree)
    nodes = {node.name: node for node in tree.nodes() if node.name}

    assert index.lca(nodes["A"], nodes["A"]) is nodes["A"]
    assert index.lca(nodes["R"], nodes["R"]) is nodes["R"]
    assert index.lca(nodes["Y"], nodes["A"]) is nodes["Y"]
    assert index.lca(nodes["A"], nodes["Y"]) is nodes["Y"]


def test_indexed_tree_rejects_topology_mutation():
    tree = parse_newick("(A,B)R;")
    index = LCAIndex(tree)
    with pytest.raises(RuntimeError, match="frozen"):
        tree.root.add_child(tree.node_by_id(1))
    with pytest.raises(AttributeError):
        tree.root = tree.node_by_id(1)
