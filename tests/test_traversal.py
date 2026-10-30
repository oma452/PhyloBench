from phylo_bench.algorithms import postorder, preorder
from phylo_bench.io import parse_newick


def test_iterative_preorder_and_postorder():
    tree = parse_newick("((A,B)X,C)R;")
    assert [node.name for node in preorder(tree)] == ["R", "X", "A", "B", "C"]
    assert [node.name for node in postorder(tree)] == ["A", "B", "X", "C", "R"]
