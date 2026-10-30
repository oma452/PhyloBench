import pytest

from phylo_bench.io import parse_newick


@pytest.mark.parametrize("newick", ["", "   "])
def test_empty_newick_is_rejected(newick):
    with pytest.raises(ValueError, match="must not be empty"):
        parse_newick(newick)


def test_duplicate_leaf_names_are_rejected():
    with pytest.raises(ValueError, match="duplicate leaf name"):
        parse_newick("(A,A)R;")


@pytest.mark.parametrize("newick", ["(A,B)", "(A,B); trailing", "(A,B);(C,D);"])
def test_semicolon_policy_is_strict(newick):
    with pytest.raises(ValueError):
        parse_newick(newick)


@pytest.mark.parametrize(
    "newick",
    ["(A:-1,B)R;", "(A:nan,B)R;", "(A:inf,B)R;", "(A:abc,B)R;", "(A:,B)R;"],
)
def test_invalid_branch_lengths_are_rejected(newick):
    with pytest.raises(ValueError):
        parse_newick(newick)


def test_colon_in_quoted_label_is_not_a_branch_length():
    tree = parse_newick("('A:foo',B)R;")

    assert [node.name for node in tree.nodes() if node.is_leaf] == ["A:foo", "B"]


def test_quoted_labels_support_punctuation_spaces_and_escaped_quotes():
    tree = parse_newick("('A, (B):foo; [x] ''quoted''',B)R;")

    assert [node.name for node in tree.nodes() if node.is_leaf] == [
        "A, (B):foo; [x] 'quoted'",
        "B",
    ]


def test_colons_and_semicolons_in_comments_are_not_newick_syntax():
    tree = parse_newick("(A[comment: not a branch; [nested: x]],B)R;")

    assert [node.name for node in tree.nodes() if node.is_leaf] == ["A", "B"]


def test_deep_newick_parses_iteratively_and_supports_lca():
    depth = 2000
    newick = "A"
    for index in range(depth):
        newick = f"({newick},B{index})N{index}"
    tree = parse_newick(newick + ";")
    from phylo_bench.algorithms import LCAIndex, postorder

    nodes = {node.name: node for node in tree.nodes() if node.name}
    assert sum(1 for _ in postorder(tree)) == 2 * depth + 1
    assert LCAIndex(tree).lca(nodes["A"], nodes["B1999"]).name == "N1999"
