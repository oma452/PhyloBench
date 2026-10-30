import pytest

from phylo_bench.io import parse_newick
from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


def test_parse_preserves_names_lengths_and_parent_links():
    tree = parse_newick("(A:0.1,(B:0.2,C:0.3)internal:0.4)root;")
    assert tree.root.name == "root"
    assert tree.root.children[0].name == "A"
    assert tree.root.children[0].branch_length == 0.1
    assert tree.root.children[1].parent is tree.root
    assert tree.root.children[1].name == "internal"
    assert tree.root.children[1].branch_length == 0.4


def test_validation_rejects_duplicate_leaf_names():
    root = TreeNode(0)
    root.add_child(TreeNode(1, "A"))
    root.add_child(TreeNode(2, "A"))
    with pytest.raises(ValueError, match="duplicate leaf name"):
        validate_tree(RootedTree(root))


def test_validation_rejects_unnamed_leaf_by_default():
    with pytest.raises(ValueError, match="has no name"):
        validate_tree(RootedTree(TreeNode(0)))


def test_validation_allows_unnamed_leaf_when_requested():
    validate_tree(RootedTree(TreeNode(0)), require_named_leaves=False)


def test_duplicate_child_attachment_is_rejected():
    root = TreeNode(0)
    child = TreeNode(1, "A")
    root.add_child(child)
    with pytest.raises(ValueError, match="already attached"):
        root.add_child(child)


def test_children_and_parent_are_read_only():
    root = TreeNode(0)
    child = TreeNode(1, "A")
    root.add_child(child)
    with pytest.raises(AttributeError):
        root.children.append(TreeNode(2, "B"))
    with pytest.raises(AttributeError):
        child.parent = None


def test_binary_validation_is_opt_in():
    tree = parse_newick("(A,B,C)R;")
    validate_tree(tree)
    with pytest.raises(ValueError, match="not binary"):
        validate_tree(tree, require_binary=True)


def test_invalid_node_id_and_branch_length_are_rejected():
    root = TreeNode(-1, "A", -0.1)
    with pytest.raises(ValueError, match="negative ID"):
        validate_tree(RootedTree(root))
