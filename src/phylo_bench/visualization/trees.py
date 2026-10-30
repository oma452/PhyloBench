"""Deterministic Matplotlib renderers for rooted trees."""

from __future__ import annotations

from collections.abc import Iterable

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from phylo_bench.model.tree import RootedTree, TreeNode


def _layout(tree: RootedTree) -> dict[int, tuple[float, float]]:
    positions: dict[int, tuple[float, float]] = {}
    next_y = 0

    def visit(node: TreeNode, depth: int) -> float:
        nonlocal next_y
        if node.is_leaf:
            y = float(next_y)
            next_y += 1
        else:
            ys = [visit(child, depth + 1) for child in node.children]
            y = sum(ys) / len(ys)
        positions[node.node_id] = (float(depth), y)
        return y

    visit(tree.root, 0)
    return positions


def _draw_tree(ax: Axes, tree: RootedTree, title: str, colors: dict[int, str] | None = None) -> None:
    positions = _layout(tree)
    for node in tree.nodes():
        x, y = positions[node.node_id]
        for child in node.children:
            cx, cy = positions[child.node_id]
            ax.plot((x, cx), (cy, cy), color="0.35", linewidth=1)
            ax.plot((x, x), (y, cy), color="0.35", linewidth=1)
        color = (colors or {}).get(node.node_id, "tab:blue")
        ax.scatter([x], [y], color=color, zorder=3)
        label = node.name or f"node_{node.node_id}"
        ax.text(x + 0.05, y, label, va="center", fontsize=9)
    ax.set_title(title)
    ax.set_axis_off()


def plot_species_tree(tree: RootedTree, *, ax: Axes | None = None) -> Figure:
    """Render a rooted species tree and return its owning figure."""

    if ax is None:
        figure, ax = __import__("matplotlib.pyplot", fromlist=["subplots"]).subplots()
    else:
        figure = ax.figure
    _draw_tree(ax, tree, "Species tree")
    return figure


def plot_gene_tree(tree: RootedTree, *, ax: Axes | None = None) -> Figure:
    """Render a rooted observed gene tree and return its owning figure."""

    if ax is None:
        figure, ax = __import__("matplotlib.pyplot", fromlist=["subplots"]).subplots()
    else:
        figure = ax.figure
    _draw_tree(ax, tree, "Gene tree")
    return figure
