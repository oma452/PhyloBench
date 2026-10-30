"""Reconciliation diagrams with event colors and inferred losses."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib.figure import Figure
from phylo_bench.model.reconciliation import EventType, ReconciliationResult
from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.visualization.trees import _layout


def plot_reconciliation(
    gene_tree: RootedTree,
    species_tree: RootedTree,
    result: ReconciliationResult,
) -> Figure:
    """Render gene/species trees, mapping lines, event colors, and losses."""

    figure, (gene_ax, species_ax) = plt.subplots(1, 2, figsize=(12, 6))
    gene_pos = _layout(gene_tree)
    species_pos = _layout(species_tree)
    event_colors = {
        EventType.DUPLICATION: "tab:red",
        EventType.SPECIATION: "tab:green",
        EventType.LEAF: "tab:blue",
    }

    def draw_tree(ax, tree, positions, title):
        for node in tree.nodes():
            x, y = positions[node.node_id]
            for child in node.children:
                cx, cy = positions[child.node_id]
                ax.plot((x, cx), (cy, cy), color="0.6")
                ax.plot((x, x), (y, cy), color="0.6")
            ax.scatter([x], [y], color="0.3", zorder=3)
            ax.text(x + 0.05, y, node.name or f"node_{node.node_id}", va="center")
        ax.set_title(title)
        ax.set_axis_off()

    draw_tree(gene_ax, gene_tree, gene_pos, "Gene tree / inferred events")
    draw_tree(species_ax, species_tree, species_pos, "Species tree")
    max_depth = max((x for x, _ in gene_pos.values()), default=1)
    for node_id, mapped in result.node_mapping.items():
        if node_id not in gene_pos:
            continue
        gx, gy = gene_pos[node_id]
        sx, sy = species_pos[mapped.node_id]
        event = result.events[node_id].event_type
        gene_ax.scatter([gx], [gy], color=event_colors[event], zorder=4)
        gene_ax.add_artist(
            ConnectionPatch(
                xyA=(gx, gy),
                coordsA=gene_ax.transData,
                xyB=(sx, sy),
                coordsB=species_ax.transData,
                color=event_colors[event],
                alpha=0.35,
                linewidth=0.8,
                arrowstyle="-",
            )
        )
    species_by_id = {node.node_id: node for node in species_tree.nodes()}
    for loss in result.losses:
        parent = species_by_id[loss.species_parent_id]
        child = species_by_id[loss.species_child_id]
        px, py = species_pos[parent.node_id]
        cx, cy = species_pos[child.node_id]
        species_ax.plot((px, cx), (py, cy), color="black", linewidth=3, alpha=0.65)
    figure.suptitle(
        f"Reconciliation: {result.duplication_count} duplications, "
        f"{result.loss_count} losses, cost {result.total_cost:g}"
    )
    return figure
