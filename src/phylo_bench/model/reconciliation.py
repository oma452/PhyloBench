"""Typed results for classical duplication-loss reconciliation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from phylo_bench.model.tree import TreeNode


class EventType(str, Enum):
    """Events represented on internal gene-tree nodes."""

    LEAF = "leaf"
    SPECIATION = "speciation"
    DUPLICATION = "duplication"


@dataclass(frozen=True, slots=True)
class ReconciliationEvent:
    """The event assigned to one gene-tree node."""

    gene_node_id: int
    event_type: EventType
    species_node_id: int


@dataclass(frozen=True, slots=True)
class LossEvent:
    """One loss of one gene lineage on one species-tree branch."""

    species_parent_id: int
    species_child_id: int
    gene_parent_id: int | None
    gene_child_id: int
    transition: str


@dataclass(frozen=True, slots=True)
class ReconciliationResult:
    """Complete output of a classical LCA duplication-loss reconciliation."""

    leaf_mapping: dict[int, TreeNode]
    node_mapping: dict[int, TreeNode]
    events: dict[int, ReconciliationEvent]
    losses: tuple[LossEvent, ...]
    duplication_count: int
    loss_count: int
    duplication_cost: float
    loss_cost: float
    total_cost: float
