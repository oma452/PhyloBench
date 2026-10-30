"""Typed outputs for deterministic gene-family simulation and evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from collections.abc import Mapping

from phylo_bench.model.reconciliation import ReconciliationResult
from phylo_bench.model.tree import RootedTree


class SimulationEventType(str, Enum):
    """Events retained in a complete simulated history."""

    ORIGIN = "origin"
    SPECIATION = "speciation"
    DUPLICATION = "duplication"
    LOSS = "loss"
    EXTANT_GENE = "extant_gene"


class SimulationLimitExceeded(RuntimeError):
    """Raised when a configured simulation resource limit is reached."""


@dataclass(frozen=True, slots=True)
class SimulationEvent:
    """One complete-history event."""

    event_id: int
    event_type: SimulationEventType
    gene_node_id: int | None = None
    species_node_id: int | None = None
    species_parent_id: int | None = None
    species_child_id: int | None = None
    lineage_id: int | None = None
    child_lineage_ids: tuple[int, ...] = ()
    child_gene_node_ids: tuple[int, ...] = ()
    time_on_branch: float | None = None
    gene_name: str | None = None
    species_name: str | None = None


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """A simulated family, its observed projection, and complete history."""

    species_tree: RootedTree
    observed_gene_tree: RootedTree | None
    gene_to_species: Mapping[str, str]
    events: tuple[SimulationEvent, ...]
    seed: int
    duplication_rate: float
    loss_rate: float

    @property
    def extinct(self) -> bool:
        """Return whether no extant gene copy survived."""

        return self.observed_gene_tree is None


@dataclass(frozen=True, slots=True)
class EvaluationMetrics:
    """Comparison of an inferred reconciliation with simulated history."""

    inferred_duplication_count: int
    inferred_loss_count: int
    inferred_cost: float
    complete_duplication_count: int
    complete_loss_count: int
    observable_duplication_count: int
    observable_speciation_count: int
    duplication_count_error: int
    matched_event_count: int
    observable_event_count: int
    exact_event_location_agreement: float
    duplication_branch_compatibility: float
    loss_branch_precision: float
    loss_branch_recall: float
    ambiguity_note: str


@dataclass(frozen=True, slots=True)
class SimulationEvaluation:
    """Simulation output paired with reconciliation and metrics."""

    simulation: SimulationResult
    reconciliation: ReconciliationResult
    metrics: EvaluationMetrics
