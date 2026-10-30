"""Reproducible experiments over the Phase 3 simulator."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from statistics import mean, stdev

from phylo_bench.io import parse_newick
from phylo_bench.model.simulation import SimulationEventType
from phylo_bench.simulation import reconcile_simulation, simulate_gene_family


@dataclass(frozen=True, slots=True)
class ExperimentConfig:
    """Complete input configuration for a reproducible experiment."""

    species_tree_newick: str
    seeds: tuple[int, ...]
    duplication_rate: float
    loss_rate: float
    duplication_cost: float = 1.0
    loss_cost: float = 1.0
    max_events: int = 100_000
    max_lineages: int = 10_000

    def __post_init__(self) -> None:
        if not isinstance(self.species_tree_newick, str) or not self.species_tree_newick.strip():
            raise ValueError("species_tree_newick must be a non-empty string")
        try:
            seeds = tuple(self.seeds)
        except TypeError as exc:
            raise ValueError("seeds must be a non-empty iterable of integers") from exc
        object.__setattr__(self, "seeds", seeds)
        if not self.seeds:
            raise ValueError("at least one seed is required")
        if any(isinstance(seed, bool) or not isinstance(seed, int) for seed in self.seeds):
            raise ValueError("seeds must be integers")
        if len(set(self.seeds)) != len(self.seeds):
            raise ValueError("seeds must be unique")
        _validate_rate("duplication_rate", self.duplication_rate)
        _validate_rate("loss_rate", self.loss_rate)
        _validate_rate("duplication_cost", self.duplication_cost)
        _validate_rate("loss_cost", self.loss_cost)
        _validate_limit("max_events", self.max_events)
        _validate_limit("max_lineages", self.max_lineages)


@dataclass(frozen=True, slots=True)
class ExperimentRow:
    """One measured trial."""

    seed: int
    extinct: bool
    complete_duplications: int
    complete_losses: int
    inferred_duplications: int | None
    inferred_losses: int | None
    inferred_cost: float | None
    identifiable_loss_count: int | None
    exact_event_location_agreement: float | None
    duplication_branch_compatibility: float | None
    loss_branch_precision: float | None
    loss_branch_recall: float | None


@dataclass(frozen=True, slots=True)
class ExperimentSummary:
    """Measured mean and standard deviation across non-extinct trials."""

    trials: int
    extinct_trials: int
    metrics: dict[str, tuple[float | None, float | None]]
    config: ExperimentConfig
    seeds: tuple[int, ...]


def run_experiment(config: ExperimentConfig) -> tuple[tuple[ExperimentRow, ...], ExperimentSummary]:
    """Run exactly the configured seeds and summarize actual outputs."""

    rows: list[ExperimentRow] = []
    for seed in config.seeds:
        simulation = simulate_gene_family(
            parse_newick(config.species_tree_newick),
            duplication_rate=config.duplication_rate,
            loss_rate=config.loss_rate,
            seed=seed,
            max_events=config.max_events,
            max_lineages=config.max_lineages,
        )
        if simulation.extinct:
            rows.append(
                ExperimentRow(
                    seed,
                    True,
                    sum(event.event_type is SimulationEventType.DUPLICATION for event in simulation.events),
                    sum(event.event_type is SimulationEventType.LOSS for event in simulation.events),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                )
            )
            continue
        evaluation = reconcile_simulation(
            simulation,
            duplication_cost=config.duplication_cost,
            loss_cost=config.loss_cost,
        )
        metrics = evaluation.metrics
        rows.append(
            ExperimentRow(
                seed,
                False,
                metrics.complete_duplication_count,
                metrics.complete_loss_count,
                metrics.inferred_duplication_count,
                metrics.inferred_loss_count,
                metrics.inferred_cost,
                metrics.identifiable_loss_count,
                metrics.exact_event_location_agreement,
                metrics.duplication_branch_compatibility,
                metrics.loss_branch_precision,
                metrics.loss_branch_recall,
            )
        )
    usable = [row for row in rows if not row.extinct]
    names = (
        "complete_duplications", "complete_losses", "inferred_duplications",
        "inferred_losses", "inferred_cost", "exact_event_location_agreement",
        "identifiable_loss_count",
        "duplication_branch_compatibility", "loss_branch_precision",
        "loss_branch_recall",
    )
    summary: dict[str, tuple[float | None, float | None]] = {}
    for name in names:
        source = rows if name in {"complete_duplications", "complete_losses"} else usable
        values = [float(getattr(row, name)) for row in source if getattr(row, name) is not None]
        summary[name] = (
            (mean(values), stdev(values) if len(values) > 1 else 0.0)
            if values else (None, None)
        )
    return tuple(rows), ExperimentSummary(
        len(rows), len(rows) - len(usable), summary, config, config.seeds
    )


def _validate_rate(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite non-negative number")
    if not isfinite(float(value)) or value < 0:
        raise ValueError(f"{name} must be a finite non-negative number")


def _validate_limit(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
