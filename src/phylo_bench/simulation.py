"""Deterministic continuous-time duplication-loss gene-family simulation.

On a branch of length ``t`` with ``n`` active gene lineages, the next event
waits an exponential time with rate ``n(lambda_D + lambda_L)``. Conditional
on an event, duplication has probability
``lambda_D / (lambda_D + lambda_L)`` and loss has the complementary
probability. Species-tree speciation is deterministic at branch endpoints.
"""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Iterable
from dataclasses import dataclass
from math import isfinite
from random import Random
from types import MappingProxyType

from phylo_bench.algorithms.reconciliation import reconcile
from phylo_bench.model.reconciliation import EventType, ReconciliationResult
from phylo_bench.model.simulation import (
    EvaluationMetrics,
    SimulationEvaluation,
    SimulationEvent,
    SimulationEventType,
    SimulationLimitExceeded,
    SimulationResult,
)
from phylo_bench.model.tree import RootedTree, TreeNode
from phylo_bench.validation import validate_tree


@dataclass(slots=True)
class _Lineage:
    lineage_id: int
    parent_lineage_id: int | None
    gene_node: TreeNode


@dataclass(frozen=True, slots=True)
class _BranchTask:
    parent_species: TreeNode
    child_species: TreeNode
    lineages: tuple[_Lineage, ...]


def simulate_gene_family(
    species_tree: RootedTree,
    *,
    duplication_rate: float,
    loss_rate: float,
    seed: int,
    max_events: int = 100_000,
    max_lineages: int = 10_000,
) -> SimulationResult:
    """Simulate one family on a fixed rooted, binary, time-scaled tree.

    Each active lineage on a branch experiences duplication and loss as
    independent exponential processes with rates ``duplication_rate`` and
    ``loss_rate``. Species-tree speciations are deterministic. Branch lengths
    must be positive finite values on every non-root species node. The
    virtual origin and biological root speciation are retained in the complete
    history; the observed tree removes the virtual origin, extinct-only
    subtrees, and unary event nodes.
    """

    _validate_rate("duplication_rate", duplication_rate)
    _validate_rate("loss_rate", loss_rate)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    _validate_limit("max_events", max_events)
    _validate_limit("max_lineages", max_lineages)
    validate_simulation_species_tree(species_tree)

    rng = Random(seed)
    events: list[SimulationEvent] = []
    next_event_id = 0
    next_lineage_id = 0
    next_gene_node_id = 0

    def record(event_type: SimulationEventType, **kwargs: object) -> None:
        nonlocal next_event_id
        if next_event_id >= max_events:
            raise SimulationLimitExceeded("max_events limit exceeded")
        events.append(SimulationEvent(next_event_id, event_type, **kwargs))
        next_event_id += 1

    origin = TreeNode(next_gene_node_id, "ORIGIN")
    next_gene_node_id += 1
    root_speciation = TreeNode(next_gene_node_id, "ROOT_SPECIATION")
    next_gene_node_id += 1
    origin.add_child(root_speciation)
    origin_lineage = _Lineage(next_lineage_id, None, root_speciation)
    next_lineage_id += 1
    record(
        SimulationEventType.ORIGIN,
        gene_node_id=origin.node_id,
        lineage_id=origin_lineage.lineage_id,
        species_node_id=species_tree.root.node_id,
    )

    queue: deque[_BranchTask] = deque()
    root_children = tuple(species_tree.root.children)
    root_lineages: list[_Lineage] = []
    for child in root_children:
        if next_lineage_id + 1 > max_lineages:
            raise SimulationLimitExceeded("max_lineages limit exceeded")
        anchor = TreeNode(next_gene_node_id, "LINEAGE")
        next_gene_node_id += 1
        root_speciation.add_child(anchor)
        lineage = _Lineage(next_lineage_id, origin_lineage.lineage_id, anchor)
        next_lineage_id += 1
        root_lineages.append(lineage)
        queue.append(_BranchTask(species_tree.root, child, (lineage,)))
    record(
        SimulationEventType.SPECIATION,
        gene_node_id=root_speciation.node_id,
        species_node_id=species_tree.root.node_id,
        lineage_id=origin_lineage.lineage_id,
        child_lineage_ids=tuple(lineage.lineage_id for lineage in root_lineages),
        child_gene_node_ids=tuple(lineage.gene_node.node_id for lineage in root_lineages),
    )

    while queue:
        task = queue.popleft()
        survivors = _simulate_branch(
            task,
            duplication_rate,
            loss_rate,
            rng,
            events,
            record,
            next_lineage_id,
            next_gene_node_id,
            max_lineages,
        )
        next_lineage_id = survivors.next_lineage_id
        next_gene_node_id = survivors.next_gene_node_id
        live = survivors.lineages

        if task.child_species.is_leaf:
            for lineage in live:
                gene_name = f"gene_{lineage.lineage_id}"
                leaf = TreeNode(next_gene_node_id, gene_name)
                next_gene_node_id += 1
                lineage.gene_node.add_child(leaf)
                record(
                    SimulationEventType.EXTANT_GENE,
                    gene_node_id=leaf.node_id,
                    lineage_id=lineage.lineage_id,
                    species_node_id=task.child_species.node_id,
                    gene_name=gene_name,
                    species_name=task.child_species.name,
                )
        else:
            children = tuple(task.child_species.children)
            for lineage in live:
                speciation = TreeNode(next_gene_node_id, "SPECIATION")
                next_gene_node_id += 1
                lineage.gene_node.add_child(speciation)
                child_lineages = []
                for child_species in children:
                    if next_lineage_id + 1 > max_lineages:
                        raise SimulationLimitExceeded(
                            "max_lineages limit exceeded"
                        )
                    child_anchor = TreeNode(next_gene_node_id, "LINEAGE")
                    next_gene_node_id += 1
                    speciation.add_child(child_anchor)
                    child = _Lineage(
                        next_lineage_id,
                        lineage.lineage_id,
                        child_anchor,
                    )
                    next_lineage_id += 1
                    child_lineages.append(child)
                    queue.append(
                        _BranchTask(task.child_species, child_species, (child,))
                    )
                record(
                    SimulationEventType.SPECIATION,
                    gene_node_id=speciation.node_id,
                    species_node_id=task.child_species.node_id,
                    lineage_id=lineage.lineage_id,
                    child_lineage_ids=tuple(x.lineage_id for x in child_lineages),
                    child_gene_node_ids=tuple(
                        x.gene_node.node_id for x in child_lineages
                    ),
                )

    true_tree = RootedTree(origin)
    observed_tree = _project_observed_tree(true_tree)
    gene_to_species = MappingProxyType({
        event.gene_name: event.species_name
        for event in events
        if event.event_type is SimulationEventType.EXTANT_GENE
        and event.gene_name is not None
        and event.species_name is not None
    })
    return SimulationResult(
        species_tree=species_tree,
        observed_gene_tree=observed_tree,
        gene_to_species=gene_to_species,
        events=tuple(events),
        seed=seed,
        duplication_rate=float(duplication_rate),
        loss_rate=float(loss_rate),
    )


@dataclass(frozen=True, slots=True)
class _BranchResult:
    lineages: tuple[_Lineage, ...]
    next_lineage_id: int
    next_gene_node_id: int


def _simulate_branch(
    task: _BranchTask,
    duplication_rate: float,
    loss_rate: float,
    rng: Random,
    events: list[SimulationEvent],
    record,
    next_lineage_id: int,
    next_gene_node_id: int,
    max_lineages: int,
) -> _BranchResult:
    active = list(task.lineages)
    time = 0.0
    length = task.child_species.branch_length
    total_rate_per_lineage = duplication_rate + loss_rate

    while active and total_rate_per_lineage > 0:
        waiting = rng.expovariate(len(active) * total_rate_per_lineage)
        if time + waiting >= length:
            break
        time += waiting
        index = rng.randrange(len(active))
        lineage = active.pop(index)
        if rng.random() < duplication_rate / total_rate_per_lineage:
            if next_lineage_id + 2 > max_lineages:
                raise SimulationLimitExceeded("max_lineages limit exceeded")
            duplicate = TreeNode(next_gene_node_id, "DUPLICATION")
            next_gene_node_id += 1
            lineage.gene_node.add_child(duplicate)
            left_anchor = TreeNode(next_gene_node_id, "LINEAGE")
            right_anchor = TreeNode(next_gene_node_id + 1, "LINEAGE")
            next_gene_node_id += 2
            duplicate.add_child(left_anchor)
            duplicate.add_child(right_anchor)
            left = _Lineage(next_lineage_id, lineage.lineage_id, left_anchor)
            right = _Lineage(next_lineage_id + 1, lineage.lineage_id, right_anchor)
            next_lineage_id += 2
            active.extend((left, right))
            record(
                SimulationEventType.DUPLICATION,
                gene_node_id=duplicate.node_id,
                species_parent_id=task.parent_species.node_id,
                species_child_id=task.child_species.node_id,
                lineage_id=lineage.lineage_id,
                child_lineage_ids=(left.lineage_id, right.lineage_id),
                child_gene_node_ids=(left_anchor.node_id, right_anchor.node_id),
                time_on_branch=time,
            )
        else:
            record(
                SimulationEventType.LOSS,
                species_parent_id=task.parent_species.node_id,
                species_child_id=task.child_species.node_id,
                lineage_id=lineage.lineage_id,
                gene_node_id=lineage.gene_node.node_id,
                time_on_branch=time,
            )
    return _BranchResult(tuple(active), next_lineage_id, next_gene_node_id)


def _project_observed_tree(true_tree: RootedTree) -> RootedTree | None:
    """Clone the extant projection while preserving original node IDs."""

    def clone(node: TreeNode) -> TreeNode | None:
        if node.is_leaf:
            return TreeNode(node.node_id, node.name) if node.name.startswith("gene_") else None
        children = [clone(child) for child in node.children]
        surviving = [child for child in children if child is not None]
        if not surviving:
            return None
        if len(surviving) == 1:
            return surviving[0]
        projected = TreeNode(node.node_id, node.name)
        for child in surviving:
            projected.add_child(child)
        return projected

    projected_root = clone(true_tree.root)
    if projected_root is None:
        return None
    return RootedTree(projected_root)


def reconcile_simulation(
    simulation: SimulationResult,
    *,
    duplication_cost: float = 1.0,
    loss_cost: float = 1.0,
) -> SimulationEvaluation:
    """Reconcile a non-extinct simulation and compute evaluation metrics."""

    if simulation.observed_gene_tree is None:
        raise ValueError("cannot reconcile an extinct gene family")
    reconciliation = reconcile(
        simulation.observed_gene_tree,
        simulation.species_tree,
        simulation.gene_to_species,
        duplication_cost=duplication_cost,
        loss_cost=loss_cost,
    )
    metrics = evaluate_simulation(simulation, reconciliation)
    return SimulationEvaluation(simulation, reconciliation, metrics)


def evaluate_simulation(
    simulation: SimulationResult,
    reconciliation: ReconciliationResult,
) -> EvaluationMetrics:
    """Compare inferred events with visible and complete simulated history.

    Complete events hidden by extinction or unary projection are reported
    separately. Speciation location agreement is exact-node agreement.
    Duplication location agreement is branch-compatible agreement because an
    extant gene tree generally cannot identify the exact time within a branch.
    """

    visible_nodes = _visible_event_nodes(simulation)
    event_by_node = {
        event.gene_node_id: event
        for event in simulation.events
        if event.gene_node_id is not None
    }
    complete_duplications = sum(
        event.event_type is SimulationEventType.DUPLICATION
        for event in simulation.events
    )
    complete_losses = sum(
        event.event_type is SimulationEventType.LOSS
        for event in simulation.events
    )
    observable_duplications = sum(
        event_by_node[node_id].event_type is SimulationEventType.DUPLICATION
        for node_id in visible_nodes
    )
    observable_speciations = sum(
        event_by_node[node_id].event_type is SimulationEventType.SPECIATION
        for node_id in visible_nodes
    )
    matched = 0
    exact_speciation_matches = 0
    compatible_duplications = 0
    speciation_count = 0
    duplication_count = 0
    for node_id in visible_nodes:
        inferred = reconciliation.events.get(node_id)
        simulated = event_by_node[node_id]
        if inferred is None:
            continue
        expected_type = (
            EventType.DUPLICATION
            if simulated.event_type is SimulationEventType.DUPLICATION
            else EventType.SPECIATION
        )
        matched += inferred.event_type is expected_type
        if inferred.event_type is expected_type:
            if simulated.event_type is SimulationEventType.SPECIATION:
                speciation_count += 1
                exact_speciation_matches += int(
                    inferred.species_node_id == simulated.species_node_id
                )
            else:
                duplication_count += 1
                compatible_duplications += _location_match(
                    inferred, simulated, simulation.species_tree
                )
    identifiable_loss_branches = _identifiable_loss_branches(simulation)
    inferred_loss_branches = [
        (event.species_parent_id, event.species_child_id)
        for event in reconciliation.losses
    ]
    precision, recall = _multiset_precision_recall(
        inferred_loss_branches, identifiable_loss_branches
    )
    return EvaluationMetrics(
        inferred_duplication_count=reconciliation.duplication_count,
        inferred_loss_count=reconciliation.loss_count,
        inferred_cost=reconciliation.total_cost,
        complete_duplication_count=complete_duplications,
        complete_loss_count=complete_losses,
        observable_duplication_count=observable_duplications,
        observable_speciation_count=observable_speciations,
        duplication_count_error=(
            reconciliation.duplication_count - observable_duplications
        ),
        matched_event_count=matched,
        observable_event_count=len(visible_nodes),
        exact_event_location_agreement=(
            exact_speciation_matches / speciation_count
            if speciation_count else 1.0
        ),
        duplication_branch_compatibility=(
            compatible_duplications / duplication_count
            if duplication_count else 1.0
        ),
        loss_branch_precision=precision,
        loss_branch_recall=recall,
        identifiable_loss_count=len(identifiable_loss_branches),
        ambiguity_note=(
            "Complete-loss count includes every simulated loss. Loss branch "
            "precision and recall compare only losses whose split sibling has "
            "an extant descendant, making the loss structurally identifiable "
            "from the observed family. Hidden losses are excluded rather than "
            "treated as false negatives; inferred cost remains an independent "
            "reconciliation quantity."
        ),
    )


def _visible_event_nodes(simulation: SimulationResult) -> set[int]:
    if simulation.observed_gene_tree is None:
        return set()
    return {
        node.node_id
        for node in simulation.observed_gene_tree.nodes()
        if not node.is_leaf
    }


def _identifiable_loss_branches(
    simulation: SimulationResult,
) -> list[tuple[int | None, int | None]]:
    """Return loss branches supported by an extant sibling lineage.

    A simulated loss has an observable structural witness only when its
    lineage was created by a speciation or duplication split and at least one
    sibling lineage from that split has an extant descendant. Losses with no
    such witness are hidden historical events and are not exact-recovery
    targets.
    """

    parent_by_lineage: dict[int, int] = {}
    siblings_by_lineage: dict[int, tuple[int, ...]] = {}
    for event in simulation.events:
        if event.event_type not in (
            SimulationEventType.SPECIATION,
            SimulationEventType.DUPLICATION,
        ) or event.lineage_id is None:
            continue
        for child_id in event.child_lineage_ids:
            parent_by_lineage[child_id] = event.lineage_id
            siblings_by_lineage[child_id] = tuple(
                sibling for sibling in event.child_lineage_ids if sibling != child_id
            )

    extant_lineages = {
        event.lineage_id
        for event in simulation.events
        if event.event_type is SimulationEventType.EXTANT_GENE
        and event.lineage_id is not None
    }
    extant_ancestry: set[int] = set()
    for lineage_id in extant_lineages:
        current: int | None = lineage_id
        while current is not None and current not in extant_ancestry:
            extant_ancestry.add(current)
            current = parent_by_lineage.get(current)

    return [
        (event.species_parent_id, event.species_child_id)
        for event in simulation.events
        if event.event_type is SimulationEventType.LOSS
        and any(
            sibling in extant_ancestry
            for sibling in siblings_by_lineage.get(event.lineage_id, ())
        )
    ]


def _location_match(
    inferred,
    simulated: SimulationEvent,
    species_tree: RootedTree,
) -> int:
    if simulated.event_type is SimulationEventType.SPECIATION:
        return int(inferred.species_node_id == simulated.species_node_id)
    if simulated.species_child_id is None:
        return 0
    inferred_node = species_tree.node_by_id(inferred.species_node_id)
    branch_child = species_tree.node_by_id(simulated.species_child_id)
    current: TreeNode | None = inferred_node
    while current is not None:
        if current is branch_child:
            return 1
        current = current.parent
    return 0


def _multiset_precision_recall(inferred: Iterable, truth: Iterable) -> tuple[float, float]:
    from collections import Counter

    inferred_counts = Counter(inferred)
    truth_counts = Counter(truth)
    matched = sum((inferred_counts & truth_counts).values())
    precision = matched / sum(inferred_counts.values()) if inferred_counts else 1.0
    recall = matched / sum(truth_counts.values()) if truth_counts else 1.0
    return precision, recall


def _validate_rate(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite non-negative number")
    if not isfinite(float(value)) or value < 0:
        raise ValueError(f"{name} must be a finite non-negative number")


def _validate_limit(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


def _validate_branch_lengths(species_tree: RootedTree) -> None:
    for node in species_tree.nodes():
        if node is species_tree.root:
            continue
        if node.branch_length is None or not isfinite(node.branch_length) or node.branch_length <= 0:
            raise ValueError(
                "simulation requires positive finite branch lengths on all "
                "non-root species nodes"
            )


def validate_simulation_species_tree(species_tree: RootedTree) -> None:
    """Validate a species tree supported by the simulator.

    Simulation models a speciation at the root, so a tree with fewer than two
    species cannot be represented by the simulator's binary event process.
    """

    validate_tree(species_tree, require_binary=True)
    if species_tree.root.is_leaf:
        raise ValueError("simulation requires at least two species")
    _validate_branch_lengths(species_tree)
