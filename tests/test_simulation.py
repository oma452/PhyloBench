import pytest

from phylo_bench.io import parse_newick
from phylo_bench.model.simulation import (
    SimulationEventType,
    SimulationLimitExceeded,
)
from phylo_bench.simulation import reconcile_simulation, simulate_gene_family


def species_tree():
    return parse_newick("((A:1,C:1)X:1,B:1)R;")


def event_signature(result):
    return tuple(
        (
            event.event_type,
            event.gene_node_id,
            event.species_node_id,
            event.species_parent_id,
            event.species_child_id,
            event.lineage_id,
            event.child_lineage_ids,
            event.time_on_branch,
            event.gene_name,
            event.species_name,
        )
        for event in result.events
    )


def test_zero_rates_preserve_root_speciation_and_reconcile():
    result = simulate_gene_family(
        species_tree(), duplication_rate=0, loss_rate=0, seed=7
    )

    assert not result.extinct
    assert result.observed_gene_tree is not None
    assert result.observed_gene_tree.root.name == "ROOT_SPECIATION"
    assert sum(
        event.event_type is SimulationEventType.SPECIATION
        for event in result.events
    ) == 2
    assert len(result.gene_to_species) == 3

    evaluation = reconcile_simulation(result)
    assert evaluation.reconciliation.duplication_count == 0
    assert evaluation.reconciliation.loss_count == 0
    assert evaluation.metrics.duplication_count_error == 0
    assert evaluation.metrics.duplication_branch_compatibility == 1.0
    assert evaluation.metrics.exact_event_location_agreement == 1.0


def test_fixed_seed_is_reproducible():
    first = simulate_gene_family(
        species_tree(), duplication_rate=0.5, loss_rate=0.2, seed=19
    )
    second = simulate_gene_family(
        species_tree(), duplication_rate=0.5, loss_rate=0.2, seed=19
    )

    assert event_signature(first) == event_signature(second)
    assert first.gene_to_species == second.gene_to_species
    assert first.observed_gene_tree is not None
    assert second.observed_gene_tree is not None
    assert [
        (node.node_id, node.name)
        for node in first.observed_gene_tree.nodes()
    ] == [
        (node.node_id, node.name)
        for node in second.observed_gene_tree.nodes()
    ]


def test_observed_trees_remain_binary_across_many_seeds():
    for seed in range(50):
        result = simulate_gene_family(
            species_tree(), duplication_rate=0.8, loss_rate=0.4, seed=seed
        )
        if result.observed_gene_tree is not None:
            assert all(
                node.is_leaf or len(node.children) == 2
                for node in result.observed_gene_tree.nodes()
            )


def test_simulation_records_duplication_and_loss_history():
    result = simulate_gene_family(
        species_tree(), duplication_rate=0.5, loss_rate=0.2, seed=1
    )
    types = {event.event_type for event in result.events}

    assert SimulationEventType.ORIGIN in types
    assert SimulationEventType.SPECIATION in types
    assert SimulationEventType.DUPLICATION in types
    assert SimulationEventType.EXTANT_GENE in types
    assert all(
        event.species_name in {"A", "B", "C"}
        for event in result.events
        if event.event_type is SimulationEventType.EXTANT_GENE
    )
    lost_lineages = {
        event.lineage_id
        for event in result.events
        if event.event_type is SimulationEventType.LOSS
    }
    assert all(
        event.lineage_id not in lost_lineages
        for event in result.events
        if event.event_type is not SimulationEventType.LOSS
    )
    for event in result.events:
        if event.event_type in {
            SimulationEventType.SPECIATION,
            SimulationEventType.DUPLICATION,
        }:
            assert len(event.child_lineage_ids) == 2
            assert len(event.child_gene_node_ids) == 2
    assert all(
        event.species_parent_id is not None
        and event.species_child_id is not None
        for event in result.events
        if event.event_type in {
            SimulationEventType.DUPLICATION,
            SimulationEventType.LOSS,
        }
    )


def test_high_loss_can_produce_an_extinct_family():
    result = simulate_gene_family(
        species_tree(), duplication_rate=0, loss_rate=100, seed=0
    )

    assert result.extinct
    assert result.observed_gene_tree is None
    assert result.gene_to_species == {}
    with pytest.raises(ValueError, match="extinct"):
        reconcile_simulation(result)


@pytest.mark.parametrize(
    "duplication_rate,loss_rate",
    [(-1, 0), (0, -1), (float("inf"), 0), (0, float("nan")), (True, 0)],
)
def test_invalid_rates_are_rejected(duplication_rate, loss_rate):
    with pytest.raises(ValueError, match="finite non-negative"):
        simulate_gene_family(
            species_tree(),
            duplication_rate=duplication_rate,
            loss_rate=loss_rate,
            seed=0,
        )


def test_invalid_seed_is_rejected():
    with pytest.raises(ValueError, match="seed must be an integer"):
        simulate_gene_family(
            species_tree(),
            duplication_rate=0,
            loss_rate=0,
            seed=1.5,
        )


def test_result_mapping_is_immutable():
    result = simulate_gene_family(
        species_tree(), duplication_rate=0, loss_rate=0, seed=0
    )
    with pytest.raises(TypeError):
        result.gene_to_species["new"] = "A"


def test_event_times_stay_inside_their_species_branches():
    result = simulate_gene_family(
        species_tree(), duplication_rate=0.8, loss_rate=0.4, seed=3
    )
    species_nodes = {node.node_id: node for node in species_tree().nodes()}
    for event in result.events:
        if event.event_type not in {
            SimulationEventType.DUPLICATION,
            SimulationEventType.LOSS,
        }:
            continue
        branch_child = species_nodes[event.species_child_id]
        assert event.time_on_branch is not None
        assert 0 <= event.time_on_branch < branch_child.branch_length


def test_configured_limits_stop_high_rate_simulations():
    with pytest.raises(SimulationLimitExceeded, match="max_events"):
        simulate_gene_family(
            species_tree(),
            duplication_rate=100,
            loss_rate=0,
            seed=0,
            max_events=20,
        )
    with pytest.raises(SimulationLimitExceeded, match="max_lineages"):
        simulate_gene_family(
            species_tree(),
            duplication_rate=100,
            loss_rate=0,
            seed=0,
            max_lineages=5,
        )


@pytest.mark.parametrize(
    "newick",
    ["((A,C)X,B)R;", "((A:0,C:1)X:1,B:1)R;", "((A:1,C:1)X:1,B:0)R;"],
)
def test_missing_or_nonpositive_branch_lengths_are_rejected(newick):
    with pytest.raises(ValueError, match="positive finite branch lengths"):
        simulate_gene_family(
            parse_newick(newick),
            duplication_rate=0,
            loss_rate=0,
            seed=0,
        )
