"""Analytical special-case checks for the stochastic simulator."""

from __future__ import annotations

import math
from statistics import mean

from phylo_bench.io import parse_newick
from phylo_bench.simulation import simulate_gene_family


def test_pure_loss_extinction_frequency_matches_two_branch_expectation():
    """For two unit branches, pure-loss extinction is (1 - exp(-lambda))^2."""

    rate = 0.5
    expected = (1.0 - math.exp(-rate)) ** 2
    observations = []
    for seed in range(3000):
        result = simulate_gene_family(
            parse_newick("(A:1,B:1)R;"),
            duplication_rate=0.0,
            loss_rate=rate,
            seed=seed,
        )
        observations.append(result.extinct)
    observed = mean(observations)
    assert abs(observed - expected) < 0.025


def test_pure_duplication_mean_matches_yule_expectation():
    """For two unit branches, expected duplications are 2 * (exp(lambda)-1)."""

    rate = 0.2
    expected = 2.0 * (math.exp(rate) - 1.0)
    observations = []
    for seed in range(3000):
        result = simulate_gene_family(
            parse_newick("(A:1,B:1)R;"),
            duplication_rate=rate,
            loss_rate=0.0,
            seed=seed,
        )
        observations.append(
            sum(
                event.event_type.value == "duplication"
                for event in result.events
            )
        )
    observed = mean(observations)
    assert abs(observed - expected) < 0.04
