"""Phase 4 visualization, experiment, export, and report tests."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from phylo_bench.algorithms.reconciliation import reconcile
from phylo_bench.export import write_figure, write_rows_csv, write_summary_csv
from phylo_bench.experiments import ExperimentConfig, run_experiment
from phylo_bench.io import parse_newick
from phylo_bench.report import write_research_report
from phylo_bench.simulation import reconcile_simulation, simulate_gene_family
from phylo_bench.visualization import (
    plot_gene_tree,
    plot_metric_summary,
    plot_reconciliation,
    plot_species_tree,
)


def _trees():
    species = parse_newick("(A:1,B:1)S;")
    gene = parse_newick("(A:1,B:1)G;")
    return species, gene


def test_experiment_is_reproducible_and_reports_measured_values():
    config = ExperimentConfig("(A:1,B:1)S;", (7, 11), 0.0, 0.0)
    first = run_experiment(config)
    second = run_experiment(config)
    assert first == second
    rows, summary = first
    assert len(rows) == 2
    assert summary.extinct_trials == 0
    assert summary.metrics["inferred_cost"] == (0.0, 0.0)
    assert summary.seeds == (7, 11)


def test_experiment_configuration_rejects_invalid_values():
    import pytest
    invalid = [
        {"seeds": 1},
        {"seeds": None},
        {"seeds": ()},
        {"seeds": (1, 1)},
        {"seeds": (True,)},
        {"duplication_rate": -1.0},
        {"loss_rate": float("nan")},
        {"max_events": 0},
        {"max_lineages": 0},
    ]
    for overrides in invalid:
        values = {
            "species_tree_newick": "(A:1,B:1)S;",
            "seeds": (1,),
            "duplication_rate": 0.0,
            "loss_rate": 0.0,
        }
        values.update(overrides)
        with pytest.raises(ValueError):
            ExperimentConfig(**values)


def test_nonzero_rate_summary_has_sample_standard_deviation():
    config = ExperimentConfig("(A:1,B:1)S;", (1, 2, 3), 0.2, 0.0)
    rows, summary = run_experiment(config)
    assert any(not row.extinct for row in rows)
    mean_value, std_value = summary.metrics["complete_duplications"]
    assert mean_value is not None
    assert std_value is not None


def test_all_extinct_trials_keep_inferred_metrics_unavailable():
    config = ExperimentConfig("(A:1,B:1)S;", (1, 2), 0.0, 100.0)
    rows, summary = run_experiment(config)
    assert all(row.extinct for row in rows)
    assert summary.metrics["inferred_cost"] == (None, None)


def test_loss_metrics_exclude_hidden_losses_from_recall():
    visible = simulate_gene_family(
        parse_newick("(A:1,B:1)S;"),
        duplication_rate=0.0,
        loss_rate=1.0,
        seed=5,
    )
    visible_evaluation = reconcile_simulation(visible)
    assert visible_evaluation.metrics.complete_loss_count == 1
    assert visible_evaluation.metrics.identifiable_loss_count == 1
    assert visible_evaluation.metrics.loss_branch_recall == 1.0

    hidden = simulate_gene_family(
        parse_newick("(A:1,B:1)S;"),
        duplication_rate=3.0,
        loss_rate=3.0,
        seed=0,
    )
    hidden_evaluation = reconcile_simulation(hidden)
    assert (
        hidden_evaluation.metrics.complete_loss_count
        > hidden_evaluation.metrics.identifiable_loss_count
    )


def test_exports_and_report_use_actual_rows(tmp_path: Path):
    config = ExperimentConfig("(A:1,B:1)S;", (3,), 0.0, 0.0)
    rows, summary = run_experiment(config)
    csv_path = tmp_path / "trials.csv"
    summary_csv_path = tmp_path / "summary.csv"
    report_path = tmp_path / "report.md"
    write_rows_csv(csv_path, rows)
    write_summary_csv(summary_csv_path, summary)
    write_research_report(report_path, config, rows, summary)
    assert "seed,extinct" in csv_path.read_text()
    assert "metric,mean,standard_deviation" in summary_csv_path.read_text()
    report = report_path.read_text()
    assert "Trials: 1" in report
    assert '"seed": 3' in report
    assert "seeds" in report


def test_visualizations_return_figures_and_can_be_saved(tmp_path: Path):
    species, gene = _trees()
    reconciliation = reconcile(gene, species, {"A": "A", "B": "B"})
    figures = [
        plot_species_tree(species),
        plot_gene_tree(gene),
        plot_reconciliation(gene, species, reconciliation),
    ]
    for index, figure in enumerate(figures):
        output = tmp_path / f"figure-{index}.png"
        write_figure(figure, output)
        assert output.stat().st_size > 0


def test_metric_chart_uses_summary_provenance(tmp_path: Path):
    config = ExperimentConfig("(A:1,B:1)S;", (5, 6), 0.0, 0.0)
    _, summary = run_experiment(config)
    figure = plot_metric_summary(summary, "inferred_cost")
    assert figure._phylo_bench_config == config
    assert figure._phylo_bench_seeds == (5, 6)
    output = tmp_path / "metric.png"
    write_figure(figure, output)
    assert output.stat().st_size > 0
    from PIL import Image
    with Image.open(output) as image:
        assert image.info["PhyloBenchSeeds"] == "[5, 6]"
        assert '"species_tree_newick": "(A:1,B:1)S;"' in image.info["PhyloBenchConfig"]


def test_reconciliation_mapping_uses_cross_axes_connections():
    from matplotlib.patches import ConnectionPatch
    species, gene = _trees()
    reconciliation = reconcile(gene, species, {"A": "A", "B": "B"})
    figure = plot_reconciliation(gene, species, reconciliation)
    connections = [
        artist for artist in figure.axes[0].patches
        if isinstance(artist, ConnectionPatch)
    ]
    assert len(connections) == len(reconciliation.node_mapping)
