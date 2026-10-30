"""Charts built from measured experiment summaries."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from phylo_bench.experiments.runner import ExperimentSummary


def plot_metric_summary(summary: ExperimentSummary, metric: str) -> Figure:
    """Plot a metric directly from a measured experiment summary."""

    if metric not in summary.metrics:
        raise ValueError(f"summary does not contain metric {metric!r}")
    figure, ax = plt.subplots(figsize=(7, 4))
    average, deviation = summary.metrics[metric]
    if average is None or deviation is None:
        ax.text(0.5, 0.5, "No valid observations", ha="center", va="center")
        ax.set_axis_off()
    else:
        ax.errorbar([metric], [average], yerr=[deviation], fmt="o", capsize=4)
        ax.set_ylabel(metric)
        ax.grid(axis="y", alpha=0.25)
    figure.suptitle(
        f"Measured {metric} ({summary.trials} trials; "
        f"seeds={','.join(map(str, summary.seeds))})"
    )
    figure.tight_layout()
    figure._phylo_bench_config = summary.config  # type: ignore[attr-defined]
    figure._phylo_bench_seeds = summary.seeds  # type: ignore[attr-defined]
    return figure
