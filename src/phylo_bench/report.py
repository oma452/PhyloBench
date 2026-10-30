"""Markdown research reports generated only from measured experiment outputs."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from phylo_bench.experiments.runner import ExperimentConfig, ExperimentRow, ExperimentSummary


def write_research_report(
    path: str | Path,
    config: ExperimentConfig,
    rows: tuple[ExperimentRow, ...],
    summary: ExperimentSummary,
) -> None:
    """Write a reproducibility report with actual rows and aggregate values."""

    lines = [
        "# PhyloBench experimental report",
        "",
        "## Methodology",
        "",
        "Experiments use the rooted binary duplication-loss simulator and "
        "classical LCA reconciliation. Duplication and loss events are "
        "simulated with configurable continuous-time rates on species-tree "
        "branches. Extinct trials are reported separately.",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(asdict(config), indent=2),
        "```",
        "",
        "## Results",
        "",
        f"- Trials: {summary.trials}",
        f"- Extinct trials: {summary.extinct_trials}",
        "",
        "Inference metrics use only non-extinct measured trials. Complete "
        "history counts use every simulated trial; unavailable inference "
        "metrics are not replaced with zeros.",
        "",
        "Complete simulated loss counts describe the full hidden history. "
        "Identifiable loss metrics use only losses with an extant sibling "
        "witness; hidden losses are excluded from precision and recall.",
        "",
        "| Metric | Mean | Standard deviation |",
        "|---|---:|---:|",
    ]
    for name, (average, deviation) in summary.metrics.items():
        average_text = "unavailable" if average is None else f"{average:.6g}"
        deviation_text = "unavailable" if deviation is None else f"{deviation:.6g}"
        lines.append(f"| {name} | {average_text} | {deviation_text} |")
    lines.extend([
        "",
        "## Interpretation and limitations",
        "",
        "Exact event-location agreement applies to inferred speciation "
        "locations. Duplication branch compatibility is a weaker measure "
        "because exact within-branch duplication time is generally not "
        "identifiable from extant leaves. Complete simulated events and "
        "inferred extant-tree events are reported separately.",
        "",
        "## Trial records",
        "",
        "```json",
        json.dumps([asdict(row) for row in rows], indent=2),
        "```",
        "",
    ])
    Path(path).write_text("\n".join(lines), encoding="utf-8")
