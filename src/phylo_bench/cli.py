"""Command-line interface for reproducible PhyloBench experiments."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from phylo_bench.experiments import ExperimentConfig, run_experiment
from phylo_bench.export import write_figure, write_rows_csv, write_summary_csv
from phylo_bench.report import write_research_report
from phylo_bench.visualization import plot_metric_summary


_FIGURE_METRICS = (
    "complete_duplications",
    "complete_losses",
    "inferred_cost",
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phylo-bench",
        description="Run reproducible PhyloBench simulation experiments.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    for name, help_text in (
        ("run", "Run an experiment and export measured artifacts."),
        ("inspect", "Validate and print a JSON experiment configuration."),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument(
            "config",
            type=Path,
            help="Path to a JSON configuration file.",
        )
        if name == "run":
            command.add_argument(
                "-o",
                "--output-dir",
                type=Path,
                required=True,
                help="Directory for CSV, PNG, Markdown, and metadata outputs.",
            )
    return parser


def _load_config(path: Path) -> ExperimentConfig:
    try:
        with path.open(encoding="utf-8") as handle:
            data: Any = json.load(handle)
    except FileNotFoundError as exc:
        raise ValueError(f"configuration file does not exist: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"configuration file is not valid JSON: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError("configuration JSON must contain an object")
    allowed = {
        "species_tree_newick",
        "seeds",
        "duplication_rate",
        "loss_rate",
        "duplication_cost",
        "loss_cost",
        "max_events",
        "max_lineages",
    }
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"unknown configuration field(s): {', '.join(unknown)}")
    required = {"species_tree_newick", "seeds", "duplication_rate", "loss_rate"}
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"missing required configuration field(s): {', '.join(missing)}")
    try:
        return ExperimentConfig(**data)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid experiment configuration: {exc}") from exc


def _run(config: ExperimentConfig, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows, summary = run_experiment(config)
    write_rows_csv(output_dir / "trials.csv", rows)
    write_summary_csv(output_dir / "summary.csv", summary)

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    for metric in _FIGURE_METRICS:
        figure = plot_metric_summary(summary, metric)
        try:
            write_figure(figure, figures_dir / f"{metric}.png")
        finally:
            figure.clf()

    write_research_report(output_dir / "report.md", config, rows, summary)
    metadata = {
        "configuration": asdict(config),
        "summary": {
            "trials": summary.trials,
            "extinct_trials": summary.extinct_trials,
            "metrics": summary.metrics,
        },
        "artifacts": {
            "trials_csv": "trials.csv",
            "summary_csv": "summary.csv",
            "report": "report.md",
            "figures": [f"figures/{metric}.png" for metric in _FIGURE_METRICS],
        },
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"completed {summary.trials} trials ({summary.extinct_trials} extinct)")
    print(f"artifacts written to {output_dir}")


def main(argv: list[str] | None = None) -> int:
    """Run the CLI and return a process exit code."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        config = _load_config(args.config)
        if args.command == "inspect":
            print(json.dumps(asdict(config), indent=2, sort_keys=True))
        else:
            _run(config, args.output_dir)
    except ValueError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
