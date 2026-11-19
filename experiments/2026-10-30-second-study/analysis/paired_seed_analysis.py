#!/usr/bin/env python3
"""Compute descriptive paired contrasts for the common-seed factorial design."""

from __future__ import annotations

import argparse
import csv
import math
from itertools import combinations
from pathlib import Path
from statistics import mean, stdev


METRICS = (
    "complete_duplications",
    "complete_losses",
    "inferred_cost",
    "loss_branch_precision",
    "loss_branch_recall",
)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _cell_key(row: dict[str, str]) -> tuple[float, float]:
    return float(row["duplication_rate"]), float(row["loss_rate"])


def _contrast(
    left: list[dict[str, str]],
    right: list[dict[str, str]],
    metric: str,
) -> tuple[int, float | None, float | None, float | None]:
    left_by_seed = {int(row["seed"]): row for row in left}
    right_by_seed = {int(row["seed"]): row for row in right}
    differences = []
    for seed in sorted(left_by_seed.keys() & right_by_seed.keys()):
        left_value = left_by_seed[seed][metric]
        right_value = right_by_seed[seed][metric]
        if left_value == "" or right_value == "":
            continue
        differences.append(float(right_value) - float(left_value))
    if not differences:
        return 0, None, None, None
    average = mean(differences)
    sample_sd = stdev(differences) if len(differences) > 1 else 0.0
    return len(differences), average, sample_sd, sample_sd / math.sqrt(len(differences))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute common-seed paired contrasts without changing archived data."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "raw" / "all_trials.csv",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = _read_rows(args.input)
    cells: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in rows:
        cells.setdefault(_cell_key(row), []).append(row)

    records = []
    for (left_key, left_rows), (right_key, right_rows) in combinations(
        sorted(cells.items()), 2
    ):
        for metric in METRICS:
            count, average, sample_sd, mcse = _contrast(
                left_rows,
                right_rows,
                metric,
            )
            records.append(
                {
                    "left_duplication_rate": left_key[0],
                    "left_loss_rate": left_key[1],
                    "right_duplication_rate": right_key[0],
                    "right_loss_rate": right_key[1],
                    "metric": metric,
                    "paired_observations": count,
                    "right_minus_left_mean": average,
                    "paired_sample_sd": sample_sd,
                    "paired_mcse": mcse,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
