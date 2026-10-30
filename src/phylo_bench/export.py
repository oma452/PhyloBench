"""Standard-library CSV and figure export helpers."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Iterable


def write_rows_csv(path: str | Path, rows: Iterable[object]) -> None:
    """Write dataclass rows to CSV using their declared fields."""

    rows = list(rows)
    if not rows:
        raise ValueError("cannot export an empty row collection")
    records = [asdict(row) for row in rows]
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def write_summary_csv(path: str | Path, summary: object) -> None:
    """Write a summary object's measured mean and standard deviation."""

    metrics = getattr(summary, "metrics", None)
    if not isinstance(metrics, dict):
        raise ValueError("summary must expose a metrics dictionary")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "mean", "standard_deviation"])
        for name, (average, deviation) in metrics.items():
            writer.writerow([name, average, deviation])


def write_figure(figure, path: str | Path, *, dpi: int = 180) -> None:
    """Save a Matplotlib figure in a user-selected format."""

    metadata = {}
    config = getattr(figure, "_phylo_bench_config", None)
    seeds = getattr(figure, "_phylo_bench_seeds", None)
    if config is not None:
        metadata["PhyloBenchConfig"] = json.dumps(asdict(config), sort_keys=True)
    if seeds is not None:
        metadata["PhyloBenchSeeds"] = json.dumps(list(seeds))
    figure.savefig(path, dpi=dpi, bbox_inches="tight", metadata=metadata or None)
