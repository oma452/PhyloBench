"""Standard-library CSV and figure export helpers."""

from __future__ import annotations

import csv
import io
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

from phylo_bench.io import atomic_write_text


def write_rows_csv(path: str | Path, rows: Iterable[Any]) -> None:
    """Write dataclass rows to CSV using their declared fields."""

    rows = list(rows)
    if not rows:
        raise ValueError("cannot export an empty row collection")
    records = [asdict(row) for row in rows]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
    atomic_write_text(path, buffer.getvalue())


def write_summary_csv(path: str | Path, summary: object) -> None:
    """Write a summary object's measured mean and standard deviation."""

    metrics = getattr(summary, "metrics", None)
    if not isinstance(metrics, dict):
        raise ValueError("summary must expose a metrics dictionary")
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    writer.writerow(["metric", "mean", "standard_deviation"])
    for name, (average, deviation) in metrics.items():
        writer.writerow([name, average, deviation])
    atomic_write_text(path, buffer.getvalue())


def write_figure(figure, path: str | Path, *, dpi: int = 180) -> None:
    """Save a Matplotlib figure in a user-selected format."""

    metadata = {}
    config = getattr(figure, "_phylo_bench_config", None)
    seeds = getattr(figure, "_phylo_bench_seeds", None)
    if config is not None:
        metadata["PhyloBenchConfig"] = json.dumps(asdict(config), sort_keys=True)
    if seeds is not None:
        metadata["PhyloBenchSeeds"] = json.dumps(list(seeds))
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.stem}.figure.tmp{output.suffix}")
    try:
        figure.savefig(
            temporary,
            dpi=dpi,
            bbox_inches="tight",
            metadata=metadata or None,
        )
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
