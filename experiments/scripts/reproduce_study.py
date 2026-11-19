#!/usr/bin/env python3
"""Regenerate new experiment artifacts from an archived study configuration.

This script intentionally writes to a caller-selected output directory. It
does not overwrite the archived CSV, figure, or report artifacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_METADATA = {
    "first": ROOT / "experiments/2026-10-30-first-study/metadata/study_configuration.json",
    "second": ROOT / "experiments/2026-10-30-second-study/metadata/study_configuration.json",
}


def _configs(metadata: dict[str, object], study: str) -> dict[str, dict[str, object]]:
    if study == "first":
        names = ("low", "medium", "high")
        return {name: metadata[f"{name}_config"] for name in names}  # type: ignore[return-value]
    return {
        name: value["config"]  # type: ignore[index]
        for name, value in metadata.items()
        if name.startswith("d") and isinstance(value, dict) and "config" in value
    }


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Regenerate PhyloBench experiment artifacts in a new directory."
    )
    parser.add_argument("--study", choices=("first", "second"), required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    metadata_path = args.metadata or DEFAULT_METADATA[args.study]
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    generated: dict[str, object] = {}

    for name, config in sorted(_configs(metadata, args.study).items()):
        cell_dir = output_dir / name
        cell_dir.mkdir(parents=True, exist_ok=True)
        config_path = cell_dir / "config.json"
        config_path.write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        command = [
            sys.executable,
            "-m",
            "phylo_bench.cli",
            "run",
            str(config_path),
            "--output-dir",
            str(cell_dir),
        ]
        subprocess.run(command, cwd=ROOT, check=True)
        generated[name] = {
            "configuration_sha256": _digest(config_path),
            "manifest": str((cell_dir / "run_manifest.json").relative_to(output_dir)),
        }

    study_manifest = {
        "study": args.study,
        "source_metadata": str(metadata_path.relative_to(ROOT)),
        "source_metadata_sha256": _digest(metadata_path),
        "generated_cells": generated,
        "note": (
            "This is an independently executable regeneration workflow. "
            "It does not overwrite or claim byte-for-byte identity with the "
            "archived study artifacts because their original generator scripts "
            "are not present."
        ),
    }
    (output_dir / "study_manifest.json").write_text(
        json.dumps(study_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
