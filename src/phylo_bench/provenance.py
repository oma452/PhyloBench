"""Run provenance and artifact checksum helpers."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from phylo_bench.experiments.runner import ExperimentConfig
from phylo_bench.io import atomic_write_text


_RUNTIME_DISTRIBUTIONS = ("phylo-bench", "biopython", "matplotlib")


def _git_provenance() -> dict[str, str | bool | None]:
    package_root = Path(__file__).resolve().parents[2]
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=package_root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=package_root,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        return {"revision": None, "working_tree_dirty": None}
    return {"revision": revision, "working_tree_dirty": dirty}


def _dependency_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for distribution in _RUNTIME_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def build_run_manifest(
    config: ExperimentConfig,
    *,
    command_line: Iterable[str],
    artifacts: dict[str, str] | None = None,
) -> dict[str, object]:
    """Build machine-readable provenance for one experiment execution."""

    return {
        "configuration": asdict(config),
        "command_line": list(command_line),
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
        },
        "dependencies": _dependency_versions(),
        "git": _git_provenance(),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "artifacts": artifacts or {},
    }


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(path: str | Path, manifest: dict[str, object]) -> None:
    """Write a manifest atomically as sorted, indented JSON."""

    atomic_write_text(
        path,
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
