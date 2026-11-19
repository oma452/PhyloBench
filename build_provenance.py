"""Helpers for recording source revision and source-tree dirty state."""

from __future__ import annotations

import subprocess
from pathlib import Path


def _git_output(root: Path, *arguments: str) -> str | None:
    try:
        return subprocess.run(
            ["git", *arguments],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _is_generated_path(path: str) -> bool:
    normalized = path.replace("\\", "/").strip()
    parts = normalized.split("/")
    if parts and parts[0] in {"build", "dist"}:
        return True
    return any(part.endswith(".egg-info") for part in parts)


def git_revision(root: Path) -> str | None:
    revision = _git_output(root, "rev-parse", "HEAD")
    return revision or None


def git_dirty(root: Path) -> bool | None:
    """Return whether tracked or untracked source files differ from HEAD."""

    staged = _git_output(root, "diff", "--cached", "--name-only", "HEAD", "--")
    unstaged = _git_output(root, "diff", "--name-only")
    untracked = _git_output(
        root,
        "ls-files",
        "--others",
        "--exclude-standard",
    )
    if staged is None or unstaged is None or untracked is None:
        return None

    changed_paths: list[str] = []
    for output in (staged, unstaged, untracked):
        changed_paths.extend(
            line.strip()
            for line in output.splitlines()
            if line.strip()
        )
    return any(not _is_generated_path(path) for path in changed_paths)
