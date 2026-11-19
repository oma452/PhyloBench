"""Setuptools build hooks for embedding repository provenance."""

from __future__ import annotations

import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


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


class build_py(_build_py):
    def run(self) -> None:
        super().run()
        root = Path(__file__).resolve().parent
        revision = _git_output(root, "rev-parse", "HEAD")
        dirty_output = _git_output(root, "status", "--porcelain")
        dirty = bool(dirty_output) if dirty_output is not None else None

        target = Path(self.build_lib) / "phylo_bench" / "_build_info.py"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            "'''Build metadata embedded for installed distributions.'''\n\n"
            "BUILD_INFO_EMBEDDED = True\n"
            "EMBEDDED_REVISION: str | None = "
            f"{revision!r}\n"
            "EMBEDDED_DIRTY: bool | None = "
            f"{dirty!r}\n",
            encoding="utf-8",
        )


setup(cmdclass={"build_py": build_py})
