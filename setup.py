"""Setuptools build hooks for embedding repository provenance."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from setuptools import setup  # noqa: E402
from setuptools.command.build_py import build_py as _build_py  # noqa: E402
from build_provenance import git_dirty, git_revision  # noqa: E402


class build_py(_build_py):
    def run(self) -> None:
        root = ROOT
        revision = git_revision(root)
        dirty = git_dirty(root)
        super().run()

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
