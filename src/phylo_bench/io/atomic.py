"""Atomic filesystem writers for experiment artifacts."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path


def atomic_write_bytes(path: str | Path, data: bytes) -> None:
    """Write bytes and replace the destination only after success."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            temporary.write(data)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_path, destination)
    finally:
        temporary_path.unlink(missing_ok=True)


def atomic_write_text(
    path: str | Path,
    text: str,
    *,
    encoding: str = "utf-8",
) -> None:
    """Write text and replace the destination only after success."""

    atomic_write_bytes(path, text.encode(encoding))
