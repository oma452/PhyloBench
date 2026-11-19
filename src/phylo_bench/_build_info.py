"""Build metadata embedded by setuptools for wheel installations."""

from __future__ import annotations

BUILD_INFO_EMBEDDED = False
EMBEDDED_REVISION: str | None = None
EMBEDDED_DIRTY: bool | None = None
