"""Shared dependencies: namer only."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from namepredict.namer import SMILESNNamer
from namepredict.types import NameResult


def get_namer() -> SMILESNNamer:
    """Return a fresh SMILESNNamer (no stale process-wide singleton)."""
    return SMILESNNamer()


def name_result_dict(result: NameResult) -> dict[str, Any]:
    """Serialize NameResult for JSON responses."""
    return asdict(result)
