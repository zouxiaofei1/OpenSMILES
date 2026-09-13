"""Shared dependencies: namer only."""

from __future__ import annotations

from dataclasses import asdict
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from namepredict.namer import SMILESNNamer
    from namepredict.types import NameResult


def get_namer() -> "SMILESNNamer":
    """Return a fresh SMILESNNamer (no stale process-wide singleton).

    namepredict 延迟到调用时导入: src/ 编译失败只让本请求 500, 不会连 app 一起崩。
    """
    from namepredict.namer import SMILESNNamer

    return SMILESNNamer()


def name_result_dict(result: NameResult) -> dict[str, Any]:
    """Serialize NameResult for JSON responses."""
    return asdict(result)
