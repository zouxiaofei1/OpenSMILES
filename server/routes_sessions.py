"""Session list/read routes over agent_loop/memory/sessions."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from server.deps import SESSIONS_DIR

router = APIRouter(prefix="/api/v1/sessions", tags=["sessions"])

_CYCLE_RE = re.compile(r"^cycle-(\d{4})\.json$")


def _sessions_dir() -> Path:
    return Path(SESSIONS_DIR)


def _iter_from_name(name: str) -> int | None:
    m = _CYCLE_RE.match(name)
    return int(m.group(1)) if m else None


def _list_session_files() -> list[tuple[int, Path]]:
    root = _sessions_dir()
    if not root.is_dir():
        return []
    items: list[tuple[int, Path]] = []
    for path in root.iterdir():
        n = _iter_from_name(path.name)
        if n is not None and path.is_file():
            items.append((n, path))
    items.sort(key=lambda t: t[0])
    return items


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise HTTPException(status_code=500, detail="session must be object")
    return data


def _session_summary(n: int, path: Path) -> dict[str, Any]:
    data = _read_json(path)
    return {
        "iter": data.get("iter", n),
        "status": data.get("status"),
        "decision": data.get("decision"),
        "path": str(path).replace("\\", "/"),
    }


@router.get("")
def list_sessions() -> list[dict[str, Any]]:
    """List cycle sessions as summary dicts."""
    return [_session_summary(n, path) for n, path in _list_session_files()]


@router.get("/{iter_n}")
def get_session(iter_n: int) -> dict[str, Any]:
    """Return one cycle-XXXX.json payload."""
    path = _sessions_dir() / f"cycle-{iter_n:04d}.json"
    if not path.is_file():
        raise HTTPException(status_code=404, detail="session not found")
    return _read_json(path)
