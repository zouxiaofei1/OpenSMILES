"""Wiki API: browse docs/wiki (markdown) with a tree navigation.

GET /api/v1/wiki        -> directory tree of docs/wiki (dirs + .md files)
GET /api/v1/wiki/doc?path=architecture/layer1-analyzer.md -> raw markdown content

Path is validated to stay inside docs/wiki and to end in .md (no traversal).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import APIRouter, Query

router = APIRouter(prefix="/api/v1", tags=["wiki"])

WIKI_DIR = Path(__file__).resolve().parents[1] / "docs" / "wiki"


def _walk(d: Path) -> dict[str, Any]:
    """Recursively build {name, path, dirs, files} for the tree."""
    entry = {
        "name": d.name,
        "path": d.relative_to(WIKI_DIR).as_posix(),
        "dirs": [],
        "files": [],
    }
    for p in sorted(d.iterdir()):
        if p.is_dir():
            entry["dirs"].append(_walk(p))
        elif p.suffix == ".md":
            entry["files"].append({
                "name": p.stem,
                "path": p.relative_to(WIKI_DIR).as_posix(),
            })
    return entry


@router.get("/wiki")
def wiki_tree() -> dict[str, Any]:
    return {"ok": True, "tree": _walk(WIKI_DIR)}


@router.get("/wiki/doc")
def wiki_doc(path: str = Query(..., min_length=1)) -> dict[str, Any]:
    resolved = (WIKI_DIR / path).resolve()
    base = WIKI_DIR.resolve()
    if resolved.suffix != ".md" or not str(resolved).startswith(str(base)):
        return {"ok": False, "error": "invalid path"}
    if not resolved.is_file():
        return {"ok": False, "error": "not found"}
    try:
        text = resolved.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return {"ok": False, "error": "unreadable file"}
    return {
        "ok": True,
        "path": resolved.relative_to(WIKI_DIR).as_posix(),
        "name": resolved.stem,
        "content": text,
    }
