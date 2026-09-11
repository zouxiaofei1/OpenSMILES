"""IUPAC 中文翻译 API: browse docs/iupac/cn_translated (markdown).

GET /api/v1/iupac        -> flat list of translated sections (name + parsed title, 目录.md first)
GET /api/v1/iupac/doc?path=P-66_中文翻译.md -> raw markdown; relative image links are rewritten
                         to /iupac-cn/<path> (figures are served by a static mount in app.py)

Path is validated to stay inside docs/iupac/cn_translated and to end in .md (no traversal).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from fastapi import APIRouter, Query

router = APIRouter(prefix="/api/v1", tags=["iupac"])

IUPAC_DIR = Path(__file__).resolve().parents[2] / "docs" / "iupac" / "cn_translated"

# ![](fig_p101/cholestane.png) -> ![](/iupac-cn/fig_p101/cholestane.png)
_IMG_RE = re.compile(r"(!\[[^\]]*\]\()([^)\s'\"#)]+)")


def _first_title(text: str) -> str:
    """Parse the leading '# ' heading of a translation file, fall back to ''."""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("# "):
            return line[2:].strip()
    return ""


@router.get("/iupac")
def iupac_index() -> dict[str, Any]:
    entries = []
    for p in sorted(IUPAC_DIR.iterdir(), key=lambda f: (f.name != "目录.md", f.name.lower())):
        if p.is_file() and p.suffix == ".md":
            title = ""
            try:
                title = _first_title(p.read_text(encoding="utf-8", errors="replace"))
            except OSError:
                title = ""
            entries.append({
                "name": p.stem,
                "path": p.relative_to(IUPAC_DIR).as_posix(),
                "title": title,
            })
    return {"ok": True, "entries": entries}


@router.get("/iupac/doc")
def iupac_doc(path: str = Query(..., min_length=1)) -> dict[str, Any]:
    resolved = (IUPAC_DIR / path).resolve()
    base = IUPAC_DIR.resolve()
    if resolved.suffix != ".md" or not str(resolved).startswith(str(base)):
        return {"ok": False, "error": "invalid path"}
    if not resolved.is_file():
        return {"ok": False, "error": "not found"}
    try:
        text = resolved.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return {"ok": False, "error": "unreadable file"}

    # Rewrite relative image links so they resolve against the static mount in app.py.
    def _rewrite_img(m: re.Match) -> str:
        prefix, path = m.group(1), m.group(2)
        if path.startswith("/") or path.startswith(("http://", "https://", "data:")):
            return m.group(0)
        return prefix + "/iupac-cn/" + path

    text = _IMG_RE.sub(_rewrite_img, text)
    return {
        "ok": True,
        "path": resolved.relative_to(IUPAC_DIR).as_posix(),
        "name": resolved.stem,
        "title": _first_title(text),
        "content": text,
    }
