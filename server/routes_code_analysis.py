"""Code analysis API: LOC / file counts per pipeline layer (0–5).

GET /api/v1/code-analysis — per-layer code/comment/blank lines, file counts, share.

Lines are classified coarsely:
  blank    → empty after strip
  comment  → stripped line starts with "#"
  code     → everything else (incl. docstrings)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fastapi import APIRouter

router = APIRouter(prefix="/api/v1", tags=["code-analysis"])

ROOT = Path(__file__).resolve().parents[1]
LAYER_DIR = ROOT / "src" / "namepredict"
LAYERS = list(range(0, 6))


def _count_lines(text: str) -> tuple[int, int, int]:
    code = comment = blank = 0
    for line in text.splitlines():
        s = line.strip()
        if not s:
            blank += 1
        elif s.startswith("#"):
            comment += 1
        else:
            code += 1
    return code, comment, blank


def _stat_layer(layer: int) -> dict[str, Any]:
    d = LAYER_DIR / f"layer{layer}"
    files: list[dict[str, Any]] = []
    counts = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
    if d.is_dir():
        for py in sorted(d.glob("*.py")):
            try:
                code, comment, blank = _count_lines(py.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError):
                continue
            files.append({"name": py.name, "code": code, "comment": comment, "blank": blank})
            counts["file_count"] += 1
            counts["code"] += code
            counts["comment"] += comment
            counts["blank"] += blank
    counts["lines"] = counts["code"] + counts["comment"] + counts["blank"]
    # 文件按代码行数降序（供前端展开明细）
    files.sort(key=lambda f: f["code"], reverse=True)
    return {"layer": layer, "path": f"src/namepredict/layer{layer}", "files": files, **counts}


@router.get("/code-analysis")
def code_analysis() -> dict[str, Any]:
    """Per-layer code statistics + aggregate, sorted by layer order."""
    layers = [_stat_layer(i) for i in LAYERS]
    agg = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
    for l in layers:
        for k in ("file_count", "code", "comment", "blank"):
            agg[k] += l[k]
    agg["lines"] = agg["code"] + agg["comment"] + agg["blank"]
    # share = code-line proportion of this layer across the pipeline
    for l in layers:
        l["share"] = round(100.0 * l["code"] / agg["code"], 1) if agg["code"] else 0.0
    return {"ok": True, "total": agg, "layers": layers}
