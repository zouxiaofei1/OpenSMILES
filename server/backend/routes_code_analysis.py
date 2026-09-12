"""Code analysis API: LOC / file counts per pipeline layer (0–5).

GET /api/v1/code-analysis — per-layer code/comment/blank lines, file counts, share.

Lines are classified:
  blank    → empty after strip
  comment  → stripped line starts with "#", or falls inside a docstring
             (module / function / class triple-quoted block)
  code     → everything else
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from server.backend import history_store

router = APIRouter(prefix="/api/v1", tags=["code-analysis"])

ROOT = Path(__file__).resolve().parents[2]
LAYER_DIR = ROOT / "src" / "namepredict"
# (gid, dirname, label, recursive) — tools is the shared layer-agnostic package;
# the empty dirname is the package root itself (__init__/constants/namer/types),
# counted non-recursively so it does not swallow the layer dirs below it.
GROUPS = [
    (0, "layer0", "Layer 0", True), (1, "layer1", "Layer 1", True),
    (2, "layer2", "Layer 2", True), (3, "layer3", "Layer 3", True),
    (4, "layer4", "Layer 4", True), (5, "layer5", "Layer 5", True),
    (6, "tools", "Tools", True), (7, "", "Core", False),
]


def _docstring_lines(tree: ast.Module) -> set[int]:
    """收集所有 docstring（模块/函数/类首语句字符串）占用的行号。"""
    lines: set[int] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        body = node.body
        if not body or not isinstance(body[0], ast.Expr):
            continue
        value = body[0].value
        if isinstance(value, ast.Constant) and isinstance(value.value, str):
            end = getattr(value, "end_lineno", value.lineno)
            lines.update(range(value.lineno, end + 1))
    return lines


def _count_lines(text: str) -> tuple[int, int, int]:
    code = comment = blank = 0
    doc_lines: set[int] = set()
    try:
        doc_lines = _docstring_lines(ast.parse(text))
    except SyntaxError:
        pass  # 无法解析时退回逐行分类，docstring 不计入注释
    for lineno, line in enumerate(text.splitlines(), start=1):
        s = line.strip()
        if not s:
            blank += 1
        elif lineno in doc_lines or s.startswith("#"):
            comment += 1
        else:
            code += 1
    return code, comment, blank


def _stat_group(gid: int, dirname: str, label: str, recursive: bool) -> dict[str, Any]:
    """统计单个分组（层/工具/包根）的 code/comment/blank 行数与文件明细。"""
    d = LAYER_DIR / dirname
    glob = d.rglob if recursive else d.glob  # 非递归用于包根，避免吞掉下层目录
    files: list[dict[str, Any]] = []
    counts = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
    if d.is_dir():
        for py in sorted(glob("*.py")):  # 递归统计子目录（layer2、tools/leaves 等）
            if "__pycache__" in py.parts:
                continue
            try:
                code, comment, blank = _count_lines(py.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError):
                continue
            rel = py.relative_to(d).as_posix()
            files.append({"name": rel, "code": code, "comment": comment, "blank": blank})
            counts["file_count"] += 1
            counts["code"] += code
            counts["comment"] += comment
            counts["blank"] += blank
    counts["lines"] = counts["code"] + counts["comment"] + counts["blank"]
    # 文件按代码行数降序（供前端展开明细）
    files.sort(key=lambda f: f["code"], reverse=True)
    return {
        "layer": gid, "label": label, "path": f"src/namepredict/{dirname}".rstrip("/"),
        "files": files, **counts,
    }


def _code_analysis_live() -> dict[str, Any]:
    """Per-layer/tools code statistics + aggregate, sorted by group order."""
    layers = [_stat_group(*g) for g in GROUPS]
    agg = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
    for l in layers:
        for k in ("file_count", "code", "comment", "blank"):
            agg[k] += l[k]
    agg["lines"] = agg["code"] + agg["comment"] + agg["blank"]
    # share = code-line proportion of this layer across the pipeline
    for l in layers:
        l["share"] = round(100.0 * l["code"] / agg["code"], 1) if agg["code"] else 0.0
    return {"ok": True, "total": agg, "layers": layers}


def _code_analysis_history(commit: str) -> dict[str, Any]:
    """Code statistics for a past commit, read statically from git objects."""
    cache = history_store.read_cache(commit, "code-analysis")
    if cache and cache.get("_data_sig") == history_store.data_sig():
        return cache["payload"]

    names = history_store.list_tree(commit, "src/namepredict")
    layers: list[dict[str, Any]] = []
    for gid, dirname, label, recursive in GROUPS:
        prefix = f"src/namepredict/{dirname}/" if dirname else "src/namepredict/"
        counts = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
        files: list[dict[str, Any]] = []
        for p in names:
            if not p.startswith(prefix):
                continue
            rel = p[len(prefix):]
            if not recursive and "/" in rel:
                continue  # 包根组只收直属文件
            text = history_store.show_file(commit, p)
            if text is None:
                continue
            code, comment, blank = _count_lines(text)
            files.append({"name": rel, "code": code, "comment": comment, "blank": blank})
            counts["file_count"] += 1
            counts["code"] += code
            counts["comment"] += comment
            counts["blank"] += blank
        counts["lines"] = counts["code"] + counts["comment"] + counts["blank"]
        files.sort(key=lambda f: f["code"], reverse=True)
        layers.append({
            "layer": gid, "label": label, "path": f"src/namepredict/{dirname}".rstrip("/"),
            "files": files, **counts,
        })

    agg = {"file_count": 0, "code": 0, "comment": 0, "blank": 0}
    for l in layers:
        for k in ("file_count", "code", "comment", "blank"):
            agg[k] += l[k]
    agg["lines"] = agg["code"] + agg["comment"] + agg["blank"]
    for l in layers:
        l["share"] = round(100.0 * l["code"] / agg["code"], 1) if agg["code"] else 0.0

    payload = {"ok": True, "total": agg, "layers": layers}
    history_store.write_cache(commit, "code-analysis", {
        "_data_sig": history_store.data_sig(), "payload": payload,
    })
    return payload


@router.get("/code-analysis")
def code_analysis(commit: str | None = None) -> dict[str, Any]:
    """Per-layer code stats for the current code, or a past git commit."""
    if not commit:
        return _code_analysis_live()
    try:
        full = history_store.resolve_commit(commit)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"invalid commit: {commit}")
    if full is None:
        return _code_analysis_live()
    return _code_analysis_history(full)
