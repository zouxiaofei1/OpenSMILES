"""函数清单提取 + 源码级 `return None` 变异。

变异用文本插入（不动其余字节），因为：
- read_text/write_text 在 Windows 会把 LF 翻成 CRLF，还原后与原文件不再逐字节相同；
- ast.unparse 会丢注释，而这个项目大量依赖注释承载规则出处。
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sweepcfg

PKG_DIR = sweepcfg.SRC / sweepcfg.PKG


def _qual_walk(node: ast.AST, prefix: str):
    for child in ast.iter_child_nodes(node):
        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
            qual = f"{prefix}{child.name}"
            yield child, qual
            yield from _qual_walk(child, qual + ".")
        elif isinstance(child, ast.ClassDef):
            yield from _qual_walk(child, f"{prefix}{child.name}.")
        else:
            yield from _qual_walk(child, prefix)


def _has_own_yield(node: ast.AST) -> bool:
    """该函数自身作用域（不含嵌套函数）里是否有 yield——决定能否 `return None`。"""
    stack = list(node.body)
    while stack:
        n = stack.pop()
        if isinstance(n, (ast.Yield, ast.YieldFrom)):
            return True
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            continue
        stack.extend(ast.iter_child_nodes(n))
    return False


def _rel(p: Path) -> str:
    return p.relative_to(sweepcfg.ROOT).as_posix()


def collect() -> list[dict]:
    """列出全部函数。key 形如 `src/namepredict/layer1/analyzer.py::ClassName.method`。"""
    seen_keys: dict[str, int] = {}
    out: list[dict] = []
    for p in sorted(PKG_DIR.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        try:
            tree = ast.parse(p.read_bytes().decode("utf-8"))
        except SyntaxError:
            continue
        for node, qual in _qual_walk(tree, ""):
            rec = {
                "file": _rel(p),
                "qual": qual,
                "key": f"{_rel(p)}::{qual}",
                "lineno": node.lineno,
                "col": node.col_offset,
                "async": isinstance(node, ast.AsyncFunctionDef),
                "async_gen": isinstance(node, ast.AsyncFunctionDef) and _has_own_yield(node),
            }
            if rec["key"] in seen_keys:  # 同名重定义：保留最后一个（后者生效），标重复
                rec["dup"] = True
            seen_keys[rec["key"]] = 1
            out.append(rec)
    out.sort(key=lambda d: (d["file"], d["lineno"]))
    return out


def _path(rec: dict) -> Path:
    return sweepcfg.ROOT / rec["file"]


def _insert_offset(text: str, first: ast.stmt) -> int:
    """函数体首语句在 text 里的字符偏移。col_offset 是 UTF-8 字节偏移，但语句起始前只有
    ASCII 缩进空白，故可当字符偏移直接用。"""
    return sum(len(x) for x in text.splitlines(keepends=True)[: first.lineno - 1]) + first.col_offset


def mutate(text: str, rec: dict, stmt: str | None = None) -> str | None:
    """在 rec 所指函数体首语句前插入 stmt；无法定位返回 None。"""
    tree = ast.parse(text)
    node = None
    for n in ast.walk(tree):
        if (isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
                and n.lineno == rec["lineno"] and n.col_offset == rec["col"]
                and n.name == rec["qual"].rsplit(".", 1)[-1]):
            node = n
            break
    if node is None or not node.body:
        return None
    first = node.body[0]
    body = stmt if stmt is not None else ("return" if rec.get("async_gen") else "return None")
    off = _insert_offset(text, first)
    if first.lineno == node.lineno:  # 单行函数体：分号衔接
        return text[:off] + body + "; " + text[off:]
    return text[:off] + body + "\n" + " " * first.col_offset + text[off:]


def apply(rec: dict, stmt: str | None = None) -> bytes | None:
    """写入变异；返回原文件原始字节，供 restore 逐字节还原。"""
    p = _path(rec)
    raw = p.read_bytes()
    new = mutate(raw.decode("utf-8"), rec, stmt)
    if new is None:
        return None
    compile(new, str(p), "exec")  # 语法自检，失败则抛出
    p.write_bytes(new.encode("utf-8"))
    return raw


def restore(rec: dict, raw: bytes) -> None:
    _path(rec).write_bytes(raw)


def snapshot() -> dict[str, str]:
    """全包 .py 的 md5，用于确认还原彻底。"""
    out = {}
    for p in PKG_DIR.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        out[_rel(p)] = hashlib.md5(p.read_bytes()).hexdigest()
    return out


def verify_clean(before: dict[str, str]) -> list[str]:
    after = snapshot()
    return sorted(set(before) ^ set(after)) + sorted(k for k in before if after.get(k) != before[k])
