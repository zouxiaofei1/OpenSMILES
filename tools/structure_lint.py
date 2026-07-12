"""Static structure lint for namepredict layers and hard constraints."""

from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path

MAX_FILE_LINES = 500
MAX_FUNC_BODY = 10
MAX_CACHE_ENTRIES = 100
SMILES_EQ_RE = re.compile(r"if\s+smiles\s*==")
LAYER_RE = re.compile(r"layer([0-5])")
LAYER_IMPORT_RE = re.compile(
    r"(?:from\s+namepredict\.layer([0-5])|import\s+namepredict\.layer([0-5]))"
)


def _iter_py_files(root: Path) -> list[Path]:
    return sorted(p for p in root.rglob("*.py") if p.is_file())


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _check_file_lines(path: Path, root: Path) -> list[str]:
    n = len(path.read_text(encoding="utf-8").splitlines())
    if n > MAX_FILE_LINES:
        return [f"{_rel(path, root)}: file has {n} lines (max {MAX_FILE_LINES})"]
    return []


def _body_nonempty_lines(source: str, node: ast.AST) -> int:
    lines = source.splitlines()
    start = getattr(node, "lineno", 1)
    end = getattr(node, "end_lineno", start)
    body = lines[start:end]
    return sum(1 for ln in body if ln.strip())


def _func_nodes(tree: ast.AST) -> list[ast.AST]:
    out: list[ast.AST] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append(node)
    return out


def _func_issue(path: Path, root: Path, node: ast.AST, source: str) -> str | None:
    n = _body_nonempty_lines(source, node)
    if n <= MAX_FUNC_BODY:
        return None
    name = getattr(node, "name", "<fn>")
    line = getattr(node, "lineno", 0)
    return (
        f"{_rel(path, root)}:{line}: function {name!r} body has "
        f"{n} non-empty lines (max {MAX_FUNC_BODY})"
    )


def _check_func_bodies(path: Path, root: Path, source: str, tree: ast.AST) -> list[str]:
    out: list[str] = []
    for node in _func_nodes(tree):
        msg = _func_issue(path, root, node, source)
        if msg:
            out.append(msg)
    return out


def _check_smiles_eq(path: Path, root: Path, source: str) -> list[str]:
    issues: list[str] = []
    for i, line in enumerate(source.splitlines(), 1):
        if SMILES_EQ_RE.search(line):
            issues.append(f"{_rel(path, root)}:{i}: banned smiles equality special-case")
    return issues


def _layer_of(path: Path, root: Path) -> int | None:
    parts = path.relative_to(root).parts
    for part in parts:
        m = LAYER_RE.fullmatch(part)
        if m:
            return int(m.group(1))
    return None


def _imported_layers(source: str) -> list[int]:
    found: list[int] = []
    for m in LAYER_IMPORT_RE.finditer(source):
        g = m.group(1) or m.group(2)
        found.append(int(g))
    return found


def _check_layer_imports(path: Path, root: Path, source: str) -> list[str]:
    layer = _layer_of(path, root)
    if layer is None:
        return []
    issues: list[str] = []
    for other in _imported_layers(source):
        if other > layer:
            issues.append(
                f"{_rel(path, root)}: layer{layer} imports higher layer{other}"
            )
    return issues


def _dict_key_count(node: ast.Dict) -> int:
    return sum(1 for k in node.keys if k is not None)


def _module_level_dict_keys(tree: ast.AST) -> int:
    if not isinstance(tree, ast.Module):
        return 0
    total = 0
    for stmt in tree.body:
        total += _stmt_dict_keys(stmt)
    return total


def _stmt_dict_keys(stmt: ast.stmt) -> int:
    if isinstance(stmt, ast.Assign):
        return _value_dict_keys(stmt.value)
    if isinstance(stmt, ast.AnnAssign) and stmt.value is not None:
        return _value_dict_keys(stmt.value)
    if isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Dict):
        return _dict_key_count(stmt.value)
    return 0


def _value_dict_keys(value: ast.AST) -> int:
    if isinstance(value, ast.Dict):
        return _dict_key_count(value)
    if isinstance(value, ast.Call):
        return _call_dict_keys(value)
    return 0


def _call_dict_keys(call: ast.Call) -> int:
    total = 0
    for arg in call.args:
        if isinstance(arg, ast.Dict):
            total += _dict_key_count(arg)
    for kw in call.keywords:
        if isinstance(kw.value, ast.Dict):
            total += _dict_key_count(kw.value)
    return total


def _is_common_names(path: Path) -> bool:
    return path.name == "common_names.py"


def _check_cache_entries(path: Path, root: Path, tree: ast.AST) -> list[str]:
    if not _is_common_names(path):
        return []
    n = _module_level_dict_keys(tree)
    if n > MAX_CACHE_ENTRIES:
        return [
            f"{_rel(path, root)}: common_names preloaded entries "
            f"{n} (max {MAX_CACHE_ENTRIES})"
        ]
    return []


def _parse(path: Path) -> tuple[str, ast.AST | None, list[str]]:
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as exc:
        return source, None, [f"{path}: syntax error: {exc}"]
    return source, tree, []


def lint_file(path: Path, root: Path) -> list[str]:
    issues = _check_file_lines(path, root)
    source, tree, parse_issues = _parse(path)
    issues.extend(parse_issues)
    if tree is None:
        return issues
    issues.extend(_check_func_bodies(path, root, source, tree))
    issues.extend(_check_smiles_eq(path, root, source))
    issues.extend(_check_layer_imports(path, root, source))
    issues.extend(_check_cache_entries(path, root, tree))
    return issues


def lint_tree(root: Path) -> list[str]:
    issues: list[str] = []
    for path in _iter_py_files(root):
        issues.extend(lint_file(path, root))
    return issues


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Structure lint for namepredict")
    p.add_argument("--root", type=Path, required=True, help="Package root")
    return p


def _report(issues: list[str]) -> int:
    for msg in issues:
        print(msg)
    if issues:
        print(f"{len(issues)} issue(s)", file=sys.stderr)
        return 1
    print("structure_lint: ok")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        print(f"error: root not a directory: {root}", file=sys.stderr)
        return 1
    return _report(lint_tree(root))


if __name__ == "__main__":
    raise SystemExit(main())
