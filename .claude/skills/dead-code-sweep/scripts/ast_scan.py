"""选项 1：常规静态分析（快，秒级）。AST 引用图 + 调用图闭包判不可达。

与变异法的分工：静态法抓"没人调用"，变异法抓"有人调用但结果没用"。静态法会漏掉
动态分派/注册表里的成员（那正是变异法的强项），所以两者的差集最有信息量。

用法: SWEEP_ROOT=tmp/deadcode python ast_scan.py [--json out.json]
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from collections import defaultdict

import sweepcfg


def refs_of(node: ast.AST) -> set[str]:
    """节点内引用到的所有名字：Name、属性末段、import 别名。

    注意 `ast.iter_child_nodes` 是浅遍历，必须 ast.walk 递归下探，否则嵌套在
    Expr/Assign/Call 里的引用全漏。
    """
    out: set[str] = set()
    for n in ast.walk(node):
        if isinstance(n, ast.Name):
            out.add(n.id)
        elif isinstance(n, ast.Attribute):
            out.add(n.attr)
        elif isinstance(n, ast.alias):          # from X import y -> y 在 alias.name，不是 Name
            out.add((n.asname or n.name).split(".")[0])
        elif isinstance(n, ast.Constant) and isinstance(n.value, str):
            out.add(n.value)                    # 注册表/契约里的字符串名
    return out


def module_refs(tree: ast.Module) -> set[str]:
    """模块顶层（含类体）的引用。`_HANDLERS = [FnHandler(topo.match_halo)]` 这类
    顶层赋值引用的 def 是活的，只在 def body 里找会误判死。"""
    out: set[str] = set()
    stack = list(tree.body)
    while stack:
        n = stack.pop()
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue                                # 函数体不算模块顶层
        if isinstance(n, ast.ClassDef):
            stack.extend(n.body)                    # 类体是顶层执行
            continue
        out |= refs_of(n)
    return out


def build(public_as_root: bool = True):
    funcs: dict[str, ast.AST] = {}                  # qual key -> node
    mod_refs: dict[str, set[str]] = {}
    func_refs: dict[str, set[str]] = {}
    decorated: set[str] = set()

    for p in sorted((sweepcfg.SRC / sweepcfg.PKG).rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        rel = p.relative_to(sweepcfg.ROOT).as_posix()
        tree = ast.parse(p.read_bytes().decode("utf-8"))
        mod_refs[rel] = module_refs(tree)

        def walk(node: ast.AST, prefix: str) -> None:
            for child in ast.iter_child_nodes(node):
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    qual = f"{prefix}{child.name}"
                    key = f"{rel}::{qual}"
                    funcs[key] = child
                    func_refs[key] = refs_of(child)
                    if child.decorator_list:
                        decorated.add(key)
                    walk(child, qual + ".")
                elif isinstance(child, ast.ClassDef):
                    walk(child, f"{prefix}{child.name}.")
                else:
                    walk(child, prefix)

        walk(tree, "")

    # 名字 -> 函数 key（同名跨模块全部登记，过近似是安全方向：宁可漏判死，不可误删活）
    by_name: dict[str, list[str]] = defaultdict(list)
    for key in funcs:
        by_name[key.split("::", 1)[1].split(".")[-1]].append(key)

    downstream_text: dict[str, str] = {}
    for d in sweepcfg.DOWNSTREAM:
        base = sweepcfg.REPO / d
        if not base.is_dir():
            continue
        it = base.glob("*.py") if d == "tools" else base.rglob("*.py")
        for p in it:
            if "__pycache__" in p.parts:
                continue
            rel = p.relative_to(sweepcfg.REPO).as_posix()
            if any(rel.startswith(x.rstrip("/") + "/") for x in sweepcfg.DOWNSTREAM_EXCLUDE):
                continue
            downstream_text[rel] = p.read_text(encoding="utf-8", errors="replace")

    roots: dict[str, str] = {}
    for key in funcs:
        name = key.split("::", 1)[1].split(".")[-1]
        if public_as_root and not name.startswith("_"):
            roots[key] = "公开 def（可能被包外调用）"
        if key in decorated:
            roots[key] = "带装饰器（@_register 等 = 注册即使用）"

    for key, node in funcs.items():
        name = key.split("::", 1)[1].split(".")[-1]
        for rel, refs in mod_refs.items():
            if rel != key.split("::", 1)[0] and name in refs:
                roots.setdefault(key, "模块顶层表达式引用（注册表/probe 元组）")
        for rel, text in downstream_text.items():
            if re.search(r"\b" + re.escape(name) + r"\b", text):
                roots.setdefault(key, f"下游引用（{rel} 等）")
                break

    # 建边：引用按名字解析到同名函数（过近似）
    edges: dict[str, set[str]] = defaultdict(set)
    for key, refs in func_refs.items():
        for r in refs:
            for tgt in by_name.get(r, ()):
                if tgt != key:
                    edges[key].add(tgt)
    for key, refs in mod_refs.items():
        for r in refs:
            for tgt in by_name.get(r, ()):
                if tgt.split("::", 1)[0] == key:
                    edges["::MODULE::" + key].add(tgt)

    reachable: set[str] = set()
    stack = [k for k in roots]
    while stack:
        cur = stack.pop()
        if cur in reachable:
            continue
        reachable.add(cur)
        stack.extend(edges.get(cur, ()))
    # 从模块顶层入口出发
    stack = [k for k in edges if k.startswith("::MODULE::")]
    while stack:
        cur = stack.pop()
        for tgt in edges.get(cur, ()):
            if tgt not in reachable:
                reachable.add(tgt)
                stack.append(tgt)

    unreachable = sorted(set(funcs) - reachable)
    return funcs, roots, unreachable, by_name


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    loose_funcs, roots, loose, _ = build(public_as_root=True)
    _, _, strict, _ = build(public_as_root=False)
    pre = f"{sweepcfg.SRC.name}/{sweepcfg.PKG}/"
    print(f"函数总数 {len(loose_funcs)}\n")
    print(f"## 宽松档（公开 def 也算根）—— 不可达 {len(loose)} 个")
    print("   保守：非下划线开头的函数一律假定可能被包外调用。适合先摸清盘子。\n")
    for key in loose:
        print(f"  {key.replace(pre, '')}")
    extra = [k for k in strict if k not in set(loose)]
    print(f"\n## 严格档（公开 def 不算根）—— 不可达 {len(strict)} 个（比宽松档多 {len(extra)}）")
    print("   只认下游引用/装饰器/模块顶层表达式为根。更激进，但也更容易误伤包外调用。\n")
    for key in extra:
        print(f"  + {key.replace(pre, '')}")
    print("\n提示：不可达 ≠ 可以删。注册表/probe 元组/契约测试里的成员会被本工具漏判，"
          "\n      删除前用变异法（选项 2）复核，或至少全量 benchmark + pytest 对比。")
    if args.json:
        sweepcfg.ROOT.joinpath(args.json).write_text(
            json.dumps({"loose": loose, "strict": strict, "roots": roots}, ensure_ascii=False, indent=1),
            encoding="utf-8")


if __name__ == "__main__":
    main()
