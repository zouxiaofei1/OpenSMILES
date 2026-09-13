"""建立隔离工作树：把主仓库的代码副本拷进去，扫描全程只在副本上改写。

为什么要副本：主仓库有自动 commit 进程，就地改写源码有被提交的风险；且变异
扫描会反复改写/还原源码，任何其它并发读源码的进程都会读到掺杂状态。

用法:
  python setup_tree.py                    # 目标 tmp/deadcode
  python setup_tree.py --dest tmp/x --src src --pkg namepredict
  python setup_tree.py --from-archive src.7z    # 从 7z 快照建树
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


NORM_NAMES = ("normalize_en", "normalize_zh", "nospace")


def _harness_shim(dest: Path, cfg: dict) -> list[str]:
    """旧快照把判分归一化函数放在 tools.re，新版挪到了 constants。评测脚手架顶层 import
    它们，为跑通需给脚手架加回退。返回被垫的符号名——它们只服务于评测、不属于被测代码，
    要记进 shims 从候选里剔除（注入 return None 会让判分恒等成立、分数虚高到 100%）。
    """
    const = dest / cfg["src"] / cfg["pkg"] / "constants.py"
    if not const.is_file():
        return []
    if all(f"def {n}(" in const.read_text(encoding="utf-8") for n in NORM_NAMES):
        return []
    pat = re.compile(r"^(?P<ind>[ \t]*)from " + re.escape(cfg["pkg"]) + r"\.constants import (?P<names>[^\n]+)$",
                     re.M)
    shimmed: list[str] = []
    for fn in ("benchmark.py", "preview_metrics.py"):
        p = dest / "benchmarks" / fn
        if not p.is_file():
            continue
        text = p.read_text(encoding="utf-8")

        def repl(m):
            ind, names = m.group("ind"), m.group("names")
            if not any(n.strip() in NORM_NAMES for n in names.split(",")):
                return m.group(0)
            shimmed.extend(n.strip() for n in names.split(",") if n.strip() in NORM_NAMES)
            return (f"{ind}try:  # 兼容垫片：旧快照把判分归一化放在 tools.re\n"
                    f"{ind}    from {cfg['pkg']}.constants import {names}\n"
                    f"{ind}except ImportError:\n"
                    f"{ind}    from {cfg['pkg']}.tools.re import {names}")

        p.write_text(pat.sub(repl, text), encoding="utf-8")
    return sorted(set(shimmed))


def _copy(src: Path, dst: Path) -> None:
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns(  # 快照带上会打印无关的 diff_vs_last
                            "__pycache__", "*.pyc", ".last_run.json", ".last_stereo.json"))
    else:
        shutil.copy2(src, dst)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", help="主仓库根")
    ap.add_argument("--dest", default="tmp/deadcode", help="工作树位置（须在 .gitignore 覆盖范围内）")
    ap.add_argument("--src", default="src")
    ap.add_argument("--pkg", default="namepredict")
    ap.add_argument("--from-archive", default=None, help="改为从该 7z 快照取 src（解压到工作树）")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    dest = (repo / args.dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)

    if args.from_archive:
        sevenzip = shutil.which("7z") or r"C:\Program Files\7-Zip\7z.exe"
        subprocess.run([sevenzip, "x", str(repo / args.from_archive), f"-o{dest}", "-y"], check=True,
                       stdout=subprocess.DEVNULL)
    else:
        _copy(repo / args.src, dest / args.src)

    for extra in ("benchmarks", "tests", "pyproject.toml"):
        if (repo / extra).exists() and not (dest / extra).exists():
            _copy(repo / extra, dest / extra)

    for pyc in dest.rglob("__pycache__"):  # 归档里带来的陈旧 .pyc 会按 mtime+size 命中错编译
        shutil.rmtree(pyc, ignore_errors=True)

    cfg_path = dest / "sweep.json"
    if not cfg_path.is_file():
        cfg = {
            "repo": str(repo),
            "src": args.src,
            "pkg": args.pkg,
            "entry": f"{args.pkg}.namer",
            "factory": f"{args.pkg}.namer:SMILESNNamer",
            "bench_module": "benchmarks.benchmark",
            "bench_data": "benchmarks/merged_benchmark.json",
            "bench_cmd": ["-m", "benchmarks.benchmark_parallel", "--data", "benchmarks/merged_benchmark.json"],
            "shims": [],
            "downstream": ["tests", "benchmarks", "server", "tools"],
            "downstream_exclude": [f"tools/{args.pkg}-v2", f"tools/{args.pkg}-v3"],
        }
        shimmed = _harness_shim(dest, cfg)
        if shimmed:
            cfg["shims"] = [f"{cfg['src']}/{cfg['pkg']}/tools/re.py::{n}" for n in shimmed]
            print(f"已给评测脚手架加兼容垫片: {', '.join(shimmed)}（已记入 shims，不参与判死）")
        cfg_path.write_text(json.dumps(cfg, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"工作树就绪: {dest}")
    print(f"下一步: SWEEP_ROOT={dest} python <scripts>/sweep_base.py   # 记录逐行基线")
    return 0


if __name__ == "__main__":
    sys.exit(main())
