"""工作树配置。所有脚本共用：读 $SWEEP_ROOT/sweep.json（缺省用下面的默认值）。

工作树由 setup_tree.py 建立，是被分析代码的独立副本，避免动到主仓库
（项目有自动 commit 进程，就地改写源码有被提交的风险）。
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("SWEEP_ROOT") or ".").resolve()
_CFG_PATH = ROOT / "sweep.json"
CFG: dict = json.loads(_CFG_PATH.read_text(encoding="utf-8")) if _CFG_PATH.is_file() else {}

REPO = Path(CFG.get("repo") or ROOT).resolve()          # 主仓库根（扫下游引用用）
SRC = ROOT / CFG.get("src", "src")                      # sys.path 里放的那个目录
PKG = CFG.get("pkg", "opensmiles")                     # 被扫描的包名
IMPORT_ROOT = f"{SRC.name}.{PKG}"                       # 形如 src.opensmiles，用于相对路径
ENTRY = CFG.get("entry", "opensmiles.namer")           # 预导入入口（破坏它 = 导入期必需）
FACTORY = CFG.get("factory", "opensmiles.namer:SMILESNNamer")
BENCH_MODULE = CFG.get("bench_module", "benchmarks.benchmark")
BENCH_DATA = ROOT / CFG.get("bench_data", "benchmarks/merged_benchmark.json")
BENCH_CMD: list[str] = CFG.get(
    "bench_cmd",
    ["-m", "benchmarks.benchmark_parallel", "--data", "benchmarks/merged_benchmark.json"],
)
WORKERS = int(CFG.get("workers") or 0) or max(1, (os.cpu_count() or 2) - 1)
# harness 兼容垫片等：只服务于评测脚手架、不属于被测代码，必须排除在候选之外
SHIMS: set[str] = set(CFG.get("shims", []))
# 下游目录（判断"删了会不会打断调用方"时扫这些目录）
DOWNSTREAM: list[str] = CFG.get("downstream", ["tests", "benchmarks", "server", "tools"])
# 独立引擎副本所在目录，不是本包的下游，判定时排除
DOWNSTREAM_EXCLUDE: list[str] = CFG.get("downstream_exclude", ["tools/namepredict-v2", "tools/namepredict-v3"])

BASE_ROWS = ROOT / "_base_rows.json"      # 逐行预测基线
PROBE_ORDER = ROOT / "_probe_order.json"  # 探针行序（按覆盖度贪心）
COVERAGE = ROOT / "_coverage.json"
COV_ROWS = ROOT / "_cov_rows.json"
RESULTS = ROOT / "_mut_results.json"
ANALYSIS = ROOT / "_analysis.json"
ONESHOT = ROOT / "_oneshot.json"          # 只跑指定 key（调试用）
RUN_TMP = ROOT / "_mutrun_tmp.json"


def setup_path() -> None:
    for p in (ROOT, SRC):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def import_factory():
    """按配置导入 (smiles) -> 结果对象 的工厂类。"""
    setup_path()
    mod_name, _, cls = FACTORY.partition(":")
    mod = __import__(mod_name, fromlist=[cls or "x"])
    return getattr(mod, cls) if cls else mod
