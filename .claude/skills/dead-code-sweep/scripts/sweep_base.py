"""记录逐行预测基线 + 生成探针行序。

基线必须用与变异跑分**同一套 harness**算（同数据、同判分、同 factory），否则
比较的不是同一个东西。跑两遍确认逐行一致，基线才可信。
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import sweepcfg

sys.path.insert(0, str(sweepcfg.ROOT))


def _init_worker() -> None:
    sys.stdout = open(os.devnull, "w", encoding="utf-8")  # 管线里有遗留 debug print
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    globals()["_NAMER"] = sweepcfg.import_factory()()


def name_row(item):
    """(行号, smiles) -> (行号, en, zh)。异常编码进 zh，便于与基线逐字比对。"""
    idx, smiles = item
    try:
        _NAMER.cache.clear()  # 跨行状态隔离；管线依赖它保证同名结果不漂移
        r = _NAMER.name(smiles)
        return idx, r.en or "", r.zh or ""
    except Exception as exc:
        return idx, "", f"\x00EXC:{type(exc).__name__}"


def load_rows() -> list[dict]:
    with sweepcfg.BENCH_DATA.open(encoding="utf-8") as f:
        rows = json.load(f)
    if not isinstance(rows, list):
        raise ValueError(f"期望 JSON 列表: {sweepcfg.BENCH_DATA}")
    return rows


def probe_order(rows: list[dict], cov_rows: list[list[str]] | None = None) -> list[int]:
    """探针行序。有逐行覆盖数据时按"贪心最大新增覆盖"排——前若干行即踏遍几乎全部
    被执行函数，早退因此极易命中；否则按 (source, tier) 桶轮转。"""
    n = len(rows)
    if cov_rows:
        cov = [set(x) for x in cov_rows]
        remaining = set().union(*cov)
        order, used = [], set()
        while remaining and len(order) < n:
            best = max((i for i in range(n) if i not in used), key=lambda i: len(cov[i] & remaining))
            gain = cov[best] & remaining
            if not gain:
                break
            order.append(best)
            used.add(best)
            remaining -= gain
        return order + [i for i in range(n) if i not in used]
    buckets: dict[tuple, list[int]] = defaultdict(list)
    for i, r in enumerate(rows):
        buckets[(str(r.get("source") or ""), r.get("tier", 0))].append(i)
    order, keys = [], sorted(buckets)
    while any(buckets[k] for k in keys):
        for k in keys:
            if buckets[k]:
                order.append(buckets[k].pop(0))
    return order


def run(rows: list[dict]) -> list[list[str]]:
    payload = [(i, str(r.get("smiles") or "")) for i, r in enumerate(rows)]
    out = []
    with Pool(processes=sweepcfg.WORKERS, initializer=_init_worker) as pool:
        for res in pool.imap_unordered(name_row, payload, chunksize=16):
            out.append(res)
    out.sort()
    return [[en, zh] for _, en, zh in out]


def main() -> None:
    rows = load_rows()
    base = run(rows)
    again = run(rows)
    if base != again:
        bad = [i for i, (a, b) in enumerate(zip(base, again)) if a != b]
        print(f"⚠ 基线不确定：{len(bad)} 行两次结果不同，前几行 {bad[:10]}。先查清再继续。")
    sweepcfg.BASE_ROWS.write_text(json.dumps(base, ensure_ascii=False), encoding="utf-8")
    cov_rows = None
    if sweepcfg.COV_ROWS.is_file():
        cov_rows = json.loads(sweepcfg.COV_ROWS.read_text(encoding="utf-8"))
    sweepcfg.PROBE_ORDER.write_text(json.dumps(probe_order(rows, cov_rows)), encoding="utf-8")
    n_exc = sum(1 for _, zh in base if zh.startswith("\x00EXC:"))
    n_full = sum(1 for en, zh in base if en or zh)
    print(f"行数 {len(rows)} | 有输出 {n_full} | 基线即异常 {n_exc} -> {sweepcfg.BASE_ROWS.name}")


if __name__ == "__main__":
    main()
