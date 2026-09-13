"""终局确认：把候选函数**一次性全部**注入 return None，跑一次**无早退的全量**
benchmark，比对总分与逐行预测是否与基线一致。

单函数跑分只管「这一行变了没」；这里管「全删了总分会不会动」——是删除前最后一道闸。

用法:
  SWEEP_ROOT=tmp/deadcode python sweep_verify.py            # 用每次扫描判为 same 的候选
  SWEEP_ROOT=tmp/deadcode python sweep_verify.py --keys k.json
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import mutlib
import sweepcfg

sys.path.insert(0, str(sweepcfg.ROOT))


def _insertions(text: str, recs: list[dict]) -> list[tuple[int, str]]:
    """同一文件多个函数：按**原始**文本算插入点，一次性落盘（避免行号漂移）。"""
    nodes = {}
    for n in ast.walk(ast.parse(text)):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            nodes.setdefault((n.lineno, n.col_offset, n.name), n)
    out = []
    for rec in recs:
        node = nodes.get((rec["lineno"], rec["col"], rec["qual"].rsplit(".", 1)[-1]))
        if node is None or not node.body:
            return []
        first = node.body[0]
        stmt = "return" if rec.get("async_gen") else "return None"
        lines = text.splitlines(keepends=True)
        off = sum(len(x) for x in lines[: first.lineno - 1]) + first.col_offset
        ins = stmt + "; " if first.lineno == node.lineno else stmt + "\n" + " " * first.col_offset
        out.append((off, ins))
    return out


def _row_compare(rows: list[dict]) -> list[int]:
    """用与基线相同的 harness 重算一遍，返回与基线不同的行号。"""
    import sweep_base

    base = json.loads(sweepcfg.BASE_ROWS.read_text(encoding="utf-8"))
    payload = [(i, str(r.get("smiles") or "")) for i, r in enumerate(rows)]
    out = []
    with Pool(processes=sweepcfg.WORKERS, initializer=sweep_base._init_worker) as pool:
        for idx, en, zh in pool.imap_unordered(sweep_base.name_row, payload, chunksize=16):
            if [en, zh] != base[idx]:
                out.append(idx)
    return sorted(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", default=None)
    args = ap.parse_args()

    if args.keys:
        keys = json.loads(Path(args.keys).read_text(encoding="utf-8"))
    else:
        res = json.loads(sweepcfg.RESULTS.read_text(encoding="utf-8"))
        keys = [k for k, v in res.items() if v.get("status") == "same" and k not in sweepcfg.SHIMS]
    allrecs = {r["key"]: r for r in mutlib.collect()}
    by_file: dict[str, list[dict]] = defaultdict(list)
    for k in keys:
        if k in allrecs:
            by_file[allrecs[k]["file"]].append(allrecs[k])
    print(f"候选 {len(keys)} 个 -> {len(by_file)} 个文件")

    baseline_snapshot = mutlib.snapshot()  # 必须在注入之前取
    saved: dict[str, bytes] = {}
    applied = 0
    for f, recs in by_file.items():
        p = sweepcfg.ROOT / f
        raw = p.read_bytes()
        ins = _insertions(raw.decode("utf-8"), recs)
        if not ins:
            print("⚠ 插入点解析失败，跳过:", f)
            continue
        text = raw.decode("utf-8")
        for off, s in sorted(ins, reverse=True):
            text = text[:off] + s + text[off:]
        compile(text, str(p), "exec")
        saved[f] = raw
        p.write_bytes(text.encode("utf-8"))
        applied += len(ins)
    print(f"已对 {len(saved)} 个文件的 {applied} 个函数注入 return None")

    try:
        import sweep_base

        rows = sweep_base.load_rows()
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", SWEEP_ROOT=str(sweepcfg.ROOT),
                   PYTHONIOENCODING="utf-8")
        proc = subprocess.run([sys.executable, *sweepcfg.BENCH_CMD], cwd=str(sweepcfg.ROOT),
                              capture_output=True, text=True, timeout=3600, env=env,
                              encoding="utf-8", errors="replace")
        tail = [ln for ln in (proc.stdout or "").splitlines() if not ln.startswith("progress")]
        print("benchmark 退出码:", proc.returncode)
        print("\n".join(tail[-3:]))
        if not tail:
            print("stderr:", (proc.stderr or "")[-800:])
        diff = _row_compare(rows)
        print(f"逐行比对: {'完全一致 ✓' if not diff else f'⚠ {len(diff)} 行不同 {diff[:10]}'}")
    finally:
        for f, raw in saved.items():
            (sweepcfg.ROOT / f).write_bytes(raw)
    dirty = mutlib.verify_clean(baseline_snapshot)
    print("源码还原校验:", "干净 ✓" if not dirty else f"⚠ 残留 {dirty}")


if __name__ == "__main__":
    main()
