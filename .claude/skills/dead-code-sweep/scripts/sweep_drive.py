"""变异扫描驱动：逐函数注入 `return None`、跑分比对基线，支持断点续跑。

用法（先设 SWEEP_ROOT 指向工作树）:
  SWEEP_ROOT=tmp/deadcode python sweep_drive.py                 # 全量
  SWEEP_ROOT=tmp/deadcode python sweep_drive.py --start 0 --end 50
  SWEEP_ROOT=tmp/deadcode python sweep_drive.py --limit 10      # 冒烟
  SWEEP_ROOT=tmp/deadcode python sweep_drive.py --keys k.json   # 只跑指定 key
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import mutlib
import sweepcfg

sys.path.insert(0, str(sweepcfg.ROOT))
sys.path.insert(0, str(sweepcfg.SRC))
HERE = os.path.dirname(os.path.abspath(__file__))
TIMEOUT = 1800.0


def _env() -> dict:
    # .pyc 头只存 mtime(秒)+size；同一文件不同变异体插入字符数完全相同，
    # 同秒内两个变异体 size/mtime 全同 -> 第二个会复用第一个的 .pyc，静默测错。
    return dict(os.environ, PYTHONDONTWRITEBYTECODE="1", SWEEP_ROOT=str(sweepcfg.ROOT),
                PYTHONIOENCODING="utf-8")


def load_results() -> dict:
    if sweepcfg.RESULTS.is_file():
        return json.loads(sweepcfg.RESULTS.read_text(encoding="utf-8"))
    return {}


def run_mutant(rec: dict) -> dict:
    original = None
    try:
        original = mutlib.apply(rec)
    except SyntaxError as exc:
        return {"status": "mutate_syntax_error", "error": str(exc)}
    if original is None:
        return {"status": "mutate_unsupported"}
    tmp = sweepcfg.RUN_TMP
    try:
        tmp.unlink(missing_ok=True)
        proc = subprocess.run(
            [sys.executable, os.path.join(HERE, "sweep_run_one.py"), tmp.name],
            cwd=str(sweepcfg.ROOT), capture_output=True, text=True, timeout=TIMEOUT,
            env=_env(), encoding="utf-8", errors="replace",
        )
        if proc.returncode != 0 or not tmp.is_file():
            return {"status": "runner_crash", "error": (proc.stderr or "")[-800:]}
        return json.loads(tmp.read_text(encoding="utf-8"))
    except subprocess.TimeoutExpired:
        return {"status": "timeout"}
    finally:
        mutlib.restore(rec, original)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--keys", default=None, help="只跑这些 key 的 JSON 文件")
    args = ap.parse_args()

    recs = mutlib.collect()
    if args.keys or sweepcfg.ONESHOT.is_file():
        want = set(json.loads(open(args.keys or sweepcfg.ONESHOT, encoding="utf-8").read()))
        recs = [r for r in recs if r["key"] in want]
    else:
        recs = recs[args.start:args.end]
        if args.limit:
            recs = recs[:args.limit]

    res = load_results()
    todo = [r for r in recs if r["key"] not in res]
    print(f"候选 {len(recs)} 待跑 {len(todo)}", flush=True)

    before = mutlib.snapshot()
    t0 = time.perf_counter()
    for i, rec in enumerate(todo, 1):
        out = run_mutant(rec)
        out["qual"], out["file"] = rec["qual"], rec["file"]
        res[rec["key"]] = out
        if i % 10 == 0 or i == len(todo):
            sweepcfg.RESULTS.write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
            done = time.perf_counter() - t0
            rate = i / done if done else 0
            same = sum(1 for v in res.values() if v.get("status") == "same")
            diff = sum(1 for v in res.values() if v.get("status") == "diff")
            print(f"[{i}/{len(todo)}] {rec['file']}::{rec['qual']} -> {out['status']} "
                  f"(checked={out.get('checked')}) {rate:.2f} fn/s "
                  f"eta={(len(todo)-i)/rate/60 if rate else 0:.1f}min | same={same} diff={diff}",
                  flush=True)
    sweepcfg.RESULTS.write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")

    dirty = mutlib.verify_clean(before)
    print("源码还原校验:", "干净" if not dirty else f"⚠ 残留 {dirty}")


if __name__ == "__main__":
    main()
