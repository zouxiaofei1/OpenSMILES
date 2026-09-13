"""单个变异体的跑分：任何一行预测异于基线即停（早期退出）。

调用方（sweep_drive.py）负责先把变异写进源码，跑完由它还原。

status 含义:
  diff          预测改变 -> 该函数有用
  same          全量扫完预测不变 -> 候选无用
  import_error  变异让包加载不了 -> 导入期必需，有用（不是无用！）
  pool_error    子进程池崩溃
"""

from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import sweepcfg

sys.path.insert(0, str(sweepcfg.ROOT))


def main() -> int:
    out_path = sweepcfg.ROOT / (sys.argv[1] if len(sys.argv) > 1 else "_mutrun_tmp.json")
    t0 = time.perf_counter()
    result: dict = {"status": None, "checked": 0, "diff_row": None}

    # 先试导入。必须在 import benchmarks 之前：评测脚手架顶层就会 import 被测包，
    # 放后面的话 import 期崩溃会逃逸成 runner_crash 而不是 import_error。
    sys.path.insert(0, str(sweepcfg.SRC))
    try:
        __import__(sweepcfg.ENTRY)
    except Exception as exc:
        result.update(status="import_error", error=f"{type(exc).__name__}: {exc}",
                      secs=round(time.perf_counter() - t0, 3))
        out_path.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
        return 0

    import sweep_base

    rows = sweep_base.load_rows()
    base = json.loads(sweepcfg.BASE_ROWS.read_text(encoding="utf-8"))
    order = json.loads(sweepcfg.PROBE_ORDER.read_text(encoding="utf-8"))
    smiles = [str(r.get("smiles") or "") for r in rows]
    payload = ((i, smiles[i]) for i in order)

    checked = 0
    try:
        with Pool(processes=sweepcfg.WORKERS, initializer=sweep_base._init_worker) as pool:
            for idx, en, zh in pool.imap_unordered(sweep_base.name_row, payload, chunksize=1):
                checked += 1
                if [en, zh] != base[idx]:
                    result.update(status="diff", diff_row=idx, base_pred=base[idx], mut_pred=[en, zh])
                    try:
                        from importlib import import_module

                        score = import_module(sweepcfg.BENCH_MODULE).score_record
                        result["base_score"] = score(base[idx][0], base[idx][1], rows[idx])
                        result["mut_score"] = score(en, zh, rows[idx])
                    except Exception:
                        pass
                    break
            else:
                result["status"] = "same"
    except Exception as exc:
        result.update(status="pool_error", error=f"{type(exc).__name__}: {exc}")
    result["checked"] = checked
    result["secs"] = round(time.perf_counter() - t0, 3)
    out_path.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
