"""一次全量命名的函数覆盖率采集（sys.monitoring），输出每个函数被执行的行数。

只在**没有其它进程改写源码**时跑。若与变异扫描并发，读到的是掺杂了变异体的
混合状态，会得出「从未执行却改变预测」这类假矛盾。

监控必须在导入被测包**之前**装好，否则漏记"仅导入期调用"的函数（如模块级
`_TABLE = _build()` 里的 `_build`）——这类函数变异后会直接让模块导入失败。
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter

import sweepcfg

sys.path.insert(0, str(sweepcfg.ROOT))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

MARK = str(sweepcfg.SRC)
hits: Counter = Counter()
current: set = set()
_row = [-1]


def _on_start(code, offset):
    fn = code.co_filename
    if MARK in fn:
        rel = fn.split(sweepcfg.SRC.name, 1)[-1].lstrip("\\/").replace("\\", "/")
        key = f"{sweepcfg.SRC.name}/{rel}::{code.co_qualname}"
        hits[key] += 1
        current.add(key)


def main() -> None:
    mon = sys.monitoring
    tool = 2  # PROFILER 槽位
    mon.use_tool_id(tool, "sweep_cov")
    mon.register_callback(tool, mon.events.PY_START, _on_start)
    mon.set_events(tool, mon.events.PY_START)

    namer = sweepcfg.import_factory()()  # 监控装好之后再导入
    import sweep_base

    rows = sweep_base.load_rows()
    per_row: list[list[str]] = []
    t0, exc = time.perf_counter(), 0
    for r in rows:
        current.clear()
        try:
            namer.cache.clear()
            namer.name(str(r.get("smiles") or ""))
        except Exception:
            exc += 1
        per_row.append(sorted(current))
    secs = time.perf_counter() - t0

    mon.set_events(tool, 0)
    mon.free_tool_id(tool)

    sweepcfg.COVERAGE.write_text(json.dumps(dict(hits), ensure_ascii=False), encoding="utf-8")
    sweepcfg.COV_ROWS.write_text(json.dumps(per_row, ensure_ascii=False), encoding="utf-8")
    nf = sum(1 for k in hits if not k.endswith("::<module>"))
    nm = sum(1 for k in hits if k.endswith("::<module>"))
    print(f"行数 {len(rows)} 异常 {exc} 耗时 {secs:.1f}s | 执行到函数 {nf} 个、模块 {nm} 个")
    print("提示：覆盖行序变了，需重跑 sweep_base.py 以重排探针序")


if __name__ == "__main__":
    main()
