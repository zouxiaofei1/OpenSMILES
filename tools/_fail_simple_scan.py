"""Find short/simple dual-fail cases for topic selection."""
from __future__ import annotations

import re
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "src"))

from benchmarks.benchmark import _is_fail, _load_rows  # noqa: E402
from benchmarks.benchmark_parallel import (  # noqa: E402
    _chunksize,
    _default_workers,
    _init_worker,
    _score_row,
)


def main() -> None:
    rows = _load_rows(Path("data/merged_benchmark.json"), None)
    workers = _default_workers()
    cs = _chunksize(len(rows), workers)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        results = list(ex.map(_score_row, rows, chunksize=cs))

    simple = []
    for score, pred_en, pred_zh, row in results:
        if not _is_fail(score):
            continue
        en = (row.get("english_name") or "").strip()
        if not en or len(en) > 50:
            continue
        if re.search(r"[\(\)\[\]]", en):
            continue
        if en.count("-") > 3:
            continue
        smi = row.get("smiles") or ""
        if len(smi) > 60:
            continue
        simple.append(
            {
                "en": en,
                "zh": row.get("chinese_name") or "",
                "pred": pred_en or "",
                "pred_zh": pred_zh or "",
                "smi": smi,
                "eval_en": row.get("eval_en"),
                "eval_zh": row.get("eval_zh"),
                "en_ok": score.get("en_ok"),
                "zh_ok": score.get("zh_ok"),
            }
        )

    simple.sort(key=lambda x: (len(x["en"]), len(x["smi"])))
    print("simple fails", len(simple))
    for s in simple[:80]:
        print(
            f"gold={s['en']!r} zh={s['zh'][:24]!r} "
            f"pred={s['pred']!r} smi={s['smi']}"
        )

    # token prefixes
    pref = Counter()
    for s in simple:
        tok = re.split(r"[\s\-]+", s["en"].lower())[0]
        pref[tok] += 1
    print("\nprefixes:")
    for k, v in pref.most_common(30):
        print(k, v)


if __name__ == "__main__":
    main()
