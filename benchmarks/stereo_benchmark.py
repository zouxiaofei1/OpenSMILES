"""立体命名 benchmark：只比名称里的立体描述符 token（R/S/E/Z，含 fused 位次如 3aR），
忽略整名其它差异（前缀排序、氧桥括号写法），在 benchmark 长立体名上量"立体正确比例"。
可反复运行：存快照（默认 benchmarks/.last_stereo.json），下次自动 diff 报每行 REGRESS/IMPROVE。
用法：python -m benchmarks.stereo_benchmark --data benchmarks/merged_benchmark.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from collections import Counter
from contextlib import redirect_stdout
from pathlib import Path

# 未 pip install -e 时也能从仓库根直接跑：把 src/ 加进 sys.path（与 tests/conftest.py 同法）
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# 静音 src 命名管线里的调试 print（Orientation/环元组等），保持输出干净
_DEVNULL = open(os.devnull, "w", encoding="utf-8")

from namepredict.namer import SMILESNNamer

# 立体描述符 token：可带位次（可含 fused 小写位次如 3a/6a）或裸 E/Z/R/S
_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9])[A-Za-z0-9]*[RSEZ](?![A-Za-z0-9])")

_SNAPSHOT_PATH = Path(__file__).resolve().parent / ".last_stereo.json"


def stereo_tokens(name: str) -> list[str]:
    """抽取名称中全部立体 token（原序）；无立体返回空列表。"""
    return _TOKEN_RE.findall(name or "")


def _load_rows(path: str | Path) -> list[dict]:
    """读取 benchmark JSON（list[dict]），容错单 dict。"""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return list(data.values()) if isinstance(data, dict) else data


def _row_key(r: dict, i: int) -> str:
    """行主键：优先 id，缺失用 smiles，再缺失用序号。"""
    return str(r.get("id") or r.get("smiles") or f"row-{i}")


_NAMER = None  # 惰性单例命名器（进程内复用；跨根片段缓存已在 namer 内按整分子隔离）


def _name_row(r: dict) -> dict:
    """对单行命名并返回立体比对字段（异常记 errored）。"""
    ref_t = stereo_tokens(r.get("english_name") or "")
    out = {"ref_tokens": sorted(ref_t), "our_tokens": [], "pred_en": "",
           "errored": False, "stereo_equal": False, "count_equal": False, "flip_free": False}
    if not ref_t:
        return out
    try:
        rr = _namer().name(r["smiles"])
    except Exception:
        out["errored"] = True
        return out
    our_t = stereo_tokens(rr.en if rr.success else "")
    out["our_tokens"] = sorted(our_t)
    out["pred_en"] = rr.en if rr.success else ""
    if not (rr.success and rr.en and our_t):
        return out
    equal = Counter(our_t) == Counter(ref_t)
    out["stereo_equal"] = equal
    if len(our_t) == len(ref_t):
        out["count_equal"] = True
        out["flip_free"] = equal
    return out


def _namer() -> SMILESNNamer:
    """惰性单例命名器（进程内复用；跨根片段缓存已在 namer 内按整分子隔离）。"""
    global _NAMER
    if _NAMER is None:
        _NAMER = SMILESNNamer()
    return _NAMER


def _quiet_worker() -> None:
    """工作进程初始化：预建命名器并静音 src 里的调试 print，避免刷屏。"""
    global _NAMER
    _NAMER = SMILESNNamer()
    sys.stdout = open(os.devnull, "w", encoding="utf-8")


def run_rows(rows: list[dict], *, workers: int = 0) -> dict:
    """逐行命名并汇总 report：含 items{key: row_result} 供快照 diff。"""
    if workers and workers > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=workers, initializer=_quiet_worker) as ex:
            results = list(ex.map(_name_row, rows))
    else:
        with redirect_stdout(_DEVNULL):  # 静音命名管线调试 print
            results = [_name_row(r) for r in rows]
    stats = {"rows": len(rows), "named_ok": 0, "errored": 0, "our_stereo": 0,
             "stereo_equal": 0, "count_equal": 0, "flip_free": 0}
    items = {}
    for i, (r, res) in enumerate(zip(rows, results)):
        items[_row_key(r, i)] = res
        if res["errored"]:
            stats["errored"] += 1
            continue
        if res["our_tokens"]:
            stats["our_stereo"] += 1
        if res["our_tokens"]:
            stats["named_ok"] += 1
        if res["stereo_equal"]:
            stats["stereo_equal"] += 1
        if res["count_equal"]:
            stats["count_equal"] += 1
            if res["flip_free"]:
                stats["flip_free"] += 1
    stats["pure_ratio"] = stats["flip_free"] / (stats["count_equal"] or 1)
    return {"stats": stats, "items": items}


def _build_snapshot(report: dict) -> dict:
    """从 report 抽可 diff 快照（顶层计数 + 每行立体状态）。"""
    return {
        "stats": report["stats"],
        "items": {
            k: {f: v.get(f) for f in ("our_tokens", "stereo_equal", "count_equal")}
            for k, v in report["items"].items()
        },
    }


def _load_snapshot(path: Path) -> dict | None:
    """读上次快照，损坏返回 None。"""
    if not path.is_file():
        return None
    try:
        with path.open(encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def _save_snapshot(report: dict, path: Path) -> None:
    """写本次快照。"""
    with path.open("w", encoding="utf-8") as f:
        json.dump(_build_snapshot(report), f, ensure_ascii=False)


def _diff(prev: dict, cur: dict) -> list[dict]:
    """与上次快照对比，给每行 REGRESS/IMPROVE/CHANGE 分类。"""
    diffs = []
    p_items, c_items = prev.get("items") or {}, cur.get("items") or {}
    for key, c in c_items.items():
        p = p_items.get(key)
        if p is None:
            continue
        p_ok, c_ok = bool(p.get("stereo_equal")), bool(c.get("stereo_equal"))
        if p_ok and not c_ok:
            kind = "REGRESS"
        elif not p_ok and c_ok:
            kind = "IMPROVE"
        elif (p.get("our_tokens") or []) == (c.get("our_tokens") or []):
            continue
        else:
            kind = "CHANGE"
        diffs.append({"key": key, "kind": kind,
                      "prev_ok": p_ok, "cur_ok": c_ok,
                      "prev_tokens": p.get("our_tokens") or [],
                      "cur_tokens": c.get("our_tokens") or []})
    return diffs


def _print_summary(report: dict, elapsed: float) -> None:
    """打印立体 benchmark 汇总（口径与主 benchmark 对齐，见 run_rows.stats）。"""
    s = report["stats"]
    print(f"rows={s['rows']} named_ok={s['named_ok']} errored={s['errored']} "
          f"our_with_stereo={s['our_stereo']}")
    print(f"stereo_equal={s['stereo_equal']} "
          f"(/{s['our_stereo']})={s['stereo_equal']/max(s['our_stereo'],1):.1%}")
    print(f"count_equal={s['count_equal']} flip_free={s['flip_free']} "
          f"pure_stereo_ok_ratio={s['pure_ratio']:.1%}")
    print(f"elapsed={elapsed:.1f}s")


def _print_diffs(diffs: list[dict]) -> None:
    """打印与上次快照的差异（改进/回退逐行、只显示 token 数避免刷屏）。"""
    if not diffs:
        print("diff_vs_last: none (same as last run, or no previous snapshot)")
        return
    regress = [d for d in diffs if d["kind"] == "REGRESS"]
    improve = [d for d in diffs if d["kind"] == "IMPROVE"]
    change = [d for d in diffs if d["kind"] == "CHANGE"]
    print(f"diff_vs_last: improve={len(improve)} regress={len(regress)} "
          f"change={len(change)} (n={len(diffs)})")
    for d in improve[:10]:
        print(f"   IMPROVE {d['key']}  tokens {len(d['prev_tokens'])} -> {len(d['cur_tokens'])}")
    for d in regress[:15]:
        print(f"   REGRESS {d['key']}  tokens {len(d['prev_tokens'])} -> {len(d['cur_tokens'])}")
    for d in change[:10]:
        print(f"   CHANGE {d['key']}  tokens {len(d['prev_tokens'])} -> {len(d['cur_tokens'])}")


def dump_errors(rows: list[dict], report: dict, path: Path) -> None:
    """把立体不匹配行导出为精简 JSON：含 key / smiles / pred_en / english_name。

    只导出"纯立体错"行：立体不匹配（stereo_equal=False）且 pred_en 与
    english_name 末尾三字母一致（尾缀功能团相同）、长度相差不超过一倍
    （较长 ≤ 较短两倍），以排除整名/主结构就错的非立体失败；命名报错/
    无产出的行不导出。smiles 直接随行给出。
    """
    items = report.get("items") or {}
    bad: list[dict] = []
    for i, row in enumerate(rows):
        res = items.get(_row_key(row, i))
        if res is None or res.get("errored"):
            continue
        if res.get("stereo_equal"):
            continue
        pred_en = res.get("pred_en") or ""
        if not pred_en:
            continue
        english_name = row.get("english_name") or ""
        if not english_name or pred_en[-3:] != english_name[-3:]:
            continue
        # 长度悬殊（较长 ≥ 较短两倍）说明整名级偏差，非纯立体错，排除
        if max(len(pred_en), len(english_name)) >= 2 * min(len(pred_en), len(english_name)):
            continue
        bad.append({
            "key": _row_key(row, i),
            "smiles": row.get("smiles") or "",
            "pred_en": pred_en,
            "english_name": english_name,
        })
    with open(path, "w", encoding="utf-8") as f:
        json.dump(bad, f, ensure_ascii=False, indent=1)
    print(f"errors written to {path}: {len(bad)}/{len(rows)}")


def _select_rows(data: list[dict], *, min_len: int, sample: int) -> list[dict]:
    """选立体行：参考名含立体 token；min_len 过滤名长；sample>0 取最长 sample 条（确定性）。"""
    rows = [
        r for r in data
        if (r.get("english_name") or "") and stereo_tokens(r.get("english_name") or "")
        and len(r.get("english_name") or "") >= min_len
    ]
    rows.sort(key=lambda r: -len(r.get("english_name") or ""))
    return rows[:sample] if sample else rows


def main(argv: list[str] | None = None) -> None:
    """CLI：全量/样本立体行测量，可写/比快照。"""
    p = argparse.ArgumentParser(description="立体命名 benchmark（比立体 token 正确比例）")
    p.add_argument("--data", type=str, default="benchmarks/merged_benchmark.json")
    p.add_argument("--min-len", type=int, default=0, help="参考英文名最小长度")
    p.add_argument("--sample", type=int, default=0, help="取名长最长 N 条（0=全部立体行）")
    p.add_argument("--workers", type=int, default=0, help="并行 worker 数（0=串行）")
    p.add_argument("--snapshot", type=Path, default=_SNAPSHOT_PATH)
    p.add_argument("--no-snapshot", action="store_true", help="不读写快照")
    p.add_argument("--errors", type=Path, default=None,
                   help="把立体不匹配/命名报错的行导出为 JSON")
    args = p.parse_args(argv)

    t0 = time.perf_counter()
    rows = _select_rows(_load_rows(args.data), min_len=args.min_len, sample=args.sample)
    report = run_rows(rows, workers=args.workers)
    elapsed = time.perf_counter() - t0
    if args.errors is not None:
        dump_errors(rows, report, args.errors)

    if args.no_snapshot:
        _print_summary(report, elapsed)
        return
    prev = _load_snapshot(args.snapshot)
    diffs = _diff(prev, report) if prev else []
    _print_summary(report, elapsed)
    if prev is None:
        print("diff_vs_last: none (first run)")
    else:
        _print_diffs(diffs)
    _save_snapshot(report, args.snapshot)


if __name__ == "__main__":
    main()
