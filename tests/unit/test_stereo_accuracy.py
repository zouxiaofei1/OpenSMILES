# IUPAC: P-92/P-93（立体），benchmark 长名样本
"""立体命名正确比例门禁：复用 benchmarks.stereo_benchmark 抽取 benchmark 超长含立体名称
样本，只比对立体描述符 token（R/S/E/Z，含 fused 位次如 3aR），忽略整名其它差异（前缀
排序、氧桥括号写法等）。用于"即使整名不全对也能量出立体正确比例"，并拦截整体 R/S
翻转/丢立体的大回归。可加 -s 查看样本统计。
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.stereo_benchmark import run_rows, stereo_tokens

ROOT = Path(__file__).resolve().parents[2]
_DATA = ROOT / "data" / "merged_benchmark.json"
_N_SAMPLE = 20  # 样本量：太长 pytest 慢；约 8s


def _sample_rows() -> list[dict]:
    """取英文名含立体 token 且长度最长的一批（长名多为深层糖苷/稠环，立体信息量大）。"""
    data = json.loads(_DATA.read_text(encoding="utf-8"))
    rows = [
        r for r in data
        if r.get("english_name") and stereo_tokens(r.get("english_name"))
    ]
    rows.sort(key=lambda r: -len(r.get("english_name") or ""))
    return rows[:_N_SAMPLE]


def test_stereo_token_extractor() -> None:
    """token 抽取只认立体描述符，不误吞普通单词。"""
    name = "(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl"
    assert stereo_tokens(name) == ["2S", "3R", "4S", "5S", "6R"]
    assert stereo_tokens("5-[(3aR,4R,5R,6aS)-5-hydroxy-2H-cyclopenta[b]furan]") == [
        "3aR", "4R", "5R", "6aS",
    ]
    assert stereo_tokens("but-3-en-1-ol") == []


@pytest.mark.skipif(not _DATA.exists(), reason="benchmark data missing")
def test_stereo_ok_ratio_floor() -> None:
    """长立体名样本上，描述符数一致的分子中"立体集合一致"占比不低于下限。"""
    report = run_rows(_sample_rows())  # 复用进程内单例命名器（run-scoped 片段缓存）
    s = report["stats"]
    ratio = s["flip_free"] / (s["count_equal"] or 1)
    print(
        f"stereo sample: rows={s['rows']} ok={s['named_ok']} errored={s['errored']} "
        f"count_equal={s['count_equal']} pure_stereo_ok={s['flip_free']} "
        f"ratio={ratio:.1%}"
    )
    # 下限取当前能力的一个余量，只挡"立体被整体丢弃/翻转"类回归；
    # 立体能力提升时应上调该下限，用 -s 看上述统计。
    assert ratio >= 0.20, (
        f"pure stereo ok ratio {ratio:.1%} below floor; stats={s}"
    )
