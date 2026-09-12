# IUPAC: P-92 / P-93
# Layer: L5
"""饱和杂环保留母体(sp3)上的 CIP R/S：pyrrolidine/piperidine/morpholine/piperazine/oxolane/oxane 手性中心。
环母体 kind 不在链式 FG 白名单(rs_fgs)，曾整体丢掉 R/S 前缀。"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")  — RDKit CIP 校准；R/S 按系统自身位次编号
CASES = [
    # piperidine 族：2 位(杂原子邻位)与 3 位手性都覆盖
    (
        "C[C@H]1CCCCN1",
        "(2S)-2-methylpiperidine",
        "(2S)-2-甲基哌啶",
    ),
    (
        "C[C@@H]1CCCCN1",
        "(2R)-2-methylpiperidine",
        "(2R)-2-甲基哌啶",
    ),
    (
        "C[C@H]1CNCCC1",
        "(3R)-3-methylpiperidine",
        "(3R)-3-甲基哌啶",
    ),
    # 其余饱和杂环 scaffold，各取一方向
    (
        "C[C@@H]1CCCN1",
        "(2R)-2-methylpyrrolidine",
        "(2R)-2-甲基吡咯烷",
    ),
    (
        "C[C@H]1CCCO1",
        "(2S)-2-methyloxolane",
        "(2S)-2-甲基四氢呋喃",
    ),
    (
        "C[C@@H]1CCCCO1",
        "(2R)-2-methyloxane",
        "(2R)-2-甲基氧杂环己烷",
    ),
    (
        "C[C@@H]1COCCN1",
        "(3R)-3-methylmorpholine",
        "(3R)-3-甲基吗啉",
    ),
    (
        "C[C@@H]1CNCCN1",
        "(2R)-2-methylpiperazine",
        "(2R)-2-甲基哌嗪",
    ),
    # negatives — 无手性中心不伪造 R/S
    ("C[C@H]1CCCCC1", "methylcyclohexane", "甲基环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ring_rs_stereo(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
