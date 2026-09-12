# IUPAC: P-66.6.1.1.3 (-CHO 连环/环系 → -carbaldehyde)、P-66.6.1.1(3) (-CHO 作取代基前缀 → formyl)
# Layer: L2,L4,L5 (+ L3 formyl 前缀 / anchored_table)
"""环上外环 -CHO：作主官能团时 -carbaldehyde（单/多），降级为取代基前缀时 formyl。

覆盖环骨架（吡啶/苯/环烷）、多基团 -dicarbaldehyde，及保留名
benzaldehyde 不回归；formyl 前缀取代 oxomethyl（P-66.6.1.1(3)）。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")
CARBALDEHYDE = [
    # 吡啶环 + 甲氨基取代，醛 locant 由 N 起算定向（tiers-47927 黄金名）
    ("CNC1=C(C=O)C=CC=N1", "2-(methylamino)pyridine-3-carbaldehyde", "2-(甲氨基)吡啶-3-甲醛"),
    ("O=CC1=CC=CC=N1", "pyridine-2-carbaldehyde", "吡啶-2-甲醛"),
    ("O=CC1=NC(C=O)=CC=C1", "pyridine-2,6-dicarbaldehyde", "吡啶-2,6-二甲醛"),
    # 苯环多醛 → 系统名 dicarbaldehyde（非保留 phthalaldehyde）
    ("C1C(C=O)=C(C=O)C=CC=1", "benzene-1,2-dicarbaldehyde", "苯-1,2-二甲醛"),
    ("O=Cc1ccc(C=O)cc1", "benzene-1,4-dicarbaldehyde", "苯-1,4-二甲醛"),
    # 单环环烷烃：醛基位次 1 隐含省略
    ("O=CC1CCCCC1", "cyclohexanecarbaldehyde", "环己烷甲醛"),
    # 哒嗪-3-甲醛（两个相邻等价环 N，方向不能随 SMILES 书写顺序变，
    # P-14.4(c)：后缀醛占低位次 3，CF3 前缀只得 6 —— 同分子两种写法必须同名）
    ("FC(C1=CC=C(N=N1)C=O)(F)F", "6-(trifluoromethyl)pyridazine-3-carbaldehyde", "6-(三氟甲基)哒嗪-3-甲醛"),
    ("C1(C(F)(F)F)N=NC(C=O)=CC=1", "6-(trifluoromethyl)pyridazine-3-carbaldehyde", "6-(三氟甲基)哒嗪-3-甲醛"),
]

FORMYL_PREFIX = [
    ("OC(=O)c1ccc(cc1)C=O", "4-formylbenzoic acid", "4-甲酰苯甲酸"),
    ("OC(=O)c1ccccc1C=O", "2-formylbenzoic acid", "2-甲酰苯甲酸"),
    ("NC(=O)c1ccc(C=O)cc1", "4-formylbenzamide", "4-甲酰苯甲酰胺"),
]

# 保留名 / 开链醛不回归（P-66.6.1.2 苯甲醛系保留名；不误伤链醛）。
REGRESSION = [
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("Cc1ccc(C=O)cc1", "4-methylbenzaldehyde", "4-甲基苯甲醛"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CARBALDEHYDE)
def test_ring_carbaldehyde(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", FORMYL_PREFIX)
def test_formyl_prefix(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", REGRESSION)
def test_regression(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
