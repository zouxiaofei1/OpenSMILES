# IUPAC: P-22.2.2
# Layer: L2,L3,L4,L5
"""生成式 Hantzsch-Widman 母体：端到端（未命中保留模板的 3-10 元杂单环）。"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# (smiles, 期望英文, 期望中文)
RING_CASES = [
    # 小环与饱和环
    ("C1CSC1", "thietane", "硫杂环丁烷"),
    ("C1CCCCCCN1", "azocane", "氮杂环辛烷"),
    ("C1CSCCN1", "1,4-thiazinane", "1,4-噻嗪烷"),
    ("C1CCCCOCC1", "oxocane", "氧杂环辛烷"),
    # 多杂原子：位次按引用顺序（P-22.2.2.1.3）
    ("C1COCSC1", "1,3-oxathiane", "1,3-氧杂硫杂环己烷"),
    ("C1OCOO1", "1,2,4-trioxolane", "1,2,4-三氧杂环戊烷"),
    ("S1SSSSSSS1", "octathiocane", "八硫杂环辛烷"),
    # 应保留的（模板命中，生成器不得抢入）
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("C1COCO1", "1,3-dioxolane", "1,3-二氧戊环"),
    ("c1ccoc1", "furan", "呋喃"),
]

FG_CASES = [
    ("O=C(O)C1CSCCN1", "1,4-thiazinane-3-carboxylic acid", "1,4-噻嗪烷-3-羧酸"),
    ("C1CC(=O)NCCS1", "1,4-thiazepan-5-one", "1,4-硫杂氮杂环庚烷-5-酮"),
    ("CN1CCSCC1", "4-methyl-1,4-thiazinane", "4-甲基-1,4-噻嗪烷"),
]


@pytest.mark.parametrize("smiles,en,zh", RING_CASES)
def test_hw_ring_parent(smiles: str, en: str, zh: str) -> None:
    """母体词干由生成器给出，杂原子不得丢失。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", FG_CASES)
def test_hw_ring_with_principal_group(smiles: str, en: str, zh: str) -> None:
    """环上/环外主基团与生成式词干拼接。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_partial_saturation_uses_hydro_prefix() -> None:
    """部分饱和环：mancude 词干 + hydro 前缀（P-54.4.1）。"""
    r = SMILESNNamer().name("C1CC=CCP1")
    assert r.success, r
    assert normalize_en(r.en) == normalize_en("1,2,3,6-tetrahydrophosphinine")


def test_saturated_ring_has_no_hydro_prefix() -> None:
    """全饱和环取饱和词干，不得被误标 dihydro/tetrahydro。"""
    for smiles in ("C1CSC1", "C1COCSC1", "C1CSCCN1"):
        en = SMILESNNamer().name(smiles).en or ""
        assert "hydro" not in en, (smiles, en)
