# IUPAC: P-14.3.4.4 / P-66 (acetonitrile 省略例) / P-69 L136 (2-substituted ethanol)
# Layer: L5 (assembler_prefixes._omit_sub_locants)
"""C2 单取代母体 locant 省略收窄：醇/硫醇/胺 FG 端碳有可取代 H，2- 不可省略。

原逻辑对所有 C2 单取代母体一刀切省略（cyclopropylethanol、aminoethanol、chloroethanol）。
需排除 kind ∈ {alcohol, amine, thiol}；而端碳无 H 的腈/酸等仍可省略（(1H-indol-5-yl)acetonitrile）。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# 必须补上 2- 的 C2 醇/胺
FULL_CASES = [
    ("NCCO", "2-aminoethanol", "2-氨基乙醇"),
    ("ClCCO", "2-chloroethanol", "2-氯乙醇"),
    ("C1(CC1)CCO", "2-cyclopropylethanol", "2-环丙基乙醇"),
]

# 仍允许省略（端碳无 H / 对称）—— 防止过度收窄
KEEP_OMIT_CASES = [
    ("N1C=CC2=CC(=CC=C12)CC#N", "(1H-indol-5-yl)acetonitrile"),
]


@pytest.mark.parametrize("smiles,en,zh", FULL_CASES)
def test_c2_fg_keeps_locant(smiles: str, en: str, zh: str) -> None:
    """C2 醇/胺单取代母体：2- 必须保留。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en", KEEP_OMIT_CASES)
def test_c2_fg_still_omits(smiles: str, en: str) -> None:
    """端碳无可取代 H 的腈母体：2- 省略仍成立。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
