# IUPAC: P-93.4 / P-91.2
# Layer: L5 (chain_engine._KIND_TABLE ez 字段)
"""链式酮/胺/硫醇母体的烯 E/Z 立体描述符。

_KIND_TABLE 里 alcohol 配了 ez_ene/ez_ene_multi(=ez_for_parent) 而 amine/thiol 全缺、
ketone 缺 ez_ene_multi —— 结果带立体双键的烯胺/烯硫醇(连单烯)与 多烯+二官能团(二酮/
二胺/二硫醇)都丢 E/Z 前缀。本文件锁四 kind 对齐醇后的前缀行为。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh") — 正例:修复后应带 E/Z 前缀
CASES = [
    # 多烯 + 多官能团(unsat_polyol 多烯分支):四 kind
    (
        "C=CC(O)C/C=C\\C(O)/C=C\\CC",
        "(5Z,8Z)-undeca-1,5,8-triene-3,7-diol",
        "(5Z,8Z)-十一烷-1,5,8-三烯-3,7-二醇",
    ),
    (
        "C=CC(N)C/C=C\\C(N)/C=C\\CC",
        "(5Z,8Z)-undeca-1,5,8-triene-3,7-diamine",
        "(5Z,8Z)-十一烷-1,5,8-三烯-3,7-二胺",
    ),
    (
        "C=CC(S)C/C=C\\C(S)/C=C\\CC",
        "(5Z,8Z)-undeca-1,5,8-triene-3,7-dithiol",
        "(5Z,8Z)-十一烷-1,5,8-三烯-3,7-二硫醇",
    ),
    (
        "C=CC(=O)C/C=C\\C(=O)/C=C\\CC",
        "(5Z,8Z)-undeca-1,5,8-triene-3,7-dione",
        "(5Z,8Z)-十一-1,5,8-三烯-3,7-二酮",
    ),
    # 单烯 胺/硫醇:曾连单烯 E/Z 都缺(表无 ez_ene)
    (
        "C/C=C/CN",
        "(2E)-but-2-en-1-amine",
        "(2E)-丁-2-烯-1-胺",
    ),
    (
        "C/C=C/CS",
        "(2E)-but-2-ene-1-thiol",
        "(2E)-丁-2-烯-1-硫醇",
    ),
    # 融合式多烯 kind(ene_base):酰胺/腈/酰氯缺 ez_ene_multi→多烯无 E/Z(单烯早已有)
    (
        "NC(=O)/C=C/C=C/C",
        "(2E,4E)-hexa-2,4-dienamide",
        "(2E,4E)-己-2,4-二烯酰胺",
    ),
    (
        "N#CC/C=C\\C/C=C\\C",
        "(3Z,6Z)-octa-3,6-dienenitrile",
        "(3Z,6Z)-辛-3,6-二烯腈",
    ),
    (
        "O=C(Cl)/C=C/C=C/C",
        "(2E,4E)-hexa-2,4-dienoyl chloride",
        "(2E,4E)-己-2,4-二烯酰氯",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_kind_ene_ez_prefix(smiles: str, en: str, zh: str) -> None:
    """定义立体化学的双键必须带 (E)/(Z) 前缀。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# negatives — 无立体标记不得伪造 E/Z 前缀
@pytest.mark.parametrize(
    "smiles,en,zh",
    [
        ("CC=CCN", "but-2-en-1-amine", "丁-2-烯-1-胺"),
        ("CC=CC(O)CO", "pent-3-ene-1,2-diol", "戊烷-3-烯-1,2-二醇"),
        ("CC=CC(=O)CC(=O)CC", "oct-6-ene-3,5-dione", "辛-6-烯-3,5-二酮"),
    ],
)
def test_plain_alkene_no_ez(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
