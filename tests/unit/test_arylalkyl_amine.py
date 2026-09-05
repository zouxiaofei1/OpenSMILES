# IUPAC: P-62.2.2.1 / P-29.3
# Layer: L2,L3
"""Sec-amine: N-端取代基命名（简单 2° 胺）与芳烷基对照。

含芳基臂的二级胺（N-ethyl-1-phenylethanamine 等）通过 L2 P-45.2.1 流水线
选出前缀取代基团数目更多（芳环臂作为母体骨架取代基）的母体链，主链臂
优先保留为母体、另一臂作 N-取代基。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # nitriles from user set (should already work)
    (
        "C(C)C=1C=C(C=CC1)CC#N",
        "2-(3-ethylphenyl)acetonitrile",
        "2-(3-乙基苯基)乙腈",
    ),
    (
        "ClC1=C(C(=CC=C1)F)CC#N",
        "2-(2-chloro-6-fluorophenyl)acetonitrile",
        "2-(2-氯-6-氟苯基)乙腈",
    ),
    # ester: prefix inserts between alkyl and acyl (user ex7)
    (
        "COC(C(CCC1=CC=CC=C1)=O)=O",
        "methyl 2-oxo-4-phenylbutanoate",
        "2-氧代-4-苯基丁酸甲酯",
    ),
    # 二级胺：N- 取代基前缀
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
    ("NCCc1ccccc1", "2-phenylethanamine", "2-苯基乙胺"),
    # 含芳基臂的二级胺：母体选前缀取代基更多的链（P-45.2.1），芳环臂作母体骨架取代基
    (
        "C(C)NC(C)C1=CC=C(C=C1)OC",
        "N-ethyl-1-(4-methoxyphenyl)ethanamine",
        "N-乙基-1-(4-甲氧基苯基)乙胺",
    ),
    (
        "C(C)NC(C)c1ccccc1",
        "N-ethyl-1-phenylethanamine",
        "N-乙基-1-苯基乙胺",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arylalkyl_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
