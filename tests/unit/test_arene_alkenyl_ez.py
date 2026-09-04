# IUPAC: P-93.4 / P-91.2 / P-29.2
# Layer: L3 (as_substituent / submol_build), L5 (chain_engine radical spec)
"""苯环母体上烯基取代基（含立体双键）的 E/Z 描述符。

根因是双重的：(1) submol_build 切子分子只重建键型，丢双键奇偶 → 递归命名
`*` 锚定的烯基片段时 RDKit 读不到立体；(2) _KIND_TABLE['radical'] 无 E/Z 拼接。
本文件锁：带斜线几何的烯基侧链应出 (1E)/(1Z)，无立体标记则保持无前缀。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # 正例：C=C 在取代基链、苯环作母体，双键立体须带前缀
    (r"C(/C)=C/c1ccccc1", "[(1Z)-prop-1-en-1-yl]benzene", "[(1Z)-丙-1-烯-1-基]苯"),
    (r"C/C=C\c1ccccc1", "[(1Z)-prop-1-en-1-yl]benzene", None),
    (r"C/C=C/c1ccccc1", "[(1E)-prop-1-en-1-yl]benzene", "[(1E)-丙-1-烯-1-基]苯"),
    (r"CC/C=C\c1ccccc1", "[(1Z)-but-1-en-1-yl]benzene", None),
    (r"CC/C=C/c1ccccc1", "[(1E)-but-1-en-1-yl]benzene", None),
    # 负例：无立体标记不得伪造 E/Z 前缀
    ("C=CCc1ccccc1", "prop-2-en-1-ylbenzene", None),
    ("CCC(C)c1ccccc1", "butan-2-ylbenzene", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_alkenyl_ez(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
