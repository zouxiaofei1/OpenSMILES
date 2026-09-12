# IUPAC: P-16.3.2（乘法前缀选择：被取代/复合组分用 bis/tris，简单组分用 di/tri）
# Layer: L5 (assembler_prefixes)
"""相同复合取代基重复时倍增前缀须用 bis/tris 而非 di/tri。

修复前 namer 仅对含 'carboxy' 的词干用 bis，漏掉 hydroxymethyl（retained 组合叶）
与 recursive 复合前缀（trimethoxyphenyl/bromobutyl/hydroxyphenyl/hydroxyethoxy 等），
输出 di(...) 而 gold/IUPAC 为 bis(...)。以下 7 例修复后应与 gold 全串一致。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer

FULL_CASES = [
    # A12 代表：retained 组合叶 hydroxymethyl 与 recursive 复合前缀
    ("C/C(=C\\CC/C(C)=C/CC/C=C(\\C)CC/C=C(\\C)CCC=C(CO)CO)CCC=C(CO)CO",
     "(6E,10E,14E,18E)-2,23-bis(hydroxymethyl)-6,10,15,19-tetramethyltetracosa-2,6,10,14,18,22-hexaene-1,24-diol"),
    ("OCC[C@@H]1C(CO)=C(CO)C[C@H]1O",
     "(1R,2R)-2-(2-hydroxyethyl)-3,4-bis(hydroxymethyl)cyclopent-3-en-1-ol"),
    ("COC1=C(C=CC(=C1OC)OC)\\C=C\\C(CC(\\C=C\\C1=C(C(=C(C=C1)OC)OC)OC)=O)=O",
     "(1E,6E)-1,7-bis(2,3,4-trimethoxyphenyl)hepta-1,6-diene-3,5-dione"),
    ("BrCCCCC1=CC=C(C=C1)CCCCBr",
     "1,4-bis(4-bromobutyl)benzene"),
    ("O=C(/C=C/c1ccc(O)cc1)CC(=O)/C=C/c1ccc(O)cc1",
     "(1E,6E)-1,7-bis(4-hydroxyphenyl)hepta-1,6-diene-3,5-dione"),
    ("CCCCCCCCCCCCCCCCCC(=O)OCCOCC(OCCO)[C@H]1OCC(OCCO)[C@H]1OCCO",
     "2-[2-[(2R,3R)-3,4-bis(2-hydroxyethoxy)oxolan-2-yl]-2-(2-hydroxyethoxy)ethoxy]ethyl octadecanoate"),
    # A12 家族第 7 条：方括号复杂基团 1-(2-methylbutan-2-yl)indol-3-yl 双取代
    ("CCC(C)(C)n1cc(C2=C(O)C(=O)C(c3cn(C(C)(C)CC)c4ccccc34)=C(O)C2=O)c2ccccc21",
     "2,5-dihydroxy-3,6-bis[1-(2-methylbutan-2-yl)indol-3-yl]cyclohexa-2,5-diene-1,4-dione"),
]


@pytest.mark.parametrize("smiles,en", FULL_CASES)
def test_bis_compound_mult(smiles: str, en: str) -> None:
    """复合取代基倍增应输出 bis(...) 且整串与 gold 一致。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
