# IUPAC: P-65.1.7.2,P-64.5.1,P-66.6.1
# Layer: L1,L2,L5
"""环/芳羧酸的酰基残基 X-C(=O)-（X = 苯环/芳杂环/饱和碳环）必须以收缩酰基名命名，
而非展开成「yl(oxo)methyl / -oxomethyl」。

- P-65.1.7.2：酰基名由对应酸名派生——benzoic acid→benzoyl、furan-2-carboxylic acid→
  furan-2-carbonyl、cyclopropanecarboxylic acid→cyclopropanecarbonyl（把 -oic/-ic acid 词尾
  换成 -oyl/-carbonyl）；α 为环/芳碳时同样是酸衍生酰基，不能退化成逐原子 oxo 描述。
- P-64.5.1(2)：侧链 1 位羰基（自由价在羰基碳的 -CO-Ar）用相应酰基名。
- 既有 test_acyl_oyl.py 已覆盖 α 为开链烷基/烯基的酰基收缩（butanoyl/enoyl/octadecanoyl）；
  本文件补齐 α 为环/芳的缺口（A1 模式，benchmark ~50 例）。

与 test_acyl_oyl.py 同约定：N-酰基残基中文（gold 常「…酰胺基/甲酰…」与 IUPAC en 不对齐）跳过，
只精确断言英文；fragment 中文按 namer acyl 词形惯例。
"""

from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer

# 片段级正例：* 锚在酰基头羰基碳的环/芳酰基，应出环酸衍生酰基名。
POS_FRAG = [
    ("*C(=O)c1ccccc1", "benzoyl"),
    ("*C(=O)c1ccco1", "furan-2-carbonyl"),
    ("*C(=O)C1CC1", "cyclopropanecarbonyl"),
    ("*C(=O)c1ccc(Cl)cc1", "4-chlorobenzoyl"),
    ("*C(=O)C1CCCN1", "pyrrolidine-2-carbonyl"),
]


@pytest.mark.parametrize("smiles,en", POS_FRAG)
def test_ring_acyl_fragment_named_carbonyl(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例：N-酰基氨基残基（羧酸主母体 + -NH-CO-X），X 为环/芳时应收缩。
POS_WHOLE = [
    ("O=C(O)CNC(=O)c1ccco1", "2-(furan-2-carbonylamino)acetic acid"),
    ("C1CC1C(=O)NCC(=O)O", "2-(cyclopropanecarbonylamino)acetic acid"),
    ("O=C(O)CNC(=O)c1ccc(Cl)cc1", "2-(4-chlorobenzoylamino)acetic acid"),
]


@pytest.mark.parametrize("smiles,en", POS_WHOLE)
def test_ring_acyl_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
