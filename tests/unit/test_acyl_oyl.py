# IUPAC: P-64.5.1,P-65.1.7.2,P-66.1.1.4.3
# Layer: L1,L2,L3,L4,L5
"""酰基/酰胺残基（自由价锚定在羰基碳的 R-C(=O)- 片段）必须以酸衍生酰基名命名，
而非展开成「1-氧代烷基」。

- P-64.5.1(2)：侧链 1 位羰基（自由价所在碳，即 -CO-R）必须用相应酰基名；官方反例
  [not 5-(1-oxoethyl)nonane-4,6-dione]。
- P-65.1.7.2：酰基名由对应酸名把 -oic acid→-oyl / -ic acid→-yl 改词尾而派生
  （acetyl、butanoyl、2-thiophen-2-ylacetyl…）。
- P-66.1.1.4.3：R-CO-NH- 残基作取代基可用 acylamino（方法 2，本簇 gold 采用）。

说明：ChEBI gold_zh 对 N-酰基残基常写作「…酰胺基」(amido 式方法 1) 且中英不对齐，
与本规则实现的方法(2) acylamino 不同写法并存；故 N-酰基行的中文断言跳过（zh=None），
只精确断言英文。O-酰基行同理以英文为准。
"""

from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# 片段级正例：锚定酰基残基的自由名（碳锚点，α 为非环碳）。
POS_FRAG = [
    ("*C(=O)C", "acetyl", "乙酰基"),
    ("*C(=O)CCC", "butanoyl", None),
    ("*C(=O)Cc1cccs1", "2-thiophen-2-ylacetyl", "2-噻吩-2-基乙酰基"),
    ("*C(=O)N", "carbamoyl", "氨基甲酰基"),  # 氨基甲酸的酰基保留前缀（P-66.1.1.4.1），锚定表条目，不走 oxo/amino 展开
]


@pytest.mark.parametrize("smiles,en,zh", POS_FRAG)
def test_acyl_fragment_named_oyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子正例（gold/pred 仅差 1-oxo→oyl 单点；英文精确，中文跳过/待验）。
# 注：O-酰基的围栏按 P-65.6.3.2.3：retained 单词酰基名（acetyl/benzoyl/hexadecanoyl）
# 与 oxy 融合成 …oyloxy（3-(benzoyloxy)propanoic acid）；带位次的系统酰基名则闭括号在
# 前端之后、oxy 留括号外（3-[(pyridine-3-carbonyl)oxy]propanoic acid）。故带双键位次/
# 立体描述符的 (9Z)-octadec-9-enoyl 取后者，(9Z)-octadec-9-enoyl]oxy；与 benchmark gold
# （chebi-2873/chebi-776）一致。
POS_WHOLE = [
    ("CC(C)=CC(O)=NCC(=O)O", "2-(3-methylbut-2-enoylamino)acetic acid"),
    ("CCCCCCCCCCCCCCCCCC(O)=NCC(=O)O", "2-(octadecanoylamino)acetic acid"),
    ("C(C(=C)C)(=O)OC(C)OC(C(=C)C)=O",
     "1-(2-methylprop-2-enoyloxy)ethyl 2-methylprop-2-enoate"),
    ("CCCCCCCC/C=C\\CCCCCCCC(=O)OC(CCCCCCCCCCCCC)CCCC(=O)O",
     "5-[(9Z)-octadec-9-enoyl]oxyoctadecanoic acid"),
    ("CCCCCC/C=C\\CCCCCCCC(=O)OC(CCCCCCCCC)CCCCCCCC(=O)O",
     "9-[(9Z)-hexadec-9-enoyl]oxyoctadecanoic acid"),
]


@pytest.mark.parametrize("smiles,en", POS_WHOLE)
def test_acyl_whole_molecule_oyl(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 负例：acyl 规则不得误伤。
# 1) 酰胺/酮作母体的既有路径不变；2) 内部氧代（锚点不在羰基碳）不变；
# 3) 锚点羰基碳带 O 离去基（酯羰基）不判酰基，保持现行为（带 N 的氨基甲酰见上方 POS_FRAG）。
NEG = [
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("CC(=O)c1ccccc1", "1-phenylethanone", "1-苯基乙酮"),
    ("*CCC(=O)C", "3-oxobutyl", "3-氧代丁基"),
    ("*C(=O)OC", "methoxy(oxo)methyl", "甲氧基(氧代)甲基"),
]


@pytest.mark.parametrize("smiles,en,zh", NEG)
def test_acyl_negative_untouched(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
