# IUPAC: P-61.1.3,P-65.1.1.2
# Layer: L1,L2,L3
"""残基/取代基内的游离 -COOH(非主官能团)须前缀化为 carboxy,不得并入残基主链
并被读成同碳「n-hydroxy-n-oxo-」逐原子描述。

- 自由基/酰基主基团(p41=1)压制非主羧酸(p41=7)时，该 COOH 是取代基叶（P-61.1.3
  carboxy 前缀）：fragment `*CH2CH2COOH` 应命名 `2-carboxyethyl`，`*CH2CH2CH2COOH`
  的酰基形态为 `3-carboxypropanoyl`，而非 `3-hydroxy-3-oxopropyl`/`4-hydroxy-4-oxobutanoyl`
  把酸碳并进链长。
- 整分子主官能团二酸/环酸（母体后缀 -dioic/-carboxylic acid）不降级，不得 carboxy 化。
- 降级酰胺/醛/酮羰基保持链化 oxo（见 test_demoted_amide_to_amino.py），不在此列。
"""

from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer

# 片段级正例：* 锚自由基残基，尾端游离 COOH 不并入主链。
POS_FRAG = [
    ("*CCC(=O)O", "2-carboxyethyl"),
    ("*C(=O)CCC(=O)O", "3-carboxypropanoyl"),
]


@pytest.mark.parametrize("smiles,en", POS_FRAG)
def test_carboxy_fragment(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 整分子正例：N-酰基残基尾端 COOH（peptide/代谢物场景），羧酸主母体不受影响。
POS_WHOLE = [
    ("CC(C)C[C@H](N=C(O)[C@@H](N)CC(=O)O)C(=O)O",
     "(2S)-2-[[(2S)-2-amino-3-carboxypropanoyl]amino]-4-methylpentanoic acid"),
    ("O=C(O)CCC(O)=N[C@@H](CCCC(=O)C(=O)O)C(=O)O",
     "(2S)-2-(3-carboxypropanoylamino)-6-oxoheptanedioic acid"),
]


@pytest.mark.parametrize("smiles,en", POS_WHOLE)
def test_carboxy_whole_molecule(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


# 负例：carboxy 规则不得误伤。
NEG = [
    # 降级酰胺羰基仍链化 oxo（伯酰胺 N → amino、羰基 C 留链作 4-oxo），不得 carboxy 化。
    ("NC(=O)CCC(=O)O", "4-amino-4-oxobutanoic acid"),
    # 主二酸母体不被降级、不得 carboxy 化，整分子仍为 dioic 母体后缀。
    ("O=C(O)CCC(=O)O", "butanedioic acid"),
]


@pytest.mark.parametrize("smiles,en", NEG)
def test_carboxy_negative_untouched(smiles: str, en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
