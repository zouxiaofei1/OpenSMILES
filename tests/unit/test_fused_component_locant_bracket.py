# IUPAC: P-25.3.1.3
# Layer: L2,L5
"""附加组分名中的结构位次（杂原子位置）在稠合名里须置于方括号内。

P-25.3.1.3：「描述组分结构特征的位次，如杂原子位置，保留在组分名称中，并置于方
括号内」；P-25.3.2.1.2 进一步规定 isoxazole/oxazole/thiazole 在稠合名中必须改用
Hantzsch-Widman 名 1,2-oxazole / 1,3-oxazole / 1,3-thiazole，且位次置于方括号内。
故 1,2,4-triazole 作附加组分时是 [1,2,4]triazolo[1,5-a]pyridine，而非
1,2,4-triazolo[1,5-a]pyridine；该方括号位次集与母体间仍需连字符
（2,3-dihydro-[1,3]thiazolo…，同 1,3-dihydro-2H-… 的前导数字体例）。

方括号只在「附加组分/稠合母体」位置出现；同一杂环作取代基或作母体（1,3-thiazol-2-yl、
1,2-oxazole-3-carboxamide）时不加方括号，末条负例守此边界。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: 三唑并[1,5-a]吡啶（用户报例）
    (
        "BrC1=NN2C(C=CC=C2Br)=N1",
        "2,5-dibromo-[1,2,4]triazolo[1,5-a]pyridine",
        "2,5-二溴-[1,2,4]三唑并[1,5-a]吡啶",
    ),
    # positive: 三唑并[1,5-a]嘧啶，前缀（吡啶-2-基）与方括号母体间须留连字符
    (
        "Cc1nc(C)c(C(O)=Nc2ccc(F)c(-c3nc4ncc(-c5ccccn5)cn4n3)c2)o1",
        "N-[4-fluoro-3-(6-pyridin-2-yl-[1,2,4]triazolo[1,5-a]pyrimidin-2-yl)phenyl]-2,4-dimethyl-1,3-oxazole-5-carboxamide",
        "N-[4-氟-3-(6-吡啶-2-基-[1,2,4]三唑并[1,5-a]嘧啶-2-基)苯基]-2,4-二甲基-1,3-噁唑-5-甲酰胺",
    ),
    # positive: 噻唑并[3,2-a]嘧啶，二氢前缀与方括号之间连字符（中文侧另存括号差异，此处只断英文）
    (
        "O=C1C=CN=C2N1C(CS2)CC(=O)NC=2C=NC=CC2",
        "2-(5-oxo-2,3-dihydro-[1,3]thiazolo[3,2-a]pyrimidin-3-yl)-N-pyridin-3-ylacetamide",
        None,
    ),
    # negative guard: 同种杂环作取代基/母体时不得带方括号
    (
        "CC1=CC(=NO1)C(=O)NC=1SC=C(N1)C=1C(OC2=CC=CC=C2C1)=O",
        "5-methyl-N-[4-(2-oxochromen-3-yl)-1,3-thiazol-2-yl]-1,2-oxazole-3-carboxamide",
        None,
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_fused_component_locant_bracket(smiles: str, en: str, zh: str | None) -> None:
    """附加组分结构位次在稠合名中带方括号，取代基/母体位置不带。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
