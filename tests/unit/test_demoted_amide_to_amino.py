# IUPAC: P-41 / P-65.1.2
# Layer: L1,L3,L5
"""链端酰胺被羧酸（主基团）挤掉后，降级为 oxo + amino 前缀表达，
且 amino/hydroxy 是简单单核前缀，不得带括号。

回归：NC(=O)… 端酰胺 + COOH 主链，曾错误输出 5-(amino)-…（递归 claim
的 paren 规则误伤）；修复 = L1 _arbitrate_parts 把降级酰胺的伯 N 回收
为 amines，走正规 amino 前缀通道。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # 用户例（酮式）：HOOC-CH(OH)-CH2-CH2-C(=O)NH2
    ("NC(=O)CCC(O)C(=O)O", "5-amino-2-hydroxy-5-oxopentanoic acid", "5-氨基-2-羟基-5-氧代戊酸"),
    # 用户例（烯醇/亚胺酸写法，互变异构后应同名）
    ("N=C(O)CCC(O)C(=O)O", "5-amino-2-hydroxy-5-oxopentanoic acid", "5-氨基-2-羟基-5-氧代戊酸"),
    # 无 α-OH 变体
    ("NC(=O)CCCC(=O)O", "5-amino-5-oxopentanoic acid", None),
    ("NC(=O)CCC(=O)O", "4-amino-4-oxobutanoic acid", "4-氨基-4-氧代丁酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_demoted_amide_n_as_amino_no_paren(smiles: str, en: str, zh: str | None) -> None:
    """端酰胺被羧酸挤掉 → amino 走正规前缀，不带括号。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
