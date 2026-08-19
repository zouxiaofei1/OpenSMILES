# IUPAC: P-25.3.2 / P-25.3.8
# Layer: L5
"""fused_namer: 未注册稠环的稠合名称组装(benzo[a].../naphtho[...]... 类)。"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # tetracene(并四苯, 未注册) → 蒽 + 苯稠合
    ("c1ccc2cc3cc4ccccc4cc3cc2c1", "benzo[b]anthracene", "苯并[b]蒽"),
    # benchmark 参考: benzo[e]pyrene(未注册) → 芘 + 苯稠合
    ("c1ccc2c(c1)c1cccc3ccc4cccc2c4c31", "benzo[e]pyrene", "苯并[e]芘"),
    # chrysene(未注册) → 菲 + 苯稠合
    ("c1cc2c3ccccc3ccc2c2ccccc12", "benzo[b]phenanthrene", "苯并[b]菲"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_fused_namer(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
