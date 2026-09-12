# IUPAC: P-25.3.2 / P-25.3.8
# Layer: L5
"""fused_namer: 未注册稠环的稠合名称组装(benzo[a].../naphtho[...]... 类)。"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # tetracene(并四苯, 未注册) → 蒽 + 苯稠合
    ("c1ccc2cc3cc4ccccc4cc3cc2c1", "benzo[b]anthracene", "苯并[b]蒽"),
    # benchmark 参考: benzo[e]pyrene(未注册) → 芘 + 苯稠合
    ("c1ccc2c(c1)c1cccc3ccc4cccc2c4c31", "benzo[e]pyrene", "苯并[e]芘"),
    # chrysene(未注册) → 菲 + 苯稠合
    ("c1cc2c3ccccc3ccc2c2ccccc12", "benzo[a]phenanthrene", "苯并[a]菲"),
    # 杂单环稠合: pyrimidine(母体) + pyridine(附加)。碱环是嘧啶时需让稠合边落在
    # 其 C4-C5(d 侧)而非 C5-C6(e 侧)——组分编号的 locant 1 在等价双 N 间浮动后
    # 由稠合原子位次最小化决定, 描述符才对(P-25.3.1.3 位次尽可能低)。
    ("ClC=1C2=C(N=C(N1)C)N=CC=C2", "4-chloro-2-methylpyrido[2,3-d]pyrimidine", "4-氯-2-甲基吡啶并[2,3-d]嘧啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_fused_namer(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
