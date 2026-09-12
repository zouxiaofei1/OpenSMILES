# IUPAC: P-45 / P-62 / P-66 (N- 与 C 数字位次并入同一乘数前缀)
# Layer: L5 (assembler_prefixes)
"""N 型取代基与 C 型取代基同词干混合组：N 成员 locant 必须渲染为 N（而非 0/羰基碳位/环连接位假数字）。

触发：母体单胺/单酰胺，某词干同时出现在 C 位与该 N 上；合并组落入数字前缀通道，
N 成员的假 locant（胺 N=0、酰胺 N 被 remap 到羰基碳/环连接位）被当普通数字打印。
期望形如 N,6,6-trimethyl / N,N,2-trimethyl / N,2-dimethyl / N,3-dimethyl（N 在最前）。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer

# 干净全串（双语由修后输出核对，en 以 IUPAC/gold 为准）
FULL_CASES = [
    # 酰胺：C2-methyl + N,N-二甲基 → N,N,2-trimethyl（曾输出 1,1,2-）
    ("CC(C)C(=O)N(C)C", "N,N,2-trimethylpropanamide"),
    # 伯胺：C2-methyl + N,N-二甲基 → N,N,2-trimethyl（曾输出 0,0,2-）
    ("CC(C)CN(C)C", "N,N,2-trimethylpropan-1-amine"),
]

# 真实 A9 数据（分子仍夹其它未修缺陷：烯炔 (2E)、稠环位次、bis/di、N- 括号等），
# 仅断言本域 N-前缀段正确、无假数字。
SUBSTR_CASES = [
    # must_contain / must_not_contain（在 r.en 中）
    ("CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12.Cl",
     "N,6,6-trimethyl", "0,6,6"),
    ("COc1ccc2c(c1)N(C[C@H](C)CN(C)C)c1ccccc1S2",
     "N,N,2-trimethylpropan-1-amine", "0,0,2"),
    ("Cc1ccccc1-c1cc(N2CCN(C)CC2)ncc1N(C)C(=O)C(C)(C)c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
     "N,2-dimethyl", "1,2-dimethyl"),
    ("Cl.CN(C(=O)C=1SC=CC1C)C1CCNCC1",
     "N,3-dimethyl", "2,3-dimethyl"),
]


@pytest.mark.parametrize("smiles,en", FULL_CASES)
def test_mixed_n_stem_full(smiles: str, en: str) -> None:
    """纯胺/酰胺混合组整串命名：N 在数字前、无假数字。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,must,mustnot", SUBSTR_CASES)
def test_mixed_n_stem_real(smiles: str, must: str, mustnot: str) -> None:
    """真实近错分子（含旁系缺陷）仅校验 N-前缀段与假数字消失。"""
    en = SMILESNNamer().name(smiles).en
    assert must in en
    assert mustnot not in en
