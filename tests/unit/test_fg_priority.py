# IUPAC: P-41
# Layer: L1
"""组合官能团主基团优先级仲裁：更高优先级 FG（如自由基）存在时，
组合羰基 FG（酰胺/醛/羧酸/酯）退出主基团并降级为 oxo 前缀，
不再丢失羰基氧。

P-41 规定 radical 的类优先级高于酸/酯/酰胺/醛；非最高组合 FG 以
前缀表达（oxo/amino/…），而不是被丢弃。
"""

from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# 正例：radical 主基团抢占后，组合羰基 FG 降级为 oxo（+ 组成基团）前缀。
POSITIVE = [
    ("C(*)C(N)=O", "2-(amino)-2-oxoethan-1-yl", "2-(氨基)-2-氧代乙-1-基"),
    ("C(*)C=O", "2-oxoethan-1-yl", "2-氧代乙-1-基"),
    ("C(*)C(O)=O", "2-(hydroxy)-2-oxoethan-1-yl", "2-(羟基)-2-氧代乙-1-基"),
    ("C(*)C(=O)OC", "2-methoxy-2-oxoethan-1-yl", "2-甲氧基-2-氧代乙-1-基"),
]


@pytest.mark.parametrize("smiles,en,zh", POSITIVE)
def test_composite_fg_degrade_to_oxo(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 负例：组合 FG 自身是最高优先级时，保持主基团命名不变。
NEGATIVE = [
    ("CC(N)=O", "acetamide", "乙酰胺"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCC(O)=O", "propanoic acid", "丙酸"),
    ("CCC(N)=O", "propanamide", "丙酰胺"),
    ("C(*)CO", "2-hydroxyethan-1-yl", "2-羟基乙-1-基"),
    ("C(*)C(C)=O", "2-oxopropan-1-yl", "2-氧代丙-1-基"),
]


@pytest.mark.parametrize("smiles,en,zh", NEGATIVE)
def test_composite_fg_principal_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
