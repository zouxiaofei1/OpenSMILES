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
# 自由基链自由价在 1 位时省略位次（P-29.2 方法 1）：-ethan-1-yl → -ethyl。
# 组成成员（伯酰胺 N → amino、酯烷氧基 → methoxy）由 L1 回收进正规前缀，
# 简单单核前缀不带括号。
# 注：*CC(=O)O 命中 anchored 保留名 carboxymethyl（非降级路径）；*CCC(=O)O 的游离
# 中性 COOH 被 radical 压制后保持羧酸叶身份作 carboxy 前缀（P-61.1.3），不再降级成
# 同碳 oxo+hydroxy 也不占链长（A4 修复，gold 2-carboxyethyl）。
POSITIVE = [
    ("C(*)C(N)=O", "2-amino-2-oxoethyl", "2-氨基-2-氧代乙基"),
    ("C(*)C=O", "2-oxoethyl", "2-氧代乙基"),
    ("*CCC(=O)O", "2-carboxyethyl", "2-羧基乙基"),
    ("C(*)C(=O)OC", "2-methoxy-2-oxoethyl", "2-甲氧基-2-氧代乙基"),
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
    ("C(*)CO", "2-hydroxyethyl", "2-羟基乙基"),
    ("C(*)C(C)=O", "2-oxopropyl", "2-氧代丙基"),
]


@pytest.mark.parametrize("smiles,en,zh", NEGATIVE)
def test_composite_fg_principal_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
