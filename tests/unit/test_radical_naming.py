# IUPAC: P-22.2.4 / P-29.3
# Layer: L2,L3,L4,L5
"""自由基主基团统一命名 worker:开链/杂环/苯经 _KIND_TABLE["radical"] 拼 -yl 词干。"""

from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# 独立自由基命名(带 * 锚点):radical 主基团走统一 worker,输出带位次的 -yl 名。
RADICAL_CASES = [
    ("*c1ccccc1", "phenyl", "苯基"),
    ("*c1ccc(Cl)cc1", "4-chlorophenyl", "4-氯苯基"),
    ("*c1ccncc1", "pyridin-4-yl", "吡啶-4-基"),
    ("*c1cc(C)ccn1", "4-methylpyridin-2-yl", "4-甲基吡啶-2-基"),
    ("*c1ccc2ccccc2c1", "naphthalen-1-yl", "萘-1-基"),
    ("*CCCCCC", "hexan-1-yl", "己-1-基"),
    ("*CC(C)C", "2-methylpropan-1-yl", "2-甲基丙-1-基"),
]


@pytest.mark.parametrize("smiles,en,zh", RADICAL_CASES)
def test_radical_stem(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 取代基递归:表外变体经 radical worker 拼出带位次的 yl 名,不再 free_to_yl 二次转换。
SIDE_CASES = [
    ("O=C(O)Cc1cc(C)ccn1", "2-(4-methylpyridin-6-yl)acetic acid", "2-(4-甲基吡啶-6-基)乙酸"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", SIDE_CASES)
def test_radical_side_recursive(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
