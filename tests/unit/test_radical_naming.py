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


# 杂原子锚点自由基：单核氢化物母体（P-15.4.1 表 2.1 oxidane/azane/sulfane）
# + 烷基取代基 → free_to_yl 转标准名；不走碳链 -yl。
HETERO_RADICAL_CASES = [
    ("*OCC", "ethyloxy", "乙氧基"),
    ("*OCCC", "propyloxy", "丙氧基"),
    ("*NCC", "ethylamino", "乙氨基"),
    ("*NCCC", "propylamino", "丙氨基"),
    ("*OCc1ccccc1", "benzyloxy", "苄氧基"),
    ("*Oc1ccccc1", "phenyloxy", "苯氧基"),
    ("*NC(=O)C", "acetylamino", "乙酰氨基"),
    ("*S(=O)(=O)C(C)C", "isopropylsulfanyl", "异丙硫基"),
]


@pytest.mark.parametrize("smiles,en,zh", HETERO_RADICAL_CASES)
def test_hetero_radical_mononuclear(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 取代基递归:表外变体经 radical worker 拼出带位次的 yl 名,不再 free_to_yl 二次转换。
SIDE_CASES = [
    # IUPAC 低位次:CH2 连 N 邻位(2),锚定 radical 自洽输出 4-methylpyridin-2-yl;6-yl 是 free_to_yl canonical-rank 反推的旧错误。
    ("O=C(O)Cc1cc(C)ccn1", "2-(4-methylpyridin-2-yl)acetic acid", "2-(4-甲基吡啶-2-基)乙酸"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    # 碳连接点锚定优先:4-氯苯基(原 _arene_yl_from_sub 芳基短路,锚定改造后仍正确)。
    ("O=C(O)Cc1ccc(Cl)cc1", "2-(4-chlorophenyl)acetic acid", "2-(4-氯苯基)乙酸"),
    # 改进:HOCH2CH2- 侧链不再被 free_to_yl 误前缀化成 ethoxy,锚定 radical 给 2-hydroxyethan-1-yl。
    ("OCCc1ccc(C(=O)O)cc1", "4-(2-hydroxyethan-1-yl)benzoic acid", "4-(2-羟基乙-1-基)苯甲酸"),
    # 非碳连接点(醚氧)仍走 H 封端 + free_to_yl:乙氧基不变。
    ("CCOc1ccccc1", "ethoxybenzene", "乙氧基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", SIDE_CASES)
def test_radical_side_recursive(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
