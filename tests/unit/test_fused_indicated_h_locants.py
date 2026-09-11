# IUPAC: P-25.3.3.1.2
# Layer: L4
"""稠环指示氢位次（P-25.3.3.1.2(f)）与镜像方向 CIP 破局（P-14.4(j)）。

Affiliated: P-58.2.1.1 指示氢定义、P-14.4(e)(i) hydro 前缀位次、P-32.2.3 位次优先级顺序。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: (a)-(d) 全平局、指示氢环位集合是唯一判据 —— 集合最小者胜（甲基落 5 位而非 2 位，
    # 即把最低位次给带氢环位、而不是给取代基），甲基在 5 位的方向还须满足 P-14.4(j) R 取较低位次。
    ("CN1C[C@@H]2CNC[C@@H]2C1",
     "(3aR,6aS)-5-methyl-2,3,3a,4,6,6a-hexahydro-1H-pyrrolo[3,4-c]pyrrole",
     "(3aR,6aS)-5-甲基-2,3,3a,4,6,6a-六氢-1H-吡咯并[3,4-c]吡咯"),
    ("C=C(C)[C@H]1CC[C@@]2(C)CCC=C(C)[C@@H]2C1",
     "(3S,4aR,8aR)-5,8a-dimethyl-3-prop-1-en-2-yl-2,3,4,4a,7,8-hexahydro-1H-naphthalene",
     "(3S,4aR,8aR)-5,8a-二甲基-3-丙-1-烯-2-基-2,3,4,4a,7,8-六氢-1H-萘"),
    ("C=C(C)[C@@H]1CCC2=CCC[C@@H](C)[C@@]2(C)C1",
     "(3R,4aR,5R)-4a,5-dimethyl-3-prop-1-en-2-yl-2,3,4,5,6,7-hexahydro-1H-naphthalene",
     "(3R,4aR,5R)-4a,5-二甲基-3-丙-1-烯-2-基-2,3,4,5,6,7-六氢-1H-萘"),
    # negative: 同类骨架（六氢/二氢-1H-稠环）但含主特征基团或可被取代基位次定妥，指示氢层不得改写
    ("CC1(CCCC=2CCC(CC12)C=O)C",
     "8,8-dimethyl-2,3,4,5,6,7-hexahydro-1H-naphthalene-2-carbaldehyde",
     "8,8-二甲基-2,3,4,5,6,7-六氢-1H-萘-2-甲醛"),
    ("Cc1ccc2c(c1)CCCC2(C)C",
     "4,4,7-trimethyl-2,3-dihydro-1H-naphthalene",
     "4,4,7-三甲基-2,3-二氢-1H-萘"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indicated_hydrogen_locant_set(smiles: str, en: str, zh: str) -> None:
    """指示氢环位集合最小者取编号：正例改判到带氢位次更低的走向，负例名不变。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
