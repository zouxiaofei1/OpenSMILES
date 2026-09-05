# IUPAC: P-65.1 / P-29.3；Layer: L2,L3,L5
"""N-烷基取代氨基醇的递归命名：仲胺带链外烷基侧链时不认领为
"amino"，整条侧链交由递归后端命名（防截断成 aminoethanol）。"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")
CASES = [
    # 断点 case：仲胺 N 同时连链碳与侧链碳（深度递归两层）
    (
        "OCCNCCNCCC",
        "2-[2-(propylamino)ethylamino]ethanol",
        "2-[2-(丙氨基)乙氨基]乙醇",
    ),
    # 同一分子反转输入：结果必须一致（此前依赖输入原子序）
    (
        "CCCNCCNCCO",
        "2-[2-(propylamino)ethylamino]ethanol",
        "2-[2-(丙氨基)乙氨基]乙醇",
    ),
    # N-甲基 / N-氨基乙基（自由基链 1 位省略：2-aminoethylamino）
    ("OCCNCCN", "2-(2-aminoethylamino)ethanol", "2-(2-氨基乙氨基)乙醇"),
    # N-苄基
    ("OCCNCc1ccccc1", "2-(benzylamino)ethanol", "2-(苄氨基)乙醇"),
    # N 桥接二醇（母体链不含 N，两侧都是侧链）
    ("OCCNCCO", "2-(2-hydroxyethylamino)ethanol", "2-(2-羟基乙氨基)乙醇"),
    # 一级胺取代乙醇：C1 有可取代 H，2- 不可省略（P-14.3.4.4）
    ("NCCO", "2-aminoethanol", "2-氨基乙醇"),
    # 大分子保持多层递归（饱和链自由价 1 位省略：ethylamino/pentylamino）
    (
        "NCCNCCCCCNCCNCCNCCCNCCCNCCO",
        "2-[3-[3-[2-[2-[5-(2-aminoethylamino)pentylamino]"
        "ethylamino]ethylamino]propylamino]"
        "propylamino]ethanol",
        "2-[3-[3-[2-[2-[5-(2-氨基乙氨基)戊氨基]"
        "乙氨基]乙氨基]丙氨基]丙氨基]乙醇",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_alkyl_amino_alcohol(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
