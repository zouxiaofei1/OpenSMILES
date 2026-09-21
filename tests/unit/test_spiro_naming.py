# IUPAC: P-24.1 / P-24.2.1 / P-24.2.2 / P-24.2.3
# Layer: L1 + L2 + L4 + L5
#
# 螺环母体名：自由螺原子判定 + 环组分合并 + von Baeyer 螺描述符 + 烷词干。
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.layer1.analyzer import analyze
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer2.parent_select import select_parent
from namepredict.layer5.spiro_namer import spiro_descriptor_str, spiro_multiplier
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# 每条 = (SMILES, 期望英文名, 期望中文名)。金标旁证见 benchmarks/spiro_benchmark.json。
CASES = [
    ("C1CCC2(CC1)CCCC2", "spiro[4.5]decane", "螺[4.5]癸烷"),
    ("C1CCC12CCCC2", "spiro[3.4]octane", "螺[3.4]辛烷"),
    ("C1CCC2(CC1)CCCCC2", "spiro[5.5]undecane", "螺[5.5]十一烷"),
    ("C1CNC12CCCC2", "1-azaspiro[3.4]octane", "1-氮杂螺[3.4]辛烷"),
    ("C1CCC12OCC2", "1-oxaspiro[3.3]heptane", "1-氧杂螺[3.3]庚烷"),
    ("C1CN(C2(CC2)CC1)S(=O)(=O)C", "4-methylsulfonyl-4-azaspiro[2.5]octane",
     "4-甲磺酰基-4-氮杂螺[2.5]辛烷"),
]

# IUPAC 原文 P-24.2.2 / P-24.2.3 的 PIN 例子（上标按判分口径写平位数字）。
# 每条 = (SMILES, 期望英文名)。SMILES 按描述符的段序反构。
IUPAC_POLYSPIRO = [
    ("C1CCC12CCC3(CCC3)CC2", "dispiro[3.2.37.24]dodecane"),        # 线性二螺
    ("C1CC12CCC1(CC2)CCC2(CC2)CC1", "trispiro[2.2.2.29.26.23]pentadecane"),   # 线性三螺
    ("C1CC12CCC1(CC1)CCC1(CC1)CC2", "trispiro[2.2.26.2.211.23]pentadecane"),  # 支链三螺
]

# 非自由螺连接：两环除共用原子外还有桥相连（P-24.1 排除，归 P-23）
NON_FREE_SPIRO = "C1CC23CC(C2)C13"
# 稠合+螺环（P-24.5 组分式命名）：本版未实现，须显式失败而非静默错名
FUSED_BRIDGED_SPIRO = "O=C1OC2(c3ccccc31)c1ccccc1Oc1ccccc12"


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_spiro_parent_name(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en", IUPAC_POLYSPIRO)
def test_iupac_polyspiro_examples(smiles, en):
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)


def test_multi_spiro_needs_superscripts():
    """多螺描述符带重访螺原子位次上标，单螺不带（P-24.2.2 vs P-24.2.1）。"""
    assert spiro_descriptor_str((2, 2, 2, 2), (0, 0, 0, 0)) == "2.2.2.2"
    assert spiro_descriptor_str((3, 2, 3, 2), (0, 0, 7, 4)) == "3.2.37.24"


def test_spiro_multiplier_keeps_spiro_stem():
    """词头：1 用 spiro，2 起用 dispiro（不是 bi/di 裸计数词）。"""
    assert spiro_multiplier(1) == ("spiro", "螺")
    assert spiro_multiplier(2) == ("dispiro", "二螺")
    assert spiro_multiplier(3) == ("trispiro", "三螺")


def test_ring_system_merges_free_spiro():
    """L1：螺连环对并入同一环系，并给出 spiro_edges 与 free_spiro_atoms。"""
    mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCC2")
    systems = build_ring_systems(mol)
    assert len(systems) == 1  # 两环合并为一个环系
    sys_ = systems[0]
    assert sys_["topology"] == "spiro"
    assert sys_["n_rings"] == 2
    assert sys_["free_spiro_atoms"] == [3]
    assert [i for _, _, i in sys_["spiro_edges"]] == [3]
    assert sys_["fusion_edges"] == []  # 螺连接不是稠合边


def test_non_free_spiro_union_not_free():
    """两环除共用原子外还有桥相连：该原子不算自由螺原子（P-24.1）。"""
    systems = build_ring_systems(Chem.MolFromSmiles(NON_FREE_SPIRO))
    assert all(not s["free_spiro_atoms"] for s in systems)
    r = SMILESNNamer().name(NON_FREE_SPIRO)
    assert r.success and "spiro[" not in r.en  # 走 P-23 桥环，不得出现螺环名


def test_non_spiro_ring_systems_unchanged():
    """无螺环的分子不得被误标为螺环系。"""
    for smi in ("c1ccc2ccccc2c1", "c1ccccc1-c1ccccc1", "C1C2CC3CC1CC(C2)C3"):
        for s in build_ring_systems(Chem.MolFromSmiles(smi)):
            assert not s["free_spiro_atoms"]
            assert s["topology"] is None


def test_scaffold_id_is_mono_spiro():
    """L2：自由螺原子为 1 且各环组分均为单环时身份取 mono_spiro。"""
    mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCC2")
    ids = {p.get("scaffold_id") for p in select_parent(analyze(mol))}
    assert "mono_spiro" in ids
    assert "fused_hetero" not in ids  # 不得回退稠环身份


def test_mono_spiro_needs_only_monocyclic_components():
    """中心环含两个螺原子仍属 mono_spiro（全单环组分）；含稠环组分则不是。"""
    poly = Chem.MolFromSmiles("C1CC2(CC1)CCC1(CC2)CC1")
    assert {p.get("scaffold_id") for p in select_parent(analyze(poly))} == {"mono_spiro"}
    fused = Chem.MolFromSmiles(FUSED_BRIDGED_SPIRO)
    assert {p.get("scaffold_id") for p in select_parent(analyze(fused))} == \
        {"fused_bridged_spiro"}


def test_fused_bridged_spiro_fails_explicitly():
    """P-24.5 组分式螺环未实现：须显式失败，不得给静默错名。"""
    r = SMILESNNamer().name(FUSED_BRIDGED_SPIRO)
    assert not r.success
    assert r.en == ""
