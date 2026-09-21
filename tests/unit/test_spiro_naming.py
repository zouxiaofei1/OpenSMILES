# IUPAC: P-24.1 / P-24.2.1 / P-24.2.2 / P-24.2.3 / P-24.5 / P-24.6 / P-24.7
# Layer: L1 + L2 + L4 + L5
#
# 螺环母体名：自由螺原子判定 + 环组分合并 + von Baeyer 螺描述符 + 烷词干。
# 组分式螺环：多环环组分划分 + 引用序（字母数字序 / 端→心→端）+ 带撇位次。
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.layer1.analyzer import analyze
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer2.parent_select import select_parent
from namepredict.layer2.spiro_system import FbsComponent, citation_order, component_indices
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
# 含多环环组分的螺环（P-24.5）：走组分式命名，见文末组分式螺环用例
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


def test_fused_bridged_spiro_uses_component_style():
    """P-24.5：多环组分不再走 spiro[x.y]，改出组分式名（保留名 + 带撇位次）。"""
    r = SMILESNNamer().name(FUSED_BRIDGED_SPIRO)
    assert r.success
    assert normalize_en(r.en) == normalize_en("spiro[2-benzofuran-3,9'-xanthene]-1-one")


# ── 组分式螺环（P-24.5~24.7） ───────────────────────────

S, D = Chem.BondType.SINGLE, Chem.BondType.DOUBLE


def _smiles(bonds, elements=None) -> str:
    """由键表构造分子，返回规范 SMILES（原子序即键表下标）。"""
    n = max(max(i, j) for i, j, _ in bonds) + 1
    rw = Chem.RWMol()
    for z in (elements or ["C"] * n):
        rw.AddAtom(Chem.Atom(z))
    for i, j, b in bonds:
        rw.AddBond(i, j, b)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


# 直链二螺：环戊烷 – 环己烷（螺原子位于 1,4） – 茚
LINEAR_DISPIRO = [
    (0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 0, S),
    (0, 5, S), (5, 6, S), (6, 7, S), (7, 8, S), (8, 9, S), (9, 0, S),
    (7, 10, S), (10, 11, D), (11, 12, S), (12, 13, S), (13, 7, S),
    (12, 14, D), (14, 15, S), (15, 16, D), (16, 17, S), (17, 13, D),
]
# 支链三螺：环己烷中心组分的三个螺原子（1,3,5）分接环戊烷 / 茚 / 氧杂环己烷
BRANCHED_TRISPIRO = [
    (0, 1, S), (1, 2, S), (2, 3, S), (3, 4, S), (4, 5, S), (5, 0, S),
    (0, 6, S), (6, 7, S), (7, 8, S), (8, 9, S), (9, 0, S),
    (2, 10, S), (10, 11, D), (11, 12, S), (12, 13, S), (13, 2, S),
    (12, 14, D), (14, 15, S), (15, 16, D), (16, 17, S), (17, 13, D),
    (4, 18, S), (18, 19, S), (19, 20, S), (20, 21, S), (21, 22, S), (22, 4, S),
]


# IUPAC 原文 P-24.5.1 的 PIN 例（两个环组分，按字母数字序定谁列首）
P245_PIN = [
    ("c1ccc2c(c1)C1(CCNCC1)c1ccccc1O2", "spiro[piperidine-4,9'-xanthene]"),
    ("c1ccc2c(c1)C3(CCCCC3)C=C2", "spiro[cyclohexane-1,1'-1H-indene]"),
]


def _comp(index: int, name: str, spiros=(), extra=()) -> FbsComponent:
    """构造仅用于排序测试的合成组分。"""
    return FbsComponent(index, "mono_ring", (), (), tuple(spiros), name, name,
                        None, None, None, None, (), tuple(extra))


@pytest.mark.parametrize("smiles,en", P245_PIN)
def test_iupac_p245_examples(smiles, en):
    """P-24.5.1：两环组分按字母数字序列出，第二组分位次带撇号。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)


def test_iupac_p246_linear_dispiro():
    """P-24.6：直链多螺从字母序较小的端组分起，位次对夹在组分名之间。"""
    r = SMILESNNamer().name(_smiles(LINEAR_DISPIRO))
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(
        "dispiro[cyclopentane-1,1'-cyclohexane-4',1''-1H-indene]")


def test_iupac_p247_branched_trispiro():
    """P-24.7.2：端→心→端；首列端组分接中心的低位次螺原子。"""
    r = SMILESNNamer().name(_smiles(BRANCHED_TRISPIRO, ["C"] * 20 + ["O", "C", "C"]))
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(
        "trispiro[cyclopentane-1,1'-cyclohexane-3',1''-1H-indene-5',4'''-oxane]")


# 金标（benchmarks/spiro_benchmark.json）逐字吻合的组分式螺环
GOLD = [
    ("O=C1OC2(c3ccc(N=C=S)cc31)c1cc(Br)c(O)c(Br)c1Oc1c2cc(Br)c(O)c1Br",
     "2',4',5',7'-tetrabromo-3',6'-dihydroxy-6-isothiocyanatospiro"
     "[2-benzofuran-3,9'-xanthene]-1-one"),
]


@pytest.mark.parametrize("smiles,en", GOLD)
def test_gold_component_style(smiles, en):
    """保留名稠环组分（2-benzofuran / xanthene）+ 带撇取代基位次。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)


def test_citation_order_alphabetical():
    """P-24.5.1：两组分时字母数字序在前的列首（与资历无关）。"""
    comps = [_comp(0, "xanthene", (0,)), _comp(1, "piperidine", (0,))]
    edges = [[(0, 1)], [(0, 0)]]
    assert citation_order(comps, edges) == ((1,), (0,))


def test_citation_order_descriptor_tie_break():
    """P-24.5.3：字母相同时按 von Baeyer 描述符数字取小。"""
    comps = [_comp(0, "bicyclo[3.2.1]octane", (0,), (3, 2, 1)),
             _comp(1, "bicyclo[2.2.2]octane", (0,), (2, 2, 2))]
    edges = [[(0, 1)], [(0, 0)]]
    assert citation_order(comps, edges) == ((1,), (0,))


def test_citation_order_linear_starts_at_smaller_terminal():
    """P-24.6：直链从字母序较小的端组分起，沿链展开引用序。"""
    comps = [_comp(0, "indene", (0,)), _comp(1, "cyclohexane", (0, 1)),
             _comp(2, "cyclopentane", (1,))]
    edges = [[(0, 1)], [(0, 0), (1, 2)], [(1, 1)]]
    assert citation_order(comps, edges) == ((2,), (1,), (0,))


def test_citation_order_branched_terminal_center_terminal():
    """P-24.7.2：字母序最早端环 → 中心组分 → 其余端环按字母序。"""
    comps = [_comp(0, "oxane", (0,)), _comp(1, "cyclohexane", (0, 1, 2)),
             _comp(2, "indene", (0,)), _comp(3, "cyclopentane", (2,))]
    edges = [[(0, 1)], [(0, 0), (1, 2), (2, 3)], [(1, 1)], [(2, 1)]]
    assert citation_order(comps, edges) == ((3,), (1,), (2,), (0,))


def test_citation_order_groups_identical_terminals():
    """P-24.7.2：同名端环合并为一个引用项（出 bis/tris 倍增词）。"""
    comps = [_comp(0, "bicyclo[1.1.0]butane", (0,)), _comp(1, "cyclopropane", (0, 1, 2)),
             _comp(2, "cyclopropane", (0,)), _comp(3, "cyclopropane", (2,))]
    edges = [[(0, 1)], [(0, 0), (1, 2), (2, 3)], [(1, 1)], [(2, 1)]]
    assert citation_order(comps, edges) == ((0,), (1,), (2, 3))


def test_citation_order_all_terminals_identical_center_first():
    """P-24.7.1：端环全同时中心组分局首，端环整体作一个 tris(...) 组。"""
    comps = [_comp(0, "cyclopropane", (0,)), _comp(1, "cyclohexane", (0, 1, 2)),
             _comp(2, "cyclopropane", (0,)), _comp(3, "cyclopropane", (2,))]
    edges = [[(0, 1)], [(0, 0), (1, 2), (2, 3)], [(1, 1)], [(2, 1)]]
    assert citation_order(comps, edges) == ((1,), (0, 2, 3))


def test_iupac_p247_bis_group():
    """P-24.7.2：同名端环用 bis(...)，组内位次对以冒号连接、各占一个撇号层级。"""
    r = SMILESNNamer().name("C12(C3CC13)C1(CC1)C12CC1")
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(
        "trispiro[bicyclo[1.1.0]butane-2,1'-cyclopropane-2',1'':3',1'''-bis(cyclopropane)]")


def test_component_split_keeps_non_free_spiro_together():
    """非自由螺原子不切断组分；自由螺原子处才切（P-24.1）。"""
    mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCC2")
    sys_ = next(s for s in build_ring_systems(mol) if s["free_spiro_atoms"])
    assert component_indices(sys_) == [[0], [1]]


def test_multi_component_system_is_fused_bridged_spiro():
    """含多环组分的螺环系身份取 fused_bridged_spiro。"""
    mol = Chem.MolFromSmiles(_smiles(LINEAR_DISPIRO))
    ids = {p.get("scaffold_id") for p in select_parent(analyze(mol))}
    assert "fused_bridged_spiro" in ids
    assert "fused_hetero" not in ids  # 不得回退稠环身份
