# IUPAC: P-24.5 / P-24.6 / P-24.7
# Layer: L1 + L2 + L4 + L5
#
# 组分式螺环母体名：多环环组分划分 + 引用序（字母数字序 / 端→心→端）+ 带撇位次。
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.layer1.analyzer import analyze
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer2.fbs_system import FbsComponent, citation_order, component_indices
from namepredict.layer2.parent_select import select_parent
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en

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
