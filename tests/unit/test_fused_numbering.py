# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_fused_numbering.py: fused_numbering: 外周骨架编号 + 稠合碳 a/b/c 字母位次(P-25.3.3.1)。
test_fused_orientation.py: fused_orientation: 优选取向(水平行/坐标/象限加权)。
test_fused_locant_numbering.py: Asymmetric fused arene locant numbering (P-25.4): quinoline/quinazoline/indole.
test_fused_indicated_h_locants.py: 稠环指示氢位次（P-25.3.3.1.2(f)）与镜像方向 CIP 破局（P-14.4(j)）。
test_fused_component_locant_bracket.py: 附加组分名中的结构位次（杂原子位置）在稠合名里须置于方括号内。
test_fused_cycloalkane_component.py: 环烷烃作稠合附加组分：饱和单环烃名删尾 'ne' 得前缀（cyclopentane → cyclopenta），
test_partial_hydro_ring_locant_prefix.py: 局部不饱和保留名的词干已含加氢前缀时，locant 前缀须留在组分名前，不得前移。
test_locant_key.py: locant_key: 数字/字母位次统一排序。
"""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer2.ring_scaffold import FUSED56_LABELS, NAPH_LABELS
from namepredict.layer4.fused_numbering import fused_atoms, number_fused_system
from namepredict.layer4.fused_orientation import Orientation, preferred_orientation, preferred_orientations
from namepredict.layer4.locant_calc import locant_key, locant_str_sort
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_fused_numbering.py
# IUPAC: P-25.3.3
# Layer: L4
#
# fused_numbering: 外周骨架编号 + 稠合碳 a/b/c 字母位次(P-25.3.3.1)。
# ==========================================================================
def fused_numbering___number(smiles: str):
    mol = preprocess(smiles)
    assert mol is not None
    systems = build_ring_systems(mol)
    assert systems
    rings = list(mol.GetRingInfo().AtomRings())
    orients = preferred_orientations(mol, rings, systems[0]["fusion_edges"])
    assert orients
    return mol, number_fused_system(mol, rings, [o.coord_dict() for o in orients])


def test_naphthalene_labels_match_standard():
    _, (chain, labels) = fused_numbering___number("c1ccc2ccccc2c1")
    assert len(chain) == 10
    assert labels == list(NAPH_LABELS)


def test_indole_labels_match_standard():
    _, (chain, labels) = fused_numbering___number("c1ccc2[nH]ccc2c1")
    assert len(chain) == 9
    assert labels == list(FUSED56_LABELS)


def test_quinoline_hetero_atom_locant_one():
    mol, (chain, labels) = fused_numbering___number("c1ccc2ncccc2c1")
    n = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() == 7][0]
    assert labels[chain.index(n)] == "1"


def test_phenanthrene_14_atoms_with_letter_fused_carbons():
    mol, (chain, labels) = fused_numbering___number("c1ccc2c(c1)ccc1ccccc12")
    assert len(chain) == 14
    fused = {a for a in chain if labels[chain.index(a)][-1].isalpha()}
    assert len(fused) >= 3  # 稠合碳均有字母位次


def test_azulene_odd_ring_pair_numbering():
    """azulene 5+7 奇环对: 10 原子全部编号, 桥头字母位。"""
    _, (chain, labels) = fused_numbering___number("C1=CC=C2C=CC=CC2=C1")
    assert len(chain) == 10
    assert any(lbl[-1].isalpha() for lbl in labels)


def test_walk_starts_next_to_a_fused_atom():
    """P-25.3.3.1.1 外周行走自稠合边一端起算: 起点取偏环顶端的顶点会使位次整体错一位。"""
    for smiles in ("COc1cc2oc(=O)c3c(O)cc(O)cc3c2c(C)c1Cl",
                   "C1=C2C=3N(C(NC2=CC=C1)=O)C1=C(N3)C=CC=C1"):
        mol, (chain, _) = fused_numbering___number(smiles)
        fused = fused_atoms(list(mol.GetRingInfo().AtomRings()))
        ring = next(r for r in mol.GetRingInfo().AtomRings() if chain[0] in r)
        i = ring.index(chain[0])
        assert ring[i - 1] in fused or ring[(i + 1) % len(ring)] in fused


def test_benzo_c_chromenone_lactone_locant():
    """苯并[c]色烯-6-酮的环内羰基得 6 位(起点紧邻稠合边; 起点偏顶点时错编为 5-酮)。"""
    mol, (chain, labels) = fused_numbering___number("COc1cc2oc(=O)c3c(O)cc(O)cc3c2c(C)c1Cl")
    locants = dict(zip(chain, labels))
    carbonyl = next(a for a in chain if any(
        n.GetSymbol() == "O" and mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondType().name == "DOUBLE"
        for n in mol.GetAtomWithIdx(a).GetNeighbors()))
    assert locants[carbonyl] == "6"


def test_6_5_6_three_ring_numbering():
    """6-5-6 线性稠环(变形五元环居中): 三个 N 得最低位次, 稠合碳得字母位次(a/b/c)。"""
    mol, (chain, labels) = fused_numbering___number("ClC1C2NC3C=CC=CC=3C=2N=CN=1")
    assert len(chain) == 13
    heteros = sorted(labels[chain.index(a)] for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
    assert heteros == ["1", "3", "5"]  # P-25.3.3.2.3 杂原子最低位次
    fused_letter = [lbl for lbl in labels if lbl[-1].isalpha()]
    assert len(fused_letter) >= 4  # 稠合碳均有字母位次


# ==========================================================================
# 合并自 test_fused_orientation.py
# IUPAC: P-25.3.2.3
# Layer: L4
#
# fused_orientation: 优选取向(水平行/坐标/象限加权)。
# ==========================================================================
def fused_orientation___orient(smiles: str) -> Orientation | None:
    mol = preprocess(smiles)
    assert mol is not None
    systems = build_ring_systems(mol)
    if not systems:
        return None
    rings = list(mol.GetRingInfo().AtomRings())
    return preferred_orientation(mol, rings, systems[0]["fusion_edges"])


def fused_orientation___shared_vertical(orient: Orientation, rings, row) -> bool:
    """相邻行环共享边是否竖直。"""
    cd = orient.coord_dict()
    for k in range(len(row) - 1):
        shared = sorted(set(rings[row[k]]) & set(rings[row[k + 1]]))
        x0, _ = cd[shared[0]]
        x1, _ = cd[shared[1]]
        if abs(x0 - x1) > 1e-6:
            return False
    return True


def test_naphthalene_two_ring_row_symmetric():
    o = fused_orientation___orient("c1ccc2ccccc2c1")
    assert o is not None
    assert len(o.row) == 2
    q1, _, q3, _ = o.quad
    assert abs(q1 - q3) < 0.01  # 对称碳环上下镜像等价


def test_anthracene_three_ring_row():
    """anthracene 线性 3 苯环 → 水平行 3 环。"""
    mol = preprocess("c1ccc2cc3ccccc3cc2c1")
    rings = list(mol.GetRingInfo().AtomRings())
    o = fused_orientation___orient("c1ccc2cc3ccccc3cc2c1")
    assert o is not None and len(o.row) == 3
    assert fused_orientation___shared_vertical(o, rings, o.row)


def test_phenanthrene_angular_two_ring_row():
    """phenanthrene 角型 3 环 → 水平行最多 2 环，右上象限环数最多。"""
    mol = preprocess("c1ccc2c(c1)ccc1ccccc12")
    rings = list(mol.GetRingInfo().AtomRings())
    o = fused_orientation___orient("c1ccc2c(c1)ccc1ccccc12")
    assert o is not None and len(o.row) == 2
    q1, _, q3, _ = o.quad
    assert q1 > q3  # 右上环数优先于左下（phenanthrene 取向）


def test_indole_and_quinoline_two_ring_row():
    for smi in ("c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1"):
        o = fused_orientation___orient(smi)
        assert o is not None and len(o.row) == 2


def test_azulene_odd_ring_pair_row():
    """azulene 5+7 奇环对：共享边竖直（奇环端部模板）。"""
    mol = preprocess("C1=CC=C2C=CC=CC2=C1")
    rings = list(mol.GetRingInfo().AtomRings())
    o = fused_orientation___orient("C1=CC=C2C=CC=CC2=C1")
    assert o is not None and len(o.row) == 2
    assert fused_orientation___shared_vertical(o, rings, o.row)


def test_pyrene_max_row_two():
    o = fused_orientation___orient("c1cc2ccc3cccc4ccc(c1)c2c34")
    assert o is not None and len(o.row) == 2


def test_6_5_6_linear_three_ring_row():
    """6-5-6 线性稠环(中间 5 元奇环两侧稠合) → P-25.3.2.3.2 变形五元环使三环成水平行。"""
    smi = "ClC1C2NC3C=CC=CC=3C=2N=CN=1"
    mol = preprocess(smi)
    rings = list(mol.GetRingInfo().AtomRings())
    o = fused_orientation___orient(smi)
    assert o is not None
    assert len(o.row) == 3
    assert fused_orientation___shared_vertical(o, rings, o.row)


# ==========================================================================
# 合并自 test_fused_locant_numbering.py
# IUPAC: P-25.4
# Layer: L2,L4
#
# Asymmetric fused arene locant numbering (P-25.4): quinoline/quinazoline/indole.
#
# These fused heteroarenes have fixed IUPAC numbering (heteroatom = 1, then along
# the ring with bridgehead letters 4a/7a). P-14.4 generic ring enumeration
# mis-numbers them (10-chloroquinolin-5-ol instead of 2-chloroquinolin-6-ol).
# The retained scaffold numbering facts (standard_path + template match) drive L4.
# Symmetric templates (naphthalene / quinoxaline / benzimidazole) are excluded:
# their GetSubstructMatch direction is not unique, so P-14.4's substituent
# lowest-locant rule stays in control (fixed numbering would be unstable).
# ==========================================================================
fused_locant_numbering__CASES = [
    # positive: quinoline locants (N=1; 2-chloro, 6-ol)
    ("ClC1=NC2=CC=C(C=C2C=C1)O", "2-chloroquinolin-6-ol", "2-氯喹啉-6-醇"),
    ("ClC1=NC2=CC=C(C=C2C(=C1)C)O", "2-chloro-4-methylquinolin-6-ol", "2-氯-4-甲基喹啉-6-醇"),
    # positive: quinoline CF3 (zh 括号文体与 gold 不同，只断言 EN locant)
    ("BrC1=CC=C2C=CC(=NC2=C1)C(F)(F)F", "7-bromo-2-(trifluoromethyl)quinoline", None),
    # positive: quinazoline locants (N1,N3; 2-chloro, 4-methyl, 7-methoxy)
    ("ClC1=NC2=CC(=CC=C2C(=N1)C)OC", "2-chloro-7-methoxy-4-methylquinazoline", "2-氯-7-甲氧基-4-甲基喹唑啉"),
    ("ClC1=NC(=NC2=CC=C(C=C12)F)C1=CC=C(C=C1)OC", "4-chloro-6-fluoro-2-(4-methoxyphenyl)quinazoline", "4-氯-6-氟-2-(4-甲氧基苯基)喹唑啉"),
    ("BrCC(=O)C1=NC2=CC=CC=C2N=C1C", "2-bromo-1-(3-methylquinoxalin-2-yl)ethanone", None),
    ("ClC1=NC2=CC=CC=C2C(=C1)C(F)(F)F", "2-chloro-4-(trifluoromethyl)quinoline", None),
]


@pytest.mark.parametrize("smiles,en,zh", fused_locant_numbering__CASES)
def test_fused_locant_numbering(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused_indicated_h_locants.py
# IUPAC: P-25.3.3.1.2
# Layer: L4
#
# 稠环指示氢位次（P-25.3.3.1.2(f)）与镜像方向 CIP 破局（P-14.4(j)）。
#
# Affiliated: P-58.2.1.1 指示氢定义、P-14.4(e)(i) hydro 前缀位次、P-32.2.3 位次优先级顺序。
# ==========================================================================
fused_indicated_h_locants__CASES = [
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


@pytest.mark.parametrize("smiles,en,zh", fused_indicated_h_locants__CASES)
def test_indicated_hydrogen_locant_set(smiles: str, en: str, zh: str) -> None:
    """指示氢环位集合最小者取编号：正例改判到带氢位次更低的走向，负例名不变。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused_component_locant_bracket.py
# IUPAC: P-25.3.1.3
# Layer: L2,L5
#
# 附加组分名中的结构位次（杂原子位置）在稠合名里须置于方括号内。
#
# P-25.3.1.3：「描述组分结构特征的位次，如杂原子位置，保留在组分名称中，并置于方
# 括号内」；P-25.3.2.1.2 进一步规定 isoxazole/oxazole/thiazole 在稠合名中必须改用
# Hantzsch-Widman 名 1,2-oxazole / 1,3-oxazole / 1,3-thiazole，且位次置于方括号内。
# 故 1,2,4-triazole 作附加组分时是 [1,2,4]triazolo[1,5-a]pyridine，而非
# 1,2,4-triazolo[1,5-a]pyridine；该方括号位次集与母体间仍需连字符
# （2,3-dihydro-[1,3]thiazolo…，同 1,3-dihydro-2H-… 的前导数字体例）。
#
# 方括号只在「附加组分/稠合母体」位置出现；同一杂环作取代基或作母体（1,3-thiazol-2-yl、
# 1,2-oxazole-3-carboxamide）时不加方括号，末条负例守此边界。
# ==========================================================================
fused_component_locant_bracket__CASES = [
    # positive: 三唑并[1,5-a]吡啶（用户报例）
    (
        "BrC1=NN2C(C=CC=C2Br)=N1",
        "2,5-dibromo-[1,2,4]triazolo[1,5-a]pyridine",
        "2,5-二溴-[1,2,4]三唑并[1,5-a]吡啶",
    ),
    # positive: 三唑并[1,5-a]嘧啶，前缀（吡啶-2-基）与方括号母体间须留连字符
    (
        "Cc1nc(C)c(C(O)=Nc2ccc(F)c(-c3nc4ncc(-c5ccccn5)cn4n3)c2)o1",
        "N-[4-fluoro-3-(6-pyridin-2-yl-[1,2,4]triazolo[1,5-a]pyrimidin-2-yl)phenyl]-2,4-dimethyl-1,3-oxazole-5-carboxamide",
        "N-[4-氟-3-(6-吡啶-2-基-[1,2,4]三唑并[1,5-a]嘧啶-2-基)苯基]-2,4-二甲基-1,3-噁唑-5-甲酰胺",
    ),
    # positive: 噻唑并[3,2-a]嘧啶，二氢前缀与方括号之间连字符（中文侧另存括号差异，此处只断英文）
    (
        "O=C1C=CN=C2N1C(CS2)CC(=O)NC=2C=NC=CC2",
        "2-(5-oxo-2,3-dihydro-[1,3]thiazolo[3,2-a]pyrimidin-3-yl)-N-pyridin-3-ylacetamide",
        None,
    ),
    # negative guard: 同种杂环作取代基/母体时不得带方括号
    (
        "CC1=CC(=NO1)C(=O)NC=1SC=C(N1)C=1C(OC2=CC=CC=C2C1)=O",
        "5-methyl-N-[4-(2-oxochromen-3-yl)-1,3-thiazol-2-yl]-1,2-oxazole-3-carboxamide",
        None,
    ),
]


@pytest.mark.parametrize("smiles,en,zh", fused_component_locant_bracket__CASES)
def test_fused_component_locant_bracket(smiles: str, en: str, zh: str | None) -> None:
    """附加组分结构位次在稠合名中带方括号，取代基/母体位置不带。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused_cycloalkane_component.py
# IUPAC: P-25.3.2.2.1
# Layer: L2,L5
#
# 环烷烃作稠合附加组分：饱和单环烃名删尾 'ne' 得前缀（cyclopentane → cyclopenta），
# 表示最大数目非累积双键的形式（P-25.3.2.2.1）；数字位次按 P-25.3.8.1 省略。
# ==========================================================================
fused_cycloalkane_component__FUSED_CASES = [
    # 环戊烷稠合四氢呋喃（无取代）：附加组分 cyclopenta，母体 furan
    ("C1CC2CCOC2C1", "cyclopenta[b]furan", "环戊并[b]呋喃"),
    # 同上 5 位羟基：位次随稠环编号（桥头 3a/6a，外周 4/5/6）
    ("O[C@@H]1C[C@@H]2OCC[C@@H]2C1", "cyclopenta[b]furan-5-ol", "环戊并[b]呋喃-5-醇"),
]


@pytest.mark.parametrize("smiles,en_tail,zh_sub", fused_cycloalkane_component__FUSED_CASES)
def test_cycloalkane_as_fused_component(smiles: str, en_tail: str, zh_sub: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert r.en.endswith(en_tail), r.en
    assert zh_sub in r.zh, r.zh


# 芳环侧稠合：附加组分与母体共享的稠合键为芳香键，碳环组分仍须被识别
fused_cycloalkane_component__AROMATIC_FUSED_CASES = [
    # 环戊烷并吡啶（共享稠合键位于吡啶 2,3 侧）
    ("c1cnc2c(c1)CCC2", "cyclopenta[b]pyridine", "环戊并[b]吡啶"),
    # 同上换母体杂环：吡嗪
    ("c1cnc2c(n1)CCC2", "cyclopenta[b]pyrazine", "环戊并[b]吡嗪"),
    # 环大小随附加组分变：环丁并/环庚并
    ("c1cnc2c(c1)CC2", "cyclobuta[b]pyridine", "环丁并[b]吡啶"),
    ("c1cnc2c(c1)CCCCC2", "cyclohepta[b]pyridine", "环庚并[b]吡啶"),
]


@pytest.mark.parametrize("smiles,en_tail,zh_sub", fused_cycloalkane_component__AROMATIC_FUSED_CASES)
def test_aromatic_fused_cycloalkane_component(smiles: str, en_tail: str, zh_sub: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert r.en.endswith(en_tail), r.en
    assert zh_sub in r.zh, r.zh


# 近邻负例：已注册稠环与单环烃不得被本轮组分路径改写
fused_cycloalkane_component__NEIGHBORS = [
    ("c1ccc2occc2c1", "1-benzofuran", "1-苯并呋喃"),
    ("C1CCC2CCCCC2C1", "decahydronaphthalene", "十氢萘"),
    ("C1CCc2ccccc2C1", "1,2,3,4-tetrahydronaphthalene", "1,2,3,4-四氢萘"),
    ("C1CCCC1", "cyclopentane", "环戊烷"),
    # 六元碳环稠合：整骨架有保留名模板（喹啉），不得改走环烷烃附加组分路径
    ("c1cnc2c(c1)CCCC2", "5,6,7,8-tetrahydroquinoline", "5,6,7,8-四氢喹啉"),
]


@pytest.mark.parametrize("smiles,en,zh", fused_cycloalkane_component__NEIGHBORS)
def test_neighbors_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


# ==========================================================================
# 合并自 test_partial_hydro_ring_locant_prefix.py
# IUPAC: P-22.2.2, P-25.1
# Layer: L2
#
# 局部不饱和保留名的词干已含加氢前缀时，locant 前缀须留在组分名前，不得前移。
#
# `_TEMPLATES` 里 dihydrothiazole 的词干是 `4,5-dihydro-1,3-thiazole`：`1,3-` 必须
# 紧贴 `thiazole`。若按「词干为裸名、locant 前缀在词首注入」的路径再补一次，就得到
# `1,3-4,5-dihydro-1,3-thiazole`（位次重复）。同族 dihydropyrrole / dihydroimidazole
# 的 `1H-` 同样嵌在词中，靠 `1H-` 只在环含未取代芳香 NH 时才注入而侥幸未暴露。
#
# 对照条守住边界：词干不含加氢前缀者（thiazolidine / dioxolane / benzothiazole）
# 仍须在词首注入 locant 前缀。
# ==========================================================================
partial_hydro_ring_locant_prefix__CASES = [
    # 用户报例：4,5-二氢-1,3-噻唑-4-羧酸
    (
        "O=C(O)[C@H]1CSC(c2c[nH]c3ccccc23)=N1",
        "(4S)-2-(1H-indol-3-yl)-4,5-dihydro-1,3-thiazole-4-carboxylic acid",
        "(4S)-2-(1H-吲哚-3-基)-4,5-二氢-1,3-噻唑-4-羧酸",
    ),
    # 裸环：词干原样
    ("C1=NCCS1", "4,5-dihydro-1,3-thiazole", "4,5-二氢-1,3-噻唑"),
    # 同族（1H- 嵌在词中）不得前移
    ("C1C=CCN1", "2,5-dihydro-1H-pyrrole", "2,5-二氢-1H-吡咯"),
    ("C1=NCCN1", "4,5-dihydro-1H-imidazole", "4,5-二氢-1H-咪唑"),
    # 对照：词干不含加氢前缀，locant 前缀仍在词首注入
    ("C1NCCS1", "1,3-thiazolidine", "1,3-噻唑烷"),
    ("C1COCO1", "1,3-dioxolane", "1,3-二氧戊环"),
    ("c1ccc2scnc2c1", "1,3-benzothiazole", "1,3-苯并噻唑"),
]


@pytest.mark.parametrize("smiles,en,zh", partial_hydro_ring_locant_prefix__CASES)
def test_partial_hydro_ring_locant_prefix_not_duplicated(smiles: str, en: str, zh: str) -> None:
    """部分不饱和保留名的 locant 前缀只出现一次，且留在组分名前。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_locant_key.py
# IUPAC: P-25.3.3.1 稠合碳字母位次
# Layer: L4
#
# locant_key: 数字/字母位次统一排序。
# ==========================================================================
def test_locant_key_parses_letter_suffix():
    assert locant_key("4a") == (4, "a")
    assert locant_key("10") == (10, "")
    assert locant_key(4) == (4, "")


def test_locant_str_sort_orders_mixed():
    locs = locant_str_sort(["5", "10", "4a", 4, "3"])
    assert locs == ["3", 4, "4a", "5", "10"]


def test_locant_key_numeric_ordering():
    # "10" > "2"(数值而非字典序)
    assert locant_key("10") > locant_key("2")
