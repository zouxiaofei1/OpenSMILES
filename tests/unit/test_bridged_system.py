# IUPAC: P-23.2.6
# Layer: L2
#
# 桥环（扩展 von Baeyer 系统）Layer2 拆解：主环/主桥/次级桥划分、全原子位次、von Baeyer 描述符。
# 本轮只做 Layer2，因此全部直接断言 L2 节点，不测端到端名字。
from __future__ import annotations

import pytest

from opensmiles.layer0.preprocessor import preprocess
from opensmiles.layer1.analyzer import analyze
from opensmiles.layer2.bridged_system import decompose_bridged_system
from opensmiles.layer2.parent_select import select_parent


def _nodes(smiles: str) -> list:
    """对分子的每个环系跑桥环拆解，收集全部候选节点。"""
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    out: list = []
    for system in info.get("ring_systems") or []:
        out.extend(decompose_bridged_system(info, system))
    return out


# 正例：descriptor 与次级桥上标位次，均由 P-23.2.1~P-23.2.6.4 手工推导，金标仅作旁证。
# 每项 = (SMILES, 期望 descriptor, 期望 locant_pairs, 金标旁证)
POSITIVE = [
    ("C1CC2CCC1C2", (2, 2, 1), (), "bicyclo[2.2.1]heptane"),
    ("C12CNCC(CC1)CC2", (3, 2, 2), (), "3-azabicyclo[3.2.2]nonane"),
    ("CC12CCC(C1)C(C)(C)C2", (2, 2, 1), (), "1,3,3-trimethylbicyclo[2.2.1]heptane"),
    ("O=C1C2C(C2CC1)C(=O)OCC", (3, 1, 0), (), "bicyclo[3.1.0]hexane"),  # 0 原子主桥
    ("C1C2CC3CC1CC(C2)C3", (3, 3, 1, 1), ((3, 7),), "tricyclo[3.3.1.13,7]decane"),  # 金刚烷
    ("O=S1(=O)N2CN3CCN(C2)CN1C3", (4, 3, 1, 1), ((3, 8),),
     "tricyclo[4.3.1.13,8]undecane 9,9-dioxide"),
    # 扭曲烷：两个 0 原子次级桥；萘是它的非诱导子图，勿被稠环路径截走
    ("C12C3CC(C(C1)CC3)CC2", (4, 4, 0, 0), ((3, 8),), "tricyclo[4.4.0.03,8]decane"),
]

# 金标精确回归：钉死窄化次序（P-23.2.4 主桥最长先于 P-23.2.6.2.1 对称）
# 与次级桥引用顺序（P-23.2.6.2.5 位次序列最小化）。
GOLDEN = [
    # tetracyclo[10.2.1.01,9.03,8]pentadeca-…（十氢化物）；对称优先会给 (7,5,…)
    ("C=C1C[C@]23C[C@H]1CC[C@H]2[C@]1(C(=O)[O-])CCC[C@@](C)(C(=O)[O-])[C@H]1[C@@H]3C(=O)[O-]",
     (10, 2, 1, 0, 0), ((1, 9), (3, 8))),
    # 金标 tricyclo[5.3.1.01,5]undec-8-ene：次级桥位次对 (1,5)
    ("CC1=CC[C@@]23C[C@@H]1C(C)(C)[C@@H]2CC[C@H]3C", (5, 3, 1, 0), ((1, 5),)),
    # 两条 0 原子次级桥：引用顺序取位次序列更低的 (2,7) 在前
    ("C1=C2[C@H]3C[C@@H](CN2CCC1)[C@@H]1CCCCN1C3", (7, 7, 1, 0, 0), ((2, 7), (10, 15))),
]


@pytest.mark.parametrize("smiles,descriptor,locants", GOLDEN)
def test_golden_descriptor_and_locants(smiles, descriptor, locants):
    nodes = _nodes(smiles)
    assert nodes
    assert nodes[0].descriptor == descriptor
    assert nodes[0].locant_pairs == locants


# 只断言 P-23.2.6.1.4 硬自检的用例（金标与结构不自洽，见下）
INVARIANT_ONLY = [
    # 该分子环系 10 个原子（V=10, E=13, r=3），而合并金标 tricyclo[3.2.1.02,7]dec-3-en-8-ol
    # 的数字和 +2 = 8，二者矛盾 —— 金标不可信，只断言自检与覆盖完整性。
    "C12C(CC3CCC(C=C13)O)C2",
]


@pytest.mark.parametrize("smiles,descriptor,locants,golden", POSITIVE)
def test_descriptor_and_locants(smiles, descriptor, locants, golden):
    nodes = _nodes(smiles)
    assert nodes, f"{smiles} 未产出桥环节点（金标旁证 {golden}）"
    node = nodes[0]
    assert node.descriptor == descriptor
    assert node.locant_pairs == locants
    assert node.n_rings == len(descriptor) - 1


@pytest.mark.parametrize("smiles", [c[0] for c in POSITIVE] + INVARIANT_ONLY)
def test_descriptor_sums_to_ring_atoms(smiles):
    """P-23.2.6.1.4：环原子总数 = 描述符数字之和 + 2。"""
    nodes = _nodes(smiles)
    assert nodes, f"{smiles} 未产出桥环节点"
    node = nodes[0]
    assert sum(node.descriptor) + 2 == len(node.atom_ids)


@pytest.mark.parametrize("smiles", [c[0] for c in POSITIVE] + [c[0] for c in GOLDEN]
                         + INVARIANT_ONLY)
def test_main_ring_is_real_cycle(smiles):
    """主环路径必须沿真实化学键闭合（换长边时若不掉转行进方向，路径会断裂）。"""
    mol = preprocess(smiles)
    seq = list(_nodes(smiles)[0].main_ring)
    k = len(seq)
    broken = [(seq[i], seq[(i + 1) % k]) for i in range(k)
              if not mol.GetBondBetweenAtoms(seq[i], seq[(i + 1) % k])]
    assert not broken, f"{smiles}: 主环路径断裂于 {broken}"


@pytest.mark.parametrize("smiles", [c[0] for c in POSITIVE] + [c[0] for c in GOLDEN]
                         + INVARIANT_ONLY)
def test_bridgeheads_take_locants_1_and_6(smiles):
    """P-23.2.3：主桥两端桥头位次必须是 1 与长边原子数 + 2。"""
    node = _nodes(smiles)[0]
    assert sorted(node.numbering[h] for h in node.main_bridge.heads) == \
        [1, node.ring_segments[0] + 2]


def test_numbering_covers_ring_system_once():
    """位次 1..N 必须无重复无遗漏地覆盖全部环系原子。"""
    nodes = _nodes("C1C2CC3CC1CC(C2)C3")  # 金刚烷，含 1 条独立次级桥
    node = nodes[0]
    assert sorted(node.numbering) == sorted(node.atom_ids)
    assert sorted(node.numbering.values()) == list(range(1, len(node.atom_ids) + 1))


def test_main_ring_and_bridge_partition():
    """主环 + 主桥 + 次级桥必须恰好划分全部环系原子。"""
    nodes = _nodes("O=S1(=O)N2CN3CCN(C2)CN1C3")
    node = nodes[0]
    covered = set(node.main_ring) | set(node.main_bridge.atoms)
    covered |= {a for seg in node.secondary_bridges for a in seg.atoms}
    covered |= set(node.main_bridge.heads)
    covered |= {h for seg in node.secondary_bridges for h in seg.heads}
    assert covered == set(node.atom_ids)


def test_hetero_bridgehead_numbering():
    """P-23.2.3：3-氮杂双环[3.2.2]壬烷的环氮应得位次 3。"""
    mol = preprocess("C12CNCC(CC1)CC2")
    n_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    nodes = _nodes("C12CNCC(CC1)CC2")
    assert nodes[0].numbering[n_idx] == 3


def test_hetero_classification():
    """P-23.3.1：含杂原子的骨架 ring 应为 hetero；纯碳为 carbo。"""
    assert _nodes("C12CNCC(CC1)CC2")[0].ring == "hetero"
    assert _nodes("C1CC2CCC1C2")[0].ring == "carbo"


# 负例：拆不出桥环，必须不产出节点
# 注意萘/十氢萘不在此列 —— 其骨架是合法的 von Baeyer 双环[4.4.0]癸烷（4a-8a 为 0 原子主桥），
# 排除它们靠的是入口门（scaffold 非 fused_hetero 或 P25 能产名），见下面的门控用例。
NEGATIVE = [
    ("C1CCC2(CC1)CCCCC2", "螺环：去桥头分量端点数不为 2"),  # 螺环
    ("c1ccccc1", "苯：环数 < 2"),                          # 单环
    ("C1CCCCC1", "环己烷：环数 < 2"),                       # 单环
]


@pytest.mark.parametrize("smiles,why", NEGATIVE)
def test_not_bridged(smiles, why):
    assert _nodes(smiles) == [], f"{smiles} 不应被当作桥环（{why}）"


def test_fused_system_not_hijacked_into_bridged():
    """入口三重门：P25 能产名的稠环/单环不得被标成 bridged。"""
    for smiles in ("c1ccc2ccccc2c1", "C1CCC2CCCCC2C1", "c1ccc2[nH]ccc2c1", "c1ccccc1"):
        mol = preprocess(smiles)
        parents = select_parent(analyze(mol))
        assert all(p.get("scaffold_id") != "bridged" for p in parents), smiles


@pytest.mark.parametrize("smiles,descriptor", [(c[0], c[1]) for c in POSITIVE])
def test_entry_gate_marks_bridged(smiles, descriptor):
    """非 mancude 且桥不全为 0 的环系直接路由到桥环（P-23）。"""
    mol = preprocess(smiles)
    parents = [p for p in select_parent(analyze(mol)) if p.get("scaffold_id") == "bridged"]
    assert len(parents) == 1
    p = parents[0]
    # kind 可能被 FG 类正交化（如环酯走 "ester"），桥环身份一律由 scaffold_id 承载
    assert p["scaffold_id"] == "bridged"
    assert p["scaffold_identity"].naming_class == "bridged"
    assert p.get("fused_tree") is None  # 桥环身份由 bridged_node 承载，不写 fused_tree
    assert p["bridged_node"].descriptor == descriptor
