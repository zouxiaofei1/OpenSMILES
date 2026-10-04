# 合并自 3 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_ring_systems.py: Ring-system topology: fusion components for naphthalene/indole/quinoline/spiro.
test_ring_template_match.py: SMILES 模板子图同构识别（原 scaffold/retained_templates.py 移植）。
test_ring_geometry.py: ring_geometry: 3-8 元环模板坐标与平面几何原语。
"""
from __future__ import annotations

import math

from opensmiles.layer0.preprocessor import preprocess
from opensmiles.layer1.analyzer import analyze
from opensmiles.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology
from opensmiles.layer2.ring_scaffold import match_retained, resolve_ring_scaffold
from opensmiles.layer4.fused_numbering import _boundary_walk
from opensmiles.layer4.ring_geometry import RING_TEMPLATES, clip_polygon, overlap_area, polygon_area, regular_polygon, ring_cyclic, ring_shape_template
from opensmiles.namer import SMILESNNamer

# ==========================================================================
# 合并自 test_ring_systems.py
# IUPAC: P-25 / P-22
# Layer: L1
#
# Ring-system topology: fusion components for naphthalene/indole/quinoline/spiro.
# ==========================================================================
def ring_systems___systems(smiles: str) -> list[dict]:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)["ring_systems"]


def test_benzene_mono():
    syss = ring_systems___systems("c1ccccc1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 1
    assert s["n_atoms"] == 6
    assert s["hetero_atoms"] == []


def test_naphthalene_fused2():
    syss = ring_systems___systems("c1ccc2ccccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["n_atoms"] == 10
    assert len(s["fusion_edges"]) == 1
    assert len(s["fusion_edges"][0][2]) == 2  # shared atoms


def test_indole_fused56():
    syss = ring_systems___systems("c1ccc2[nH]ccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["n_atoms"] == 9
    zs = {h["Z"] for h in s["hetero_atoms"]}
    assert 7 in zs


def test_quinoline_fused66():
    syss = ring_systems___systems("c1ccc2ncccc2c1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 2
    assert s["n_atoms"] == 10
    assert any(h["Z"] == 7 for h in s["hetero_atoms"])


def test_biphenyl_two_systems():
    """Two unfused rings → two mono systems (not fused)."""
    syss = ring_systems___systems("c1ccc(-c2ccccc2)cc1")
    assert len(syss) == 2
    assert all(s["n_rings"] == 1 for s in syss)


def test_non_spiro_unchanged():
    """Fused rings should NOT become spiro."""
    syss = ring_systems___systems("c1ccc2ccccc2c1")  # naphthalene
    assert len(syss) == 1
    assert syss[0]["topology"] != "spiro"


def test_biphenyl_unchanged():
    """Two separate rings should remain separate."""
    syss = ring_systems___systems("c1ccc(-c2ccccc2)cc1")
    assert len(syss) == 2
    assert all(s["topology"] != "spiro" for s in syss)


def test_pyridine_hetero():
    syss = ring_systems___systems("c1ccncc1")
    assert len(syss) == 1
    s = syss[0]
    assert s["n_rings"] == 1
    assert any(h["Z"] == 7 for h in s["hetero_atoms"])


def test_open_chain_empty():
    assert ring_systems___systems("CCCC") == []


# ==========================================================================
# 合并自 test_ring_template_match.py
# IUPAC: P-22 / P-25 template-based retained ring identification
# Layer: L2
#
# SMILES 模板子图同构识别（原 scaffold/retained_templates.py 移植）。
#
# 五元组 _TOPOLOGY 无法区分位置异构体（azulene/naphthalene、
# isoindole/indole 等字段全等），模板方案修复这些误配。
# ==========================================================================
def ring_template_match___resolve(smiles: str) -> str | None:
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    system = info["ring_systems"][0]
    skeleton = ParentSkeleton(
        SkeletonTopology.RING_SYSTEM, tuple(system["atom_ids"]), frozenset(),
    )
    identity = resolve_ring_scaffold(info, skeleton)
    return identity.id if identity else None


def ring_template_match___template_id(smiles: str) -> str | None:
    mol = preprocess(smiles)
    assert mol is not None
    info = analyze(mol)
    return match_retained(info, info["ring_systems"][0]["atom_ids"])


def test_registered_templates_resolve():
    assert ring_template_match___resolve("c1ccccc1") == "benzene"
    assert ring_template_match___resolve("c1ccncc1") == "pyridine"
    assert ring_template_match___resolve("c1ccc2ccccc2c1") == "naphthalene"
    assert ring_template_match___resolve("c1ccc2[nH]ccc2c1") == "indole"


def test_azulene_no_longer_misidentified_as_naphthalene():
    # 5+7 稠合全碳环：五元组与 naphthalene 字段全等曾误配为 naphthalene。
    assert ring_template_match___resolve("C1=CC2=CC=CC=CC2=CC1") == "carbocycle"


def test_twistane_no_longer_misidentified_as_naphthalene():
    # 扭曲烷三环骨架：naphthalene 是其非诱导子图（多余环闭合键不可见），曾误配为 naphthalene。
    assert ring_template_match___template_id("C12C3CC(C(C1)CC3)CC2") is None
    assert ring_template_match___resolve("C12C3CC(C(C1)CC3)CC2") == "fused_hetero"


def test_isoindole_no_longer_misidentified_as_indole():
    # 苯并[c]吡咯：五元组与 indole 字段全等曾误配为 indole；未注册非全碳
    # 多环回落 fused_hetero（L5 按 fused_tree 组装稠合名），仍非 indole。
    assert ring_template_match___resolve("c1ccc2c(c1)c[nH]c2") == "fused_hetero"


def test_substituted_ring_still_resolves():
    assert ring_template_match___resolve("Oc1ccccc1") == "benzene"
    assert ring_template_match___resolve("Oc1cccc2ccccc12") == "naphthalene"


def test_all_templates_resolve_to_own_scaffold():
    # 唯一事实来源：模板全部派生 spec，命中即解析为自身。
    assert ring_template_match___resolve("c1ccoc1") == "furan"
    assert ring_template_match___resolve("c1ccc2cc3ccccc3cc2c1") == "anthracene"


def test_hydrogenated_fused_rings_hit_unsaturated_parent():
    # 骨架完全氢化后比较：加氢/部分加氢的稠环仍归属其未饱和保留母体
    # （P-25.3.4 加氢衍生物：四氢萘/十氢萘/二氢吲哚均由保留母体命名）。
    assert ring_template_match___resolve("C1CCc2ccccc2C1") == "naphthalene"   # 1,2,3,4-四氢萘
    assert ring_template_match___resolve("C1=CCc2ccccc2C1") == "naphthalene"  # 1,2-二氢萘
    assert ring_template_match___resolve("C1CCC2CCCCC2C1") == "naphthalene"   # 十氢萘
    assert ring_template_match___resolve("C1Cc2ccccc2N1") == "indole"         # 二氢吲哚


def test_saturated_skeleton_near_neighbours_stay_unmatched():
    # 近邻负例：环数/碳数相近但稠合拓扑不同，仍不得误配为保留稠环。
    assert ring_template_match___resolve("C1Cc2ccccc2C1") == "indene"            # 茚满 C9：归 2,3-二氢-1H-茚，非萘 C10
    assert ring_template_match___resolve("C1=CC2=CC=CC=CC2=CC1") == "carbocycle" # 薁 5+7，非萘 6+6


def test_positional_isomers_hit_own_template_only():
    # 元素标注的子图同构区分位置异构体，各自单命中。
    cases = {
        "quinoline": "c1ccc2ncccc2c1",
        "isoquinoline": "c1nccc2ccccc21",
        "pyridazine": "c1ccnnc1",
        "pyrimidine": "c1cncnc1",
        "pyrazine": "c1cnccn1",
        "imidazole": "c1cnc[nH]1",
        "pyrazole": "c1ccn[nH]1",
        "benzimidazole": "c1ccc2[nH]cnc2c1",
        "indazole": "c1ccc2cn[nH]c2c1",
        "quinazoline": "c1ccc2ncncc2c1",
        "quinoxaline": "c1ccc2nccnc2c1",
    }
    for sid, smi in cases.items():
        assert ring_template_match___template_id(smi) == sid


# ==========================================================================
# 合并自 test_ring_geometry.py
# IUPAC: P-25.3.2.3.1
# Layer: L4
#
# ring_geometry: 3-8 元环模板坐标与平面几何原语。
# ==========================================================================
def ring_geometry___edge_len(p1, p2):
    return math.hypot(p2[0] - p1[0], p2[1] - p1[1])


def test_regular_polygon_edges_equal_and_vertical():
    for n in range(3, 9):
        pts = RING_TEMPLATES[n]
        assert len(pts) == n
        lens = [ring_geometry___edge_len(pts[i], pts[(i + 1) % n]) for i in range(n)]
        assert max(lens) - min(lens) < 1e-9, n
        # 边(0,1) 竖直（右侧，x 相等）
        assert abs(pts[0][0] - pts[1][0]) < 1e-9
        assert pts[0][0] > 0


def test_regular_polygon_odd_mirror_has_two_variants():
    """奇环模板镜像后仍为正多边形且竖直边在左侧。"""
    for n in (3, 5, 7):
        pts = regular_polygon(n, right_edge_vertical=False)
        assert abs(pts[0][0] - pts[1][0]) < 1e-9
        assert pts[0][0] < 0  # 竖直边在左侧


def test_ring_shape_template_distorted_pentagon():
    """5 元行中双侧融合的变形五元环: 两共享边平行(铺行后竖直)、5 顶点、边01 单位对齐。"""
    order = [0, 1, 2, 3, 4]  # 左共享边(0,1), 右共享边(2,3)
    tpl = ring_shape_template(order, 2)
    assert tpl is not None
    assert len(tpl) == 5
    assert tpl[0] == (0.0, 0.0) and tpl[1] == (1.0, 0.0)  # 边01 为单位对齐边
    # 右共享边(2,3) 与左共享边(0,1) 平行(都沿 x) → 铺行后都竖直
    def _parallel(p0, p1, q0, q1):
        return abs((p1[1] - p0[1]) * (q1[0] - q0[0]) - (p1[0] - p0[0]) * (q1[1] - q0[1])) < 1e-9
    assert _parallel(tpl[0], tpl[1], tpl[2], tpl[3])
    # 非双边(偶环/奇数环端部)返回 None, 调用方保持原弃行逻辑
    assert ring_shape_template(list(range(6)), 2) is None


def test_ring_cyclic_forward_and_reverse():
    ring = (0, 1, 2, 3, 4, 5)  # 六元环
    # b 是 a 的后继 → 正序
    assert ring_cyclic(ring, 0, 1) == [0, 1, 2, 3, 4, 5]
    # b 是 a 的前驱 → 反序
    assert ring_cyclic(ring, 0, 5) == [0, 5, 4, 3, 2, 1]


def test_clip_polygon_disjoint_returns_empty():
    a = [(0, 0), (1, 0), (1, 1), (0, 1)]
    b = [(3, 3), (4, 3), (4, 4), (3, 4)]
    assert clip_polygon(a, b) == []


def test_clip_polygon_intersect_positive_area():
    a = [(0, 0), (2, 0), (2, 2), (0, 2)]
    b = [(1, 1), (3, 1), (3, 3), (1, 3)]
    inter = clip_polygon(a, b)
    assert polygon_area(inter) > 0


def test_polygon_area_rectangle():
    assert abs(polygon_area([(0, 0), (2, 0), (2, 1), (0, 1)]) - 2.0) < 1e-9


def test_overlap_area_known():
    a = [(0, 0), (2, 0), (2, 2), (0, 2)]
    b = [(1, 1), (3, 1), (3, 3), (1, 3)]
    assert abs(overlap_area(a, b) - 1.0) < 1e-9
    # 重叠面积方向对称
    assert abs(overlap_area(b, a) - 1.0) < 1e-9


# ==========================================================================
# 合并自 test_fused_numbering.py
# IUPAC: P-25.3.3.1.2
# Layer: L4
#
# 稠环外边界行走：外边界图须为简单环，否则拒绝而非绕内环死循环。
# ==========================================================================
def test_boundary_walk_accepts_simple_cycle():
    """正方形外边界：走满 4 个顶点并闭回起点。"""
    coords = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (1.0, 1.0), 2: (0.0, 1.0)}
    neighbors = {0: {1, 2}, 1: {0, 3}, 2: {0, 3}, 3: {1, 2}}
    exterior = frozenset(frozenset(e) for e in ((0, 1), (1, 3), (2, 3), (0, 2)))
    walk = _boundary_walk(coords, neighbors, exterior, 0)
    assert len(walk) == 4 and sorted(walk) == [0, 1, 2, 3]


def test_boundary_walk_rejects_fork():
    """中途外部度为 3：拒绝，不绕内环打转。"""
    coords = {a: (float(a), 0.0) for a in range(5)}
    neighbors = {0: {1, 2}, 1: {0, 3, 4}, 2: {0, 3}, 3: {1, 2}, 4: {1}}
    exterior = frozenset(frozenset(e) for e in ((0, 1), (1, 3), (2, 3), (0, 2), (1, 4)))
    assert _boundary_walk(coords, neighbors, exterior, 0) == []


def test_boundary_walk_rejects_forked_start():
    """起点外部度非 2：拒绝。"""
    coords = {a: (float(a), 0.0) for a in range(4)}
    neighbors = {0: {1, 2, 3}, 1: {0, 2}, 2: {0, 1}, 3: {0}}
    exterior = frozenset(frozenset(e) for e in ((0, 1), (1, 2), (0, 2), (0, 3)))
    assert _boundary_walk(coords, neighbors, exterior, 0) == []


def test_cage_carbohydride_naming_terminates():
    """高稠合笼状输入：稠环无候选后回落，命名必须终止。"""
    result = SMILESNNamer().name("C123C45C61C42C1C5C6C31")
    assert result.en and result.zh
