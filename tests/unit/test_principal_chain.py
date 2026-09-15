# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_principal_hydrocarbon.py: 
test_principal_chain_fg.py: 
test_p44_2_ring_seniority.py: 
test_p44_4_unsaturation.py: 
test_p44_principal_group_count.py: Principal-group count precedes later parent-selection criteria.
test_p45_candidate_ladder.py: 并列母体候选按 P-45.2.2 前缀取代基位次集合收窄（P-44.1.1 未决时才生效）。
test_topology_principal_gaps.py: 
test_fg_priority.py: 组合官能团主基团优先级仲裁：更高优先级 FG（如自由基）存在时，
"""
from __future__ import annotations

import hashlib
import time

import namepredict.namer as namer_module
import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2 import kind_registry
from namepredict.layer2.parent_select import _collect_candidates
from namepredict.layer2.kind_registry import pack_parent_stem
from namepredict.layer2.parent_select import select_parent
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, keep_p44_2, keep_p44_4_unsaturation
from namepredict.layer2.principal_expression import (
    PrincipalChargeState,
    PrincipalExpressionFacts,
    PrincipalRelation,
    express_chain_principal,
)
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_select import rule_driven_parent_candidates, select_principal_parent_skeletons
from namepredict.layer4.candidate_keys import prefix_locant_set, suffix_locant_set
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.types import NameResult
from rdkit import Chem

# ==========================================================================
# 合并自 test_principal_hydrocarbon.py
# IUPAC: P-44.1.2 / P-44.3 / P-44.4
# Layer: L2
# 
# 验证 rule_driven_parent_candidates 在【无主官能团】时独立产出纯烃候选，
# 不依赖经典生产者通道（ring/unsat/benzene/alkane fallback）。
# 通过 _principal_name 直接走 principal 管线 + L3-L5 装配，避免经典通道假绿。
# ==========================================================================
principal_hydrocarbon__CASES = [
    ("CC", "ethane", "乙烷"),
    ("CCCCC", "pentane", "戊烷"),
    ("C=C", "ethene", "乙烯"),
    ("C#CC", "propyne", "丙炔"),
    ("C=CC=C", "buta-1,3-diene", "丁-1,3-二烯"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("C1CC=CCC1", "cyclohexene", "环己烯"),
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("c1ccccc1CC", "ethylbenzene", "乙基苯"),
]

# 负例：含主官能团，必须继续走主官能团管线，不被纯烃逻辑误伤
principal_hydrocarbon__NEGATIVE = [
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
]


def principal_hydrocarbon___principal_name(smiles: str):
    """Only the rule-driven principal pipeline + L3-L5, no legacy producers."""
    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    t0 = time.perf_counter()
    for parent in rule_driven_parent_candidates(info):
        # 模拟完整管线的词干填充（iter_parent_candidates 内部会 pack）。
        packed = pack_parent_stem(parent, info["mol"])
        parent, subst, complete = namer_module._prepare_candidate(info, packed)
        if not complete:
            continue
        hit = namer_module._assemble_candidate(parent, subst, t0=t0)
        if hit is not None and hit.success:
            return hit
    return None


def test_principal_hydrocarbon_negative_fg_untouched():
    for smiles, en, zh in principal_hydrocarbon__NEGATIVE:
        r = principal_hydrocarbon___principal_name(smiles)
        assert r is not None, f"principal 管线未产出: {smiles}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


# ==========================================================================
# 合并自 test_principal_chain_fg.py
# IUPAC: P-44 / P-65.6 / P-66
# Layer: L2,L3,L4,L5
# 
# 验证 rule-driven principal 管线对无环单官能团 ester/amide/aldehyde/nitrile
# 也产出信息完整的 chain parent dict（kind + 单数 *_c_idx + 不饱和 + 酯烷氧基）。
# 经典 builder 路径 _open_chain_expression 已删除，骨架表达为唯一来源。
# ==========================================================================
principal_chain_fg__CASES = [
    # esters
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("CCCC(=O)OC", "methyl butanoate", "丁酸甲酯"),
    ("C=CC(=O)OC", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
    ("CC=CC(=O)OC", "methyl but-2-enoate", "丁-2-烯酸甲酯"),
    ("C#CCC(=O)OC", "methyl but-3-ynoate", "丁-3-炔酸甲酯"),
    # amides
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCCC(=O)N", "butanamide", "丁酰胺"),
    ("C=CC(=O)N", "acrylamide", "丙烯酰胺"),
    ("CC=CC(=O)N", "but-2-enamide", "丁-2-烯酰胺"),
    # aldehydes
    ("CCCCC=O", "pentanal", "戊醛"),
    ("CC(C)C=O", "2-methylpropanal", "2-甲基丙醛"),
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("C=CCC=O", "but-3-enal", "丁-3-烯醛"),
    # nitriles
    ("CCCC#N", "butanenitrile", "丁腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
    ("C=CC#N", "prop-2-enenitrile", "丙-2-烯腈"),
    ("CC=CC#N", "but-2-enenitrile", "丁-2-烯腈"),
]

# 负例：有环 / 二酯 / 高优先级官能团，必须走既有路径不被误伤
principal_chain_fg__NEGATIVE = [
    ("CC(=O)O", "acetic acid", "乙酸"),                        # 羧酸 > 酯
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),       # 苯环酯 → benzoate
    ("COC(=O)CCC(=O)OC", "dimethyl butanedioate", "丁二酸二甲酯"),  # 二酯
    ("O=C(O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),  # 环酸
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),                 # 环腈
    ("C1CCCCC1", "cyclohexane", "环己烷"),                     # 纯烃环
]


# ── 表达层字段完整性：express_chain_principal 直接产出完整 parent dict ──

def principal_chain_fg___chain_parent(smiles: str, kind: str):
    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    selection = select_principal_parent_skeletons(info)
    if selection.principal is None or selection.skeletons is None:
        return None
    for skeleton in selection.skeletons.candidates:
        if skeleton.topology is not SkeletonTopology.ACYCLIC:
            continue
        parent = express_chain_principal(info, selection.principal, skeleton)
        if parent is not None and parent["kind"] == kind:
            return parent
    return None


def test_chain_amide_aldehyde_nitrile_typed_anchors():
    """链式酰胺/醛/腈的 principal 锚点由 typed facts 承载（单一格式，不再平行写扁平 *_c_idx）。"""
    for smiles, kind in [
        ("CC=CC(=O)N", "amide"),
        ("CC=CC=O", "aldehyde"),
        ("CC=CC#N", "nitrile"),
    ]:
        p = principal_chain_fg___chain_parent(smiles, kind)
        assert p is not None, f"{smiles} 无 {kind} 表达候选"
        facts = p["principal_expression_facts"]
        assert facts.group_class.value == kind
        assert facts.anchor_atoms, f"{smiles} 缺 principal 锚点"
        assert p["n_carbons"] == 4
        assert isinstance(p["double_bond"], tuple), f"{smiles} 缺 double_bond"


def test_chain_nitrile_yne_field():
    p = principal_chain_fg___chain_parent("CC#CC#N", "nitrile")
    assert p is not None
    assert isinstance(p["triple_bond"], tuple) and len(p["triple_bond"]) == 2


# ==========================================================================
# 合并自 test_p44_2_ring_seniority.py
# IUPAC: P-44.2
# Layer: L2
# ==========================================================================
def p44_2_ring_seniority___ring(atoms):
    return ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(atoms), frozenset({"fg:0"}))


def test_p44_2_prefers_heterocycle_over_carbocycle():
    mol = Chem.MolFromSmiles("C1CC1.N1CC1")
    assert keep_p44_2(mol, (p44_2_ring_seniority___ring((0, 1, 2)), p44_2_ring_seniority___ring((3, 4, 5)))) == (p44_2_ring_seniority___ring((3, 4, 5)),)


def test_p44_2_prefers_more_nitrogens():
    mol = Chem.MolFromSmiles("N1CCCC1.N1CCNC1")
    one = p44_2_ring_seniority___ring(range(5))
    two = p44_2_ring_seniority___ring(range(5, 10))
    assert keep_p44_2(mol, (one, two)) == (two,)


def test_p44_2_without_n_prefers_senior_heteroatom():
    mol = Chem.MolFromSmiles("O1CCC1.S1CCC1")
    oxygen = p44_2_ring_seniority___ring(range(4))
    sulfur = p44_2_ring_seniority___ring(range(4, 8))
    assert keep_p44_2(mol, (sulfur, oxygen)) == (oxygen,)


def test_p44_2_prefers_more_rings_before_more_atoms():
    mol = Chem.MolFromSmiles("N1CCC2CC12.N1CCCCC1")
    bicyclic = p44_2_ring_seniority___ring(range(6))
    monocyclic = p44_2_ring_seniority___ring(range(6, 12))
    assert keep_p44_2(mol, (monocyclic, bicyclic)) == (bicyclic,)


def test_p44_2_prefers_more_ring_atoms_after_hetero_ties():
    mol = Chem.MolFromSmiles("N1CCC1.N1CCCC1")
    small = p44_2_ring_seniority___ring(range(4))
    large = p44_2_ring_seniority___ring(range(4, 9))
    assert keep_p44_2(mol, (small, large)) == (large,)


# ==========================================================================
# 合并自 test_p44_4_unsaturation.py
# IUPAC: P-44.4.1.1 / P-44.4.1.2
# Layer: L2
# ==========================================================================
def p44_4_unsaturation___chain(atoms):
    return ParentSkeleton(SkeletonTopology.ACYCLIC, tuple(atoms), frozenset({"fg:0"}))


def test_p44_4_prefers_more_multiple_bonds():
    mol = Chem.MolFromSmiles("C=CC.C=CC#C")
    one, two = p44_4_unsaturation___chain(range(3)), p44_4_unsaturation___chain(range(3, 7))
    assert keep_p44_4_unsaturation(mol, (one, two)) == (two,)


def test_p44_4_prefers_more_double_bonds_when_total_ties():
    mol = Chem.MolFromSmiles("C=CC=C.C=CC#C")
    diene, enyne = p44_4_unsaturation___chain(range(4)), p44_4_unsaturation___chain(range(4, 8))
    assert keep_p44_4_unsaturation(mol, (enyne, diene)) == (diene,)


def test_p44_4_aromatic_counts_as_kekule_multiple_bonds():
    mol = Chem.MolFromSmiles("c1ccccc1.C1CCCCC1")
    aromatic, saturated = p44_4_unsaturation___chain(range(6)), p44_4_unsaturation___chain(range(6, 12))
    assert keep_p44_4_unsaturation(mol, (aromatic, saturated)) == (aromatic,)


def test_p44_4_pyrimidine_beats_saturated_ring():
    # 哌嗪-嘧啶环环相连：嘧啶(芳香,3 多重键当量)应优先于哌嗪(饱和,0)。
    mol = Chem.MolFromSmiles("C1C(N2C(C)CNCC2)=NC=NC=1")
    pymidine = ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple([0, 1, 9, 10, 11, 12]), frozenset())
    piperazine = ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple([2, 3, 5, 6, 7, 8]), frozenset())
    assert keep_p44_4_unsaturation(mol, (pymidine, piperazine)) == (pymidine,)


# ==========================================================================
# 合并自 test_p44_principal_group_count.py
# IUPAC: P-44.1.1
# Layer: L2
#
# Principal-group count precedes later parent-selection criteria.
# ==========================================================================
def p44_principal_group_count___result(en):
    return NameResult(en=en, zh=en, success=True, source="iupac", time_ms=0.0)


def p44_principal_group_count___prepared(kind, complete):
    return ({"kind": kind}, [], complete)


def test_higher_phase_partial_preserves_seniority(monkeypatch):
    phases = [[{"kind": "acid"}], [{"kind": "alcohol"}]]
    monkeypatch.setattr(namer_module, "_candidate_phases", lambda info: phases)
    monkeypatch.setattr(
        namer_module, "_prepare_candidate",
        lambda info, parent, **kwargs: p44_principal_group_count___prepared(parent["kind"], parent["kind"] == "alcohol"),
    )
    monkeypatch.setattr(namer_module, "_assemble_candidate", lambda p, *a, **k: p44_principal_group_count___result(p["kind"]))
    assert namer_module._run_candidates({}, t0=0).en == "acid"


def test_higher_phase_complete_never_downgrades(monkeypatch):
    phases = [[{"kind": "acid"}], [{"kind": "alcohol"}]]
    monkeypatch.setattr(namer_module, "_candidate_phases", lambda info: phases)
    monkeypatch.setattr(
        namer_module, "_prepare_candidate",
        lambda info, parent, **kwargs: p44_principal_group_count___prepared(parent["kind"], True),
    )
    monkeypatch.setattr(namer_module, "_assemble_candidate", lambda p, *a, **k: p44_principal_group_count___result(p["kind"]))
    assert namer_module._run_candidates({}, t0=0).en == "acid"

p44_principal_group_count__CASES = [
    # Near-neighbour negative: the already claimable C1 arm stays unchanged.
    ("OCC(O)CO", "propane-1,2,3-triol", "丙烷-1,2,3-三醇"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
]


@pytest.mark.parametrize("smiles,en,zh", p44_principal_group_count__CASES)
def test_principal_group_count_controls_entry_parent(smiles, en, zh):
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(result.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_p45_candidate_ladder.py
# IUPAC: P-45.2.2
# Layer: L2,L4
#
# 并列母体候选按 P-45.2.2 前缀取代基位次集合收窄（P-44.1.1 未决时才生效）。
# ==========================================================================
p45_candidate_ladder__PIN_CASES = [
    # 文档 P-45.2.2 邻近案例：位次集合 '2,4' 低于 '2,5'。
    ("Oc1c(CCc2cc(O)c(Cl)cc2)cc(Cl)cc1",
     "4-chloro-2-[2-(4-chloro-3-hydroxyphenyl)ethyl]phenol"),
]

# 近邻负例：后缀位次集合不同（'2,4,5' vs '3,4,5'），P-44.1.1 先决，P-45.2.2 不得越级。
p45_candidate_ladder__GUARD_CASES = [
    ("OC[C@H]1OC(O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O",
     "(3R,4S,5R,6R)-6-(hydroxymethyl)-3-[(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxane-2,4,5-triol"),
]


def test_prefix_locant_set_is_ascending_with_duplicates():
    """位次集合升序排列且保留重复（P-14.3.5）。"""
    numbered = {"substituents": [{"locant": 4}, {"locant": 2}, {"locant": 2}]}
    assert prefix_locant_set(numbered) == ((2, ""), (2, ""), (4, ""))


def test_prefix_locant_set_ignores_locantless_prefixes():
    """近邻负例：无位次前缀（N- 取代基等）不入键，缺 substituents 返回空元组。"""
    assert prefix_locant_set({"substituents": [{"locant": None}, {"locant": 3}]}) == ((3, ""),)
    assert prefix_locant_set({}) == ()


def p45_candidate_ladder___facts(attachment: set) -> PrincipalExpressionFacts:
    """构造只填附着原子的 principal facts，供 suffix_locant_set 纯函数测试。"""
    return PrincipalExpressionFacts(
        group_class=FunctionalGroupClass.ALCOHOL, multiplicity=len(attachment),
        relation=PrincipalRelation.IN_SKELETON, occurrence_ids=(),
        characteristic_atoms=frozenset(), anchor_atoms=frozenset(attachment),
        attachment_atoms=frozenset(attachment), charge_state=PrincipalChargeState.NEUTRAL)


def test_suffix_locant_set_uses_principal_attachment_atoms():
    """P-44.1.1 键取 principal 特征基团附着原子的链上位次，升序去重前保留重复。"""
    numbered = {"parent": {"chain": [10, 11, 12], "principal_expression_facts": p45_candidate_ladder___facts({12, 10})}}
    assert suffix_locant_set(numbered) == ((1, ""), (3, ""))
    assert suffix_locant_set({"parent": {"chain": [1, 2]}}) == ()


p45_candidate_ladder__TIED_SMILES = "Oc1c(CCc2cc(O)c(Cl)cc2)cc(Cl)cc1"


def test_select_parent_keeps_tied_group_finalized():
    """存在并列候选时返回整组而非单个，组内候选均已终态化（owned_atoms 固化为 frozenset）。"""
    info = analyze(preprocess(p45_candidate_ladder__TIED_SMILES))
    group = select_parent(info)
    assert len(group) > 1
    assert all(isinstance(c.get("owned_atoms"), frozenset) for c in group)
    # 并列候选是分子上不同的两个酚环，不是同一母体的重复项
    assert len({tuple(sorted(c["owned_atoms"])) for c in group}) == len(group)


@pytest.mark.parametrize("smiles,en", p45_candidate_ladder__PIN_CASES)
def test_tied_candidates_resolved_by_p45_2_2(smiles, en):
    """P-45.2.2 决定并列候选：前缀位次集合更小者胜出。"""
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,en", p45_candidate_ladder__GUARD_CASES)
def test_suffix_locant_difference_blocks_p45_2_2(smiles, en):
    """近邻负例：后缀位次集合不同时 P-44.1.1 先决，P-45.2.2 不介入，候选顺序不变。"""
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


def test_p45_2_2_tie_keeps_previous_candidate():
    """近邻负例：两侧位次集合全平（无位次前缀）时 P-45.2.2 不介入，候选顺序不变。"""
    result = SMILESNNamer().name("O=C(COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H](O)COP(=O)(O)O")
    assert result.success
    assert result.en == "(3S,4R,5R)-3,4,5-trihydroxy-2-oxo-6-phosphonooxyhexyl dihydrogen phosphate"


def test_single_candidate_unchanged():
    """近邻负例：无并列候选的普通分子行为不变。"""
    result = SMILESNNamer().name("OCCO")
    assert result.success
    assert normalize_en(result.en) == normalize_en("ethane-1,2-diol")


# ==========================================================================
# 合并自 test_topology_principal_gaps.py
# IUPAC: P-44
# Layer: L1,L2
# ==========================================================================
topology_principal_gaps__CASES = [
    # 苯环 + FG 的 kind 收敛为 FG 类别；保留名由 L5 typed_kinds 决定。
    ("O=C(O)c1ccccc1", "acid", "benzoic acid"),
    ("O=Cc1ccccc1", "aldehyde", "benzaldehyde"),
    ("N#Cc1ccccc1", "nitrile", "benzonitrile"),
    ("NC(=O)c1ccccc1", "amide", "benzamide"),
    ("Oc1ccccc1", "alcohol", "phenol"),
    ("Nc1ccccc1", "amine", "aniline"),
    ("CCCC(=O)OC", "ester", "methyl butanoate"),
    ("CCCC(=O)N", "amide", "butanamide"),
    ("CCCC=O", "aldehyde", "butanal"),
    ("CC(=O)CC", "ketone", "butan-2-one"),
    ("CCCC#N", "nitrile", "butanenitrile"),
]


@pytest.mark.parametrize("smiles,kind,en", topology_principal_gaps__CASES)
def test_topology_first_principal_path_owns_target(smiles, kind, en):
    # fg_try_fns 已随 fg 注册层删除；principal-only 是默认状态。
    info = analyze(Chem.MolFromSmiles(smiles))
    parents = _collect_candidates(info)
    result = SMILESNNamer().name(smiles)
    owned = [parent for parent in parents if parent.get("kind") == kind]
    assert owned
    count = 2 if kind == "diester" else 1
    assert owned[0].get("principal_group_count") == count
    assert len(owned[0].get("covered_principal_ids") or ()) == count
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)


def test_open_chain_monoester_is_not_claimed_as_diester():
    info = analyze(Chem.MolFromSmiles("COC(=O)CCC"))
    assert all(parent.get("kind") != "diester" for parent in _collect_candidates(info))


@pytest.mark.parametrize("smiles,kind,group_class", [
    # 环 + FG 一律收敛为 FG 类别 kind（苯/饱和环/稠环/杂环平等，词干由 scaffold 承载）。
    ("O=C(O)C1CCCCC1", "acid", "acid"),
    ("N#CC1CCCCC1", "nitrile", "nitrile"),
    ("O=C(O)c1ccc2ccccc2c1", "acid", "acid"),
])
def test_ring_principal_uses_base_scaffold_kind(smiles, kind, group_class):
    parents = _collect_candidates(analyze(Chem.MolFromSmiles(smiles)))
    matched = [parent for parent in parents if parent.get("kind") == kind]
    assert matched
    assert matched[0]["principal_expression_facts"].group_class.value == group_class
    assert "carboxylic" not in matched[0]["kind"]
    assert "carbonitrile" not in matched[0]["kind"]


# ==========================================================================
# 合并自 test_fg_priority.py
# IUPAC: P-41
# Layer: L1
#
# 组合官能团主基团优先级仲裁：更高优先级 FG（如自由基）存在时，
# 组合羰基 FG（酰胺/醛/羧酸/酯）退出主基团并降级为 oxo 前缀，
# 不再丢失羰基氧。
#
# P-41 规定 radical 的类优先级高于酸/酯/酰胺/醛；非最高组合 FG 以
# 前缀表达（oxo/amino/…），而不是被丢弃。
# ==========================================================================
fg_priority__POSITIVE = [
    ("C(*)C(N)=O", "2-amino-2-oxoethyl", "2-氨基-2-氧代乙基"),
    ("C(*)C=O", "2-oxoethyl", "2-氧代乙基"),
    ("*CCC(=O)O", "2-carboxyethyl", "2-羧基乙基"),
    ("C(*)C(=O)OC", "2-methoxy-2-oxoethyl", "2-甲氧基-2-氧代乙基"),
]


@pytest.mark.parametrize("smiles,en,zh", fg_priority__POSITIVE)
def test_composite_fg_degrade_to_oxo(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 负例：组合 FG 自身是最高优先级时，保持主基团命名不变。
fg_priority__NEGATIVE = [
    ("CC(N)=O", "acetamide", "乙酰胺"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCC(O)=O", "propanoic acid", "丙酸"),
    ("CCC(N)=O", "propanamide", "丙酰胺"),
    ("C(*)CO", "2-hydroxyethyl", "2-羟基乙基"),
    ("C(*)C(C)=O", "2-oxopropyl", "2-氧代丙基"),
]


@pytest.mark.parametrize("smiles,en,zh", fg_priority__NEGATIVE)
def test_composite_fg_principal_unchanged(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
