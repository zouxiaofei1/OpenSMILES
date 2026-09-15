"""为选中的骨架标注主基团表达式 facts。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum

from namepredict.constants import (
    HALO_Z, MONONUCLEAR_BY_ELEMENT, MONONUCLEAR_HYDRIDES, NITROGEN_STEM_BY_FREE_DOUBLE,
    O, PHOSPHORUS_STEM_BY_OXO, SULFUR_STEM_BY_OXO,
)
from namepredict.layer1.analyzer import _alkoxy_c_of, _double_bonded_o_idxs
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, _anchors
from namepredict.layer2.principal import PrincipalGroupSelection, feature_spec
from namepredict.layer1.ring_systems import sssr_rings, kekulized



class PrincipalRelation(str, Enum):
    """主基团相对骨架的位置关系（骨架内/环外）。"""
    IN_SKELETON = "in_skeleton"
    EXOCYCLIC = "exocyclic"


class PrincipalChargeState(str, Enum):
    """主基团的电荷状态（中性/全阴离子/阴离子混合）。"""
    NEUTRAL = "neutral"
    ANION = "anion"
    MIXED = "mixed"


@dataclass(frozen=True)
class PrincipalExpressionFacts:
    """主基团表达事实：类别、个数、与骨架关系、特征/锚点/附着原子集与电荷态。"""
    group_class: FunctionalGroupClass
    multiplicity: int
    relation: PrincipalRelation
    occurrence_ids: tuple[str, ...]
    characteristic_atoms: frozenset[int]
    anchor_atoms: frozenset[int]  # 官能团原锚点（occurrence.parent_anchors）
    attachment_atoms: frozenset[int]  # 骨架内附着原子：骨架外的锚点取其骨架内邻居（exocyclic）
    charge_state: PrincipalChargeState

_FG_CLASSES = frozenset(FunctionalGroupClass) - {FunctionalGroupClass.NONE}  # 全部注册 FG 类别（NONE = 纯烃）

def _anchor_fields(group_class: FunctionalGroupClass) -> tuple[str, str] | None:
    """取基团类的 (单, 复数) anchor 字段名。"""
    if group_class is FunctionalGroupClass.NONE:
        return "none_c_idx", "none_c_idxs"
    spec = feature_spec(group_class)
    return spec.anchor_fields if spec else None


_SEMANTIC_ANCHOR_FGS = frozenset({FunctionalGroupClass.RADICAL, FunctionalGroupClass.ACYL})  # 位次由 P-14.4(a) 固定 locant 1 直读 parent 字段


def _semantic_anchor_fields(group_class: FunctionalGroupClass, anchors: list[int]) -> dict:
    """P-14.4(a) 固定 locant 1 锚点的显式字段（仅单锚点写）。"""
    if group_class not in _SEMANTIC_ANCHOR_FGS or len(anchors) != 1:
        return {}
    return {_anchor_fields(group_class)[0]: anchors[0]}


def _covered(selection: PrincipalGroupSelection, skeleton: ParentSkeleton):
    """取骨架覆盖的 occurrence 子集。"""
    ids = skeleton.covered_principal_ids
    return tuple(o for o in selection.occurrences if o.id in ids)


def _chain_kind(group_class: FunctionalGroupClass, count: int) -> str | None:
    """按链 FG 类别与个数决定 kind（不支持时 None）。"""
    if group_class is FunctionalGroupClass.NONE:
        return group_class.value if count == 0 else None
    if group_class is FunctionalGroupClass.ACYL:
        return "acyl"  # 酰基残基：羰基头为 locant 1（P-65.1.7.2）
    if group_class is FunctionalGroupClass.RADICAL:
        return "radical"  # 自由基连接点位次由 L4 radical_c_idx 承载
    return group_class.value if count >= 1 else None


def _is_anion_occurrence(occurrence, mol) -> bool:
    """occurrence 是否为阴离子：周边原子带负形式电荷。"""
    if mol is None:
        return False
    return any(mol.GetAtomWithIdx(i).GetFormalCharge() < 0 for i in occurrence.payload.get("surr_idx") or ())


def _expression_flags(selection: PrincipalGroupSelection, occurrences, mol) -> dict:
    """主基团表达标志（酸全阴离子时加 anion）。"""
    if selection.group_class is not FunctionalGroupClass.ACID:
        return {}
    return {"anion": True} if occurrences and all(_is_anion_occurrence(o, mol) for o in occurrences) else {}


def _charge_state(occurrences, mol) -> PrincipalChargeState:
    """按 occurrence 阴离子情况推断电荷状态。"""
    charges = [_is_anion_occurrence(o, mol) for o in occurrences]
    if charges and all(charges):
        return PrincipalChargeState.ANION
    return PrincipalChargeState.MIXED if any(charges) else PrincipalChargeState.NEUTRAL


def _skeletal_attachments(mol, skeleton, occurrences) -> frozenset[int]:
    """计算主官能团在骨架内的附着原子（骨架外则取邻居）。"""
    atoms = set(skeleton.atom_ids)
    anchors = {i for o in occurrences for i in o.parent_anchors}
    included = anchors & atoms
    if included or mol is None:
        return frozenset(included or anchors)
    return frozenset(n.GetIdx() for i in anchors for n in mol.GetAtomWithIdx(i).GetNeighbors()
                     if n.GetIdx() in atoms)


def _facts(selection, skeleton, occurrences, mol=None) -> PrincipalExpressionFacts:
    """汇总主基团表达式 facts（关系/锚点/附着/电荷）。"""
    characteristic = frozenset(i for o in occurrences for i in o.characteristic_atoms)
    relation = PrincipalRelation.IN_SKELETON if characteristic & set(skeleton.atom_ids) else PrincipalRelation.EXOCYCLIC
    anchors = frozenset(i for o in occurrences for i in o.parent_anchors)
    attachment = _skeletal_attachments(mol, skeleton, occurrences)
    return PrincipalExpressionFacts(selection.group_class, len(occurrences), relation,
                                    tuple(o.id for o in occurrences), characteristic, anchors, attachment,
                                    _charge_state(occurrences, mol))


def _parent_dict(kind: str, skeleton: ParentSkeleton, occurrences, fields: dict,
                 facts: PrincipalExpressionFacts, all_occurrences=()) -> dict:
    """按骨架构造母体 dict（含全部主基团 occurrence 供所有权判定）。"""
    return {"kind": kind, "chain": list(skeleton.atom_ids), "n_carbons": len(skeleton.atom_ids),
            "covered_principal_ids": tuple(o.id for o in occurrences),
            "principal_occurrences": all_occurrences or occurrences,
            "principal_group_count": len(occurrences), "principal_expression_facts": facts, **fields}


# 苯系保留名由 L5 苯 variant 注入；L2 只表达结构 kind。


def _ring_endocyclic_triple(mol: Mol, atoms: set[int]) -> bool:
    """骨架内是否有成环三键（有则不能按 mancude 芳香保留名表达）。"""
    from rdkit import Chem
    return any(b.GetBondType() == Chem.BondType.TRIPLE and b.IsInRing()
               and b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
               for b in mol.GetBonds())


def _generic_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    """无保留 scaffold 时的通用环 kind（统一收敛为 alkane）。"""
    mol = info["mol"]
    atoms = set(skeleton.atom_ids)
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms) and not _ring_endocyclic_triple(mol, atoms):
        n_rings = sum(1 for ring in sssr_rings(mol) if set(ring) <= atoms)  # 未注册芳香稠环 kind 收敛 alkane
        if n_rings >= 2:
            return "alkane"
        return None
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in atoms)
    return "alkane" if all_carbon else None  # 纯烃环统一 kind='alkane'（环系/不饱和度另由字段承载）


def _resolved_ring_kind(scaffold, info: dict, skeleton: ParentSkeleton) -> str | None:
    """按保留 scaffold 解析环 kind（稠环收敛为 alkane）。"""
    if scaffold and scaffold.id != "carbocycle":
        if scaffold.id == "benzene":
            return "alkane"  # 苯环 kind 收敛 alkane，环系由 scaffold_id 承载
        if scaffold.id == "fused_hetero":
            return "alkane"  # 未注册稠环：kind 收敛 alkane，身份由 fused_tree 承载
        return scaffold.id
    return _generic_ring_kind(info, skeleton)


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int, scaffold) -> str | None:
    """决定环骨架母体的 kind（radical/正交化 FG 类/结构 kind）。"""
    if selection.group_class is FunctionalGroupClass.RADICAL:
        if scaffold is None:  # 苯基取代基保留名由 L5 benzene variant 表达
            return None  # 未知杂环：L5 无 -yl 词干，显式失败
        return "radical"
    if scaffold is not None and selection.group_class in _FG_CLASSES:  # 环 + 主 FG → FG 类别 kind（命名由 L5 通用词干引擎拼接）
        kind = _chain_kind(selection.group_class, count)
        if kind is not None:
            return kind
    if scaffold is not None and selection.group_class is FunctionalGroupClass.ALDEHYDE:
        return "aldehyde"  # 环上外环 -CHO 可多个（P-66.6.1.1.3），仅环骨架放行
    return _resolved_ring_kind(scaffold, info, skeleton)


def _scaffold_fields(info: dict, skeleton: ParentSkeleton, facts=None, scaffold=None) -> dict:
    """解析并写入 scaffold 身份与表达能力字段。"""

    from namepredict.layer2.ring_expression_policy import supports_ring_expression
    if scaffold is None:
        from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
        scaffold = resolve_ring_scaffold(info, skeleton)  # print(scaffold)
    fields: dict = {}
    if scaffold:
        supported = supports_ring_expression(scaffold, facts) if facts else False
        match = None  # 保留 fused 模板匹配映射，供 L4 固定编号用
        from namepredict.layer2.ring_scaffold import _match_with_map, get_spec, hydrogenated_atoms
        spec = get_spec(scaffold.id)
        if spec and spec.retained:  # 保留模板都取 match，供算加氢位（P-31.2.2）
            hit = _match_with_map(info, skeleton.atom_ids)
            match = hit[1] if hit and hit[0] == scaffold.id else None

        fields = {"scaffold_id": scaffold.id, "scaffold_identity": scaffold,
                  "scaffold_match": match,
                  "typed_ring_expression_supported": supported}  # print({"scaffold_id": scaffold.id, "scaffold_identity": scaffold, "scaffold_match": match, "typed_ring_expression_supported": supported})
        if match:  # 加氢原子集：编号完成后由 L4 换算为 hydro 前缀位次
            hydro = hydrogenated_atoms(info["mol"], scaffold.id, match)
            if hydro:
                fields["hydro_atoms"] = hydro
    system = next((s for s in info.get("ring_systems") or []
                   if (s.get("atom_ids") or []) == list(skeleton.atom_ids)), None)  # 多环骨架的稠环拆解结构，独立于 scaffold 身份（P-25.3.2）

    if system is not None and len(system.get("sssr_indices") or ()) >= 2:
        from namepredict.layer2.fused_system import decompose_fused_system
        node = decompose_fused_system(info, system)
        if node is not None:
            fields["fused_tree"] = node
    return fields


def _ester_o_idx(mol, e: dict) -> int | None:
    """由酯条目现算酯氧索引：羰基碳上另连烷氧基碳的单键氧。"""
    center = mol.GetAtomWithIdx(int(e["center_idx"]))
    return next((n.GetIdx() for n in center.GetNeighbors()
                 if n.GetAtomicNum() == O and _alkoxy_c_of(n, center) is not None), None)


def ester_fields(info: dict, occurrences, fields: dict) -> dict:
    """取首个酯 occurrence 并写入酯字段（alkoxy_n 恒 0）。"""
    o_idx = _ester_o_idx(info["mol"], occurrences[0].payload)
    return {**fields, "o_idx": o_idx, "alkoxy_n": 0} if o_idx is not None else fields


def express_ring_principal(info: dict, selection: PrincipalGroupSelection,
                           skeleton: ParentSkeleton) -> dict | None:
    """环骨架：表达主基团并生成母体 dict（不支持返回 None）。"""
    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold  # 骨架原子集已定：一次识别 scaffold 供下游复用
    scaffold = resolve_ring_scaffold(info, skeleton)
    occurrences = _covered(selection, skeleton)
    kind = _ring_kind(info, selection, skeleton, len(occurrences), scaffold)
    if kind is None:
        return None
    facts = _facts(selection, skeleton, occurrences, info["mol"])
    fields = {**_semantic_anchor_fields(selection.group_class, _anchors(occurrences)),  # 固定 locant 1 锚点字段
              **_chain_unsat_fields(info, skeleton,
                                    _scaffold_fields(info, skeleton, facts, scaffold))}  # 补环内不饱和字段（烯/炔由 double_bond 等承载）
    if kind == "radical" and _radical_ylidene(info, occurrences):  # 环上碳锚点自由价双键（*=C1CCCC1）：链引擎出 -ylidene
        fields = {**fields, "radical_ylidene": True}
    if facts.group_class is FunctionalGroupClass.ACID:
        fields = {**fields, **_expression_flags(selection, occurrences, info.get("mol"))}  # 环酸全阴离子补 anion 标志，L5 据此转 -ate
    if facts.group_class is FunctionalGroupClass.ESTER and facts.multiplicity == 1:
        fields = ester_fields(info, occurrences, fields)
    if facts.group_class is FunctionalGroupClass.ACYL_HALIDE and facts.multiplicity == 1:
        fields = _chain_acyl_halide_fields(info, occurrences, fields)  # 环外酰卤同样要卤素字段（hal_z/hal_idx）
    return _parent_dict(kind, skeleton, occurrences, fields, facts, selection.occurrences)


def _chain_fields(selection, occurrences, mol) -> dict:
    """构造链主基团的 anchor 与表达标志字段。"""
    return {**_semantic_anchor_fields(selection.group_class, _anchors(occurrences)),
            **_expression_flags(selection, occurrences, mol)}


def _unsat_bond_fields(dbs: list[dict], tbs: list[dict]) -> dict:
    """不饱和度 → 双键/三键字段字典（单数单键、列表多键，烯/炔可共存）。"""
    fields = {}
    if len(dbs) == 1:
        fields["double_bond"] = (dbs[0]["c1"], dbs[0]["c2"])
    elif len(dbs) >= 2:
        fields["double_bonds"] = [(d["c1"], d["c2"]) for d in dbs]
    if len(tbs) == 1:
        fields["triple_bond"] = (tbs[0]["c1"], tbs[0]["c2"])
    elif len(tbs) >= 2:
        fields["triple_bonds"] = [(t["c1"], t["c2"]) for t in tbs]
    return fields


def _chain_unsat_fields(info: dict, skeleton: ParentSkeleton, fields: dict) -> dict:
    """为骨架内 C=C/C≡C 附加双/三键位次字段（mancude 环内略过）。"""
    atom_set = set(skeleton.atom_ids)
    dbs, tbs = _chain_polys(info, atom_set)
    if fields.get("scaffold_id") == "carbocycle":  # 无保留 mancude 名兜底的碳环：环内 C=C 由 Kekulé 补回
        dbs = dbs + _kekule_ring_dbs(info, atom_set, dbs)
    implied = _implied_ring_atoms(fields, atom_set)
    if implied:
        dbs = [d for d in dbs if not (d["c1"] in implied and d["c2"] in implied)]
        tbs = [t for t in tbs if not (t["c1"] in implied and t["c2"] in implied)]
    return {**fields, **_unsat_bond_fields(dbs, tbs)}


def _kekule_ring_dbs(info: dict, atom_set: set[int], known: list[dict]) -> list[dict]:
    """补回碳环环内被芳香感知剔除的 C=C（无保留名兜底）。"""
    from rdkit import Chem

    kek = kekulized(info["mol"])
    if kek is None:  # Kekulize 失败：环内 C=C 无法补回，放弃
        return []
    have = {frozenset((d["c1"], d["c2"])) for d in known}
    out = []
    for b in kek.GetBonds():
        a, z = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (b.GetBondType() == Chem.BondType.DOUBLE and b.IsInRing()
                and a in atom_set and z in atom_set and frozenset((a, z)) not in have):
            out.append({"c1": a, "c2": z})
    return out


def _implied_ring_atoms(fields: dict, atom_set: set[int]) -> frozenset[int]:
    """保留 mancude 母体覆盖的分子原子集（多重键由母体名隐含）。"""
    from namepredict.layer2.ring_scaffold import mancude_ring_atoms
    implied = mancude_ring_atoms(fields.get("scaffold_id"), fields.get("scaffold_match"))
    if implied:
        return implied
    return frozenset(atom_set) if fields.get("fused_tree") is not None else frozenset()


def _chain_phosphate_fields(info: dict, occurrences, fields: dict) -> dict | None:
    """L5 磷酸命名的计数与盐元数据，门控不过返回 None。"""
    if len(occurrences) != 1:
        return None
    payload = occurrences[0].payload
    n_oh, n_om = int(payload.get("n_oh", 0)), int(payload.get("n_om", 0))
    salt = dict(info.get("salt") or {})
    if n_om > 0:
        if salt.get("metal") and int(salt.get("n_metal") or 0) != n_om:
            return None
    elif salt.get("metal"):
        return None
    return {**fields, "n_oh": n_oh, "n_om": n_om,
            "n_arms": int(payload.get("n_arms", 0)), "salt_meta": salt or None}


def _chain_ester_fields(info: dict, occurrences, fields: dict) -> dict:
    """L5 酯命名的酯氧侧字段：o_idx 供 o_side 识别。"""
    if not occurrences:
        return fields
    o_idx = _ester_o_idx(info["mol"], occurrences[0].payload)
    if o_idx is None:
        return fields
    if len(occurrences) == 1:
        return {**fields, "o_idx": o_idx, "alkoxy_n": 0}
    return {**fields, "o_idx": o_idx}

def _anchor_free_double(mol: Mol, idx: int) -> bool:
    """锚点原子与 `*` 虚拟原子之间的键是否为双键。"""
    from rdkit.Chem import BondType

    for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
        if nb.GetAtomicNum() == 0:
            b = mol.GetBondBetweenAtoms(idx, nb.GetIdx())
            return b is not None and b.GetBondType() == BondType.DOUBLE
    return False


def _anchor_oxo_count(mol: Mol, idx: int) -> int:
    """锚点原子上双键氧（=O）的个数，用于判定高价态硫的词干。"""
    return len(_double_bonded_o_idxs(mol.GetAtomWithIdx(idx)))


def _mononuclear_radical(info: dict, skeleton: ParentSkeleton,
                         occurrences) -> tuple[ParentSkeleton, dict] | None:
    """杂原子锚点自由基收敛为单核氢化物骨架（表 2.1），仅支持单锚点。"""
    mol = info["mol"]
    anchors = sorted({i for o in occurrences for i in o.parent_anchors})
    if len(anchors) != 1:
        return None
    stem_en = MONONUCLEAR_BY_ELEMENT.get(mol.GetAtomWithIdx(anchors[0]).GetAtomicNum())
    if stem_en is None:
        return None
    element = MONONUCLEAR_HYDRIDES[stem_en][0]
    if element == "S":  # 硫的氧化态并入词干（sulfane/sulfinyl/sulfonyl）
        stem_en = SULFUR_STEM_BY_OXO.get(_anchor_oxo_count(mol, anchors[0]), stem_en)
    elif element == "N":  # 自由价键级并入词干（azane/imine）
        stem_en = NITROGEN_STEM_BY_FREE_DOUBLE[_anchor_free_double(mol, anchors[0])]
    elif element == "P":  # 磷的氧化态并入词干（P-67.1.4.1.1.2）
        stem_en = PHOSPHORUS_STEM_BY_OXO.get(_anchor_oxo_count(mol, anchors[0]))
        if stem_en is None:  # 非 0/1 个 =O（如二氧代磷烷）无对应酰基词干，明确失败
            return None
    stem_zh = MONONUCLEAR_HYDRIDES[stem_en][1]
    new = replace(skeleton, atom_ids=(anchors[0],))
    return new, {"radical_anchor_element": element,
                 "stem_en": stem_en, "stem_zh": stem_zh}


def _radical_ylidene(info: dict, occurrences) -> bool:
    """碳锚点自由价是否为双键（*=C< ylidene，出 -ylidene）。"""
    mol = info.get("mol")
    anchors = sorted({i for o in occurrences for i in o.parent_anchors})
    if mol is None or len(anchors) != 1:
        return False
    return _anchor_free_double(mol, anchors[0])


def _chain_acyl_halide_fields(info: dict, occurrences, fields: dict) -> dict:
    """酰卤的卤素字段：hal_idx 供母体纳入卤素原子。"""
    if len(occurrences) != 1:
        return fields
    mol = info.get("mol")
    if mol is None:
        return fields
    hal = next((i for i in occurrences[0].payload["surr_idx"]
                if mol.GetAtomWithIdx(i).GetAtomicNum() in HALO_Z), None)
    if hal is None:
        return fields
    return {**fields, "hal_idx": hal, "hal_z": mol.GetAtomWithIdx(hal).GetAtomicNum()}


def express_chain_principal(info: dict, selection: PrincipalGroupSelection,
                            skeleton: ParentSkeleton) -> dict | None:
    """链骨架：表达主基团并生成母体 dict（不支持返回 None）。"""
    if skeleton.topology is not SkeletonTopology.ACYCLIC:
        return None

    occurrences = _covered(selection, skeleton)
    kind = _chain_kind(selection.group_class, len(occurrences))

    if kind is None:
        return None
    extra: dict = {}
    if kind == "radical":
        mono = _mononuclear_radical(info, skeleton, occurrences)
        if mono is not None:
            skeleton, extra = mono
        elif _radical_ylidene(info, occurrences):  # 碳锚点自由价双键（*=C<）：链引擎出 -ylidene
            extra = {"radical_ylidene": True}
    fields = _chain_unsat_fields(info, skeleton, {**_chain_fields(selection, occurrences, info.get("mol")), **extra})
    if kind == "ester":
        fields = _chain_ester_fields(info, occurrences, fields)
    elif kind == "acyl_halide":
        fields = _chain_acyl_halide_fields(info, occurrences, fields)
    elif kind == "phosphate":
        fields = _chain_phosphate_fields(info, occurrences, fields)
        if fields is None:
            return None
    return _parent_dict(kind, skeleton, occurrences, fields,
                        _facts(selection, skeleton, occurrences, info.get("mol")), selection.occurrences)


# ── 无主官能团（纯烃）表达：P-44.1 缺位按拓扑分配 kind ──

def _chain_polys(info: dict, atom_set: set[int]) -> tuple[list[dict], list[dict]]:
    """骨架内的 C=C / C≡C 条目（端点都在 atom_set 中）。"""
    dbs = [d for d in info.get("double_bonds") or [] if d["c1"] in atom_set and d["c2"] in atom_set]
    tbs = [t for t in info.get("triple_bonds") or [] if t["c1"] in atom_set and t["c2"] in atom_set]
    return dbs, tbs