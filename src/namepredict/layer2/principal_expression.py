"""为选中的骨架标注主基团表达式 facts。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, _anchors
from namepredict.layer2.principal import PrincipalGroupSelection, feature_spec


class PrincipalRelation(str, Enum):
    IN_SKELETON = "in_skeleton"
    EXOCYCLIC = "exocyclic"


class PrincipalChargeState(str, Enum):
    NEUTRAL = "neutral"
    ANION = "anion"
    MIXED = "mixed"


@dataclass(frozen=True)
class PrincipalExpressionFacts:
    group_class: FunctionalGroupClass
    multiplicity: int
    relation: PrincipalRelation
    occurrence_ids: tuple[str, ...]
    characteristic_atoms: frozenset[int]
    attachment_atoms: frozenset[int]
    charge_state: PrincipalChargeState


# 链式主官能团：kind 恒为 FG 类别名；acid/alcohol/amine/ketone 任意数量恒用基团名，其余链 FG 仅单基。
_CHAIN_FG = frozenset({
    FunctionalGroupClass.ACID, FunctionalGroupClass.ALCOHOL,
    FunctionalGroupClass.AMINE, FunctionalGroupClass.KETONE,
    FunctionalGroupClass.ESTER, FunctionalGroupClass.AMIDE,
    FunctionalGroupClass.NITRILE, FunctionalGroupClass.ALDEHYDE,
})
_MULTI_FG = frozenset({
    FunctionalGroupClass.ACID, FunctionalGroupClass.ALCOHOL,
    FunctionalGroupClass.AMINE, FunctionalGroupClass.KETONE,
})
def _anchor_fields(group_class: FunctionalGroupClass) -> tuple[str, str] | None:
    if group_class is FunctionalGroupClass.NONE:
        return "none_c_idx", "none_c_idxs"
    spec = feature_spec(group_class)
    return spec.anchor_fields if spec else None


def _covered(selection: PrincipalGroupSelection, skeleton: ParentSkeleton):
    ids = skeleton.covered_principal_ids
    return tuple(o for o in selection.occurrences if o.id in ids)


def _chain_kind(group_class: FunctionalGroupClass, count: int) -> str | None:
    if group_class is FunctionalGroupClass.NONE:
        return group_class.value if count == 0 else None
    if group_class not in _CHAIN_FG:
        return None
    if group_class not in _MULTI_FG and count != 1:
        return None
    return group_class.value if count >= 1 else None


def _principal_fields(group_class: FunctionalGroupClass, anchors: list[int]) -> dict:
    single, plural = _anchor_fields(group_class)
    return {single: anchors[0], plural: anchors} if len(anchors) == 1 else {plural: anchors}


def _expression_flags(selection: PrincipalGroupSelection, occurrences) -> dict:
    if selection.group_class is not FunctionalGroupClass.ACID:
        return {}
    return {"anion": True} if occurrences and all(o.payload.get("anion") for o in occurrences) else {}


def _charge_state(occurrences) -> PrincipalChargeState:
    charges = [bool(o.payload.get("anion")) for o in occurrences]
    if charges and all(charges):
        return PrincipalChargeState.ANION
    return PrincipalChargeState.MIXED if any(charges) else PrincipalChargeState.NEUTRAL


def _skeletal_attachments(mol, skeleton, occurrences) -> frozenset[int]:
    atoms = set(skeleton.atom_ids)
    anchors = {i for o in occurrences for i in o.parent_anchors}
    included = anchors & atoms
    if included or mol is None:
        return frozenset(included or anchors)
    return frozenset(n.GetIdx() for i in anchors for n in mol.GetAtomWithIdx(i).GetNeighbors()
                     if n.GetIdx() in atoms)


def _facts(selection, skeleton, occurrences, mol=None) -> PrincipalExpressionFacts:
    characteristic = frozenset(i for o in occurrences for i in o.characteristic_atoms)
    relation = PrincipalRelation.IN_SKELETON if characteristic & set(skeleton.atom_ids) else PrincipalRelation.EXOCYCLIC
    attachment = _skeletal_attachments(mol, skeleton, occurrences)
    return PrincipalExpressionFacts(selection.group_class, len(occurrences), relation,
                                    tuple(o.id for o in occurrences), characteristic, attachment,
                                    _charge_state(occurrences))


def _parent_dict(kind: str, skeleton: ParentSkeleton, occurrences, fields: dict,
                 facts: PrincipalExpressionFacts) -> dict:
    return {"kind": kind, "chain": list(skeleton.atom_ids), "n_carbons": len(skeleton.atom_ids),
            "covered_principal_ids": tuple(o.id for o in occurrences),
            "principal_group_count": len(occurrences), "principal_expression_facts": facts, **fields}


# 苯系保留名已全部迁往 L5 typed_kinds（_BENZENE_RETAINED）；L2 只表达结构 kind。
def _is_benzene(info: dict, skeleton: ParentSkeleton) -> bool:
    mol = info["mol"]
    return len(skeleton.atom_ids) == 6 and all(
        mol.GetAtomWithIdx(i).GetAtomicNum() == 6 and mol.GetAtomWithIdx(i).GetIsAromatic()
        for i in skeleton.atom_ids)


def _generic_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    mol = info["mol"]
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in skeleton.atom_ids):
        return None
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids)
    # 纯烃环统一 kind='alkane'（正交化：环系由 scaffold_id 承载，不饱和度由 double_bond/double_bonds 字段承载，命名由 chain_engine 动态加 cyclo 前缀）。
    return "alkane" if all_carbon else None


def _resolved_ring_kind(scaffold, info: dict, skeleton: ParentSkeleton) -> str | None:
    if scaffold and scaffold.id != "carbocycle":
        # 苯环（纯烃芳香单环）kind 收敛为 alkane，环系由 scaffold_id="benzene" 承载（对齐环烷烃正交化）。
        return "alkane" if scaffold.id == "benzene" else scaffold.id
    return _generic_ring_kind(info, skeleton)


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int, scaffold) -> str | None:
    if selection.group_class is FunctionalGroupClass.AMINE and count != 1:
        return None
    # 苯基取代基 radical 保持 'phenyl'（P-22.2.4）；环 + 主 FG 的 kind 收敛为 FG 类别（正交化），环骨架由 scaffold_id 承载，命名 kind（cycloalcohol/cycloketone/benzoic/phenol/...）由 L5 typed_kinds 决定。
    if selection.group_class is FunctionalGroupClass.RADICAL and _is_benzene(info, skeleton):
        return "phenyl"
    # 环（饱和环/苯环）+ 主 FG → FG 类别 kind（正交化）；命名 kind（cycloalcohol/cycloketone/cycloalkanecarboxylic/benzoic/phenol/...）由 L5 typed_kinds 决定；稠环（naphthalene 等）暂保持结构 kind。
    if scaffold and scaffold.id in ("carbocycle", "benzene"):
        if selection.group_class in (FunctionalGroupClass.ALCOHOL,
                                     FunctionalGroupClass.KETONE,
                                     FunctionalGroupClass.AMINE,
                                     FunctionalGroupClass.ACID,
                                     FunctionalGroupClass.ALDEHYDE,
                                     FunctionalGroupClass.NITRILE,
                                     FunctionalGroupClass.AMIDE,
                                     FunctionalGroupClass.ESTER):
            kind = _chain_kind(selection.group_class, count)
            if kind is not None:
                return kind
    return _resolved_ring_kind(scaffold, info, skeleton)


def _ring_fields(selection: PrincipalGroupSelection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    fields = _principal_fields(selection.group_class, anchors) if _anchor_fields(selection.group_class) else {}
    characteristic = sorted({i for o in occurrences for i in o.characteristic_atoms})
    return {**fields, "principal_characteristic_atoms": characteristic,
            "principal_attachment_atoms": anchors,
            "principal_group_class": selection.group_class.value}


def _ring_fact_fields(fields: dict, facts: PrincipalExpressionFacts) -> dict:
    attachments = sorted(facts.attachment_atoms)
    extra = {"principal_attachment_atoms": attachments}
    if len(attachments) == 1 and facts.group_class in (
        FunctionalGroupClass.ACID, FunctionalGroupClass.ESTER,
    ):
        extra["ring_attach_idx"] = attachments[0]
    return {**fields, **extra}


def _scaffold_fields(info: dict, skeleton: ParentSkeleton, facts=None, scaffold=None) -> dict:
    from namepredict.layer2.ring_expression_policy import supports_ring_expression
    if scaffold is None:
        from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
        scaffold = resolve_ring_scaffold(info, skeleton)
    if not scaffold:
        return {}
    supported = supports_ring_expression(scaffold, facts) if facts else False
    return {"scaffold_id": scaffold.id, "scaffold_identity": scaffold,
            "typed_ring_expression_supported": supported}


def _benzoate_ester_fields(info: dict, fields: dict) -> dict:
    e = info["esters"][0]
    return {**fields, "o_idx": e["o_idx"], "alkoxy_n": 0}


def express_ring_principal(info: dict, selection: PrincipalGroupSelection,
                           skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    # 骨架原子集已确定：一次识别 scaffold，下游复用（不再重复调 _producer_id）。
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
    scaffold = resolve_ring_scaffold(info, skeleton)
    occurrences = _covered(selection, skeleton)
    kind = _ring_kind(info, selection, skeleton, len(occurrences), scaffold)
    if kind is None:
        return None
    facts = _facts(selection, skeleton, occurrences, info["mol"])
    # 补环内不饱和字段：kind 正交化后（醇/酮/纯烃环 → FG 类别/alkane），烯/炔由 double_bond(s)/triple_bond 字段承载（否则环烯酮/环烯醇/环烯烃烯丢失）。
    fields = {**_ring_fact_fields(_ring_fields(selection, occurrences), facts),
              **_chain_unsat_fields(info, skeleton,
                                    _scaffold_fields(info, skeleton, facts, scaffold))}
    if facts.group_class is FunctionalGroupClass.ESTER and facts.multiplicity == 1:
        fields = _benzoate_ester_fields(info, fields)
    return _parent_dict(kind, skeleton, occurrences, fields, facts)


def _chain_fields(selection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    return {**_principal_fields(selection.group_class, anchors),
            **_expression_flags(selection, occurrences)}


def _unsat_bond_fields(dbs: list[dict], tbs: list[dict]) -> dict:
    """不饱和度 → 双键/三键字段字典（混合/无时为空）。"""
    if len(tbs) == 1 and not dbs:
        return {"triple_bond": (tbs[0]["c1"], tbs[0]["c2"])}
    if len(dbs) == 1 and not tbs:
        return {"double_bond": (dbs[0]["c1"], dbs[0]["c2"])}
    if len(dbs) >= 2 and not tbs:
        return {"double_bonds": [(d["c1"], d["c2"]) for d in dbs]}
    return {}


def _chain_unsat_fields(info: dict, skeleton: ParentSkeleton, fields: dict) -> dict:
    """为骨架内的 C=C/C≡C 附加双键/三键位次字段。"""
    dbs, tbs = _chain_polys(info, set(skeleton.atom_ids))
    return {**fields, **_unsat_bond_fields(dbs, tbs)}


def _chain_ester_fields(info: dict, occurrences, fields: dict) -> dict:
    """L5 酯命名的酯烷氧基侧字段：o_idx 供 o_side 识别；仅严格线性给 alkoxy_n 保留名。"""
    if len(occurrences) != 1:
        return fields
    match = next((e for e in (info.get("esters") or [])
                  if e["c_idx"] in occurrences[0].characteristic_atoms), None)
    if match is None:
        return fields
    return {**fields, "o_idx": match["o_idx"], "alkoxy_n": 0}


def express_chain_principal(info: dict, selection: PrincipalGroupSelection,
                            skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.ACYCLIC:
        return None
    
    occurrences = _covered(selection, skeleton)
    kind = _chain_kind(selection.group_class, len(occurrences))

    if kind is None:
        return None
    fields = _chain_unsat_fields(info, skeleton, _chain_fields(selection, occurrences))
    if kind == "ester":
        fields = _chain_ester_fields(info, occurrences, fields)
    return _parent_dict(kind, skeleton, occurrences, fields,
                        _facts(selection, skeleton, occurrences))


# ── 无主官能团（纯烃）表达：P-44.1 缺位时按拓扑分配 hydrocarbon kind ──

def _system_is_aromatic(info: dict, atoms: set[int]) -> bool:
    return any(set(s.get("atom_ids") or ()) == atoms and s.get("is_aromatic_mancude")
               for s in info.get("ring_systems") or [])


def _mono_ring_chain(info: dict, atoms: set[int]) -> list[int] | None:
    """SSSR 环序（合法顺序），供单环骨架使用。"""
    return next((list(r["atom_ids"]) for r in info.get("rings") or []
                 if set(r["atom_ids"]) == atoms), None)


def _chain_polys(info: dict, atom_set: set[int]) -> tuple[list[dict], list[dict]]:
    """骨架内的 C=C / C≡C 条目（端点都在 atom_set 中）。"""
    dbs = [d for d in info.get("double_bonds") or [] if d["c1"] in atom_set and d["c2"] in atom_set]
    tbs = [t for t in info.get("triple_bonds") or [] if t["c1"] in atom_set and t["c2"] in atom_set]
    return dbs, tbs


def _hydrocarbon_chain_parent(info: dict, skeleton: ParentSkeleton) -> dict:
    chain = list(skeleton.atom_ids)
    dbs, tbs = _chain_polys(info, set(chain))
    bf = _unsat_bond_fields(dbs, tbs)
    if "triple_bond" in bf:
        kind = "alkyne"
    elif "double_bond" in bf:
        kind = "alkene"
    elif "double_bonds" in bf:
        kind = "polyene"
    else:
        kind = "alkane"
    return {"kind": kind, "chain": chain, "n_carbons": len(chain), **bf}


def _aromatic_scaffold_parent(info: dict, skeleton: ParentSkeleton, atoms: set[int]) -> dict | None:
    """芳香环：解析到保留 scaffold 则用其 id（benzene/naphthalene 等）。"""
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
    scaffold = resolve_ring_scaffold(info, skeleton)
    if scaffold is None or scaffold.id == "carbocycle":
        return None  # 芳香碳环未匹配保留 scaffold（如 anthracene）
    if scaffold.id != "naphthalene":
        return {"kind": scaffold.id, "chain": _mono_ring_chain(info, atoms) or list(skeleton.atom_ids),
                "n_carbons": len(atoms), **_scaffold_fields(info, skeleton, None, scaffold)}
    return None


def _saturated_ring_parent(info: dict, skeleton: ParentSkeleton, atoms: set[int]) -> dict | None:
    """非芳香环：kind 统一为 alkane（正交化），烯信息由 double_bond(s) 字段承载。"""
    dbs, tbs = _chain_polys(info, atoms)
    if tbs:
        return None  # 环炔暂不支持
    chain = _mono_ring_chain(info, atoms) or list(skeleton.atom_ids)
    bf = _unsat_bond_fields(dbs, tbs)  # tbs 已排除：只可能 double_bond/double_bonds/空
    return {"kind": "alkane", "chain": chain, "n_carbons": len(chain), **bf}


def _hydrocarbon_ring_parent(info: dict, skeleton: ParentSkeleton) -> dict | None:
    atoms = set(skeleton.atom_ids)
    if _system_is_aromatic(info, atoms):
        return _aromatic_scaffold_parent(info, skeleton, atoms)
    parent = _saturated_ring_parent(info, skeleton, atoms)
    return {**parent, **_scaffold_fields(info, skeleton, None)} if parent else None


def express_hydrocarbon_principal(info: dict, skeleton: ParentSkeleton) -> dict | None:
    """无主官能团时：按拓扑分配纯烃 kind（alkane/ene/yne/polyene/环/保留 scaffold）。"""
    if skeleton.topology is SkeletonTopology.RING_SYSTEM:
        return _hydrocarbon_ring_parent(info, skeleton)
    return _hydrocarbon_chain_parent(info, skeleton)
