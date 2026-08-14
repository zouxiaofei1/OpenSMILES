"""Annotate selected skeletons with principal-group expression facts."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology, _anchors
from namepredict.layer2.principal import PrincipalGroupSelection


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


# 数量派生 kind 已统一:multiplicity 由 principal_expression_facts 承载,
# acid/alcohol/amine 对任意数量恒用基团名;仅 KETONE 保留 dione 区分环酮表达。
_CHAIN_KINDS = {
    FunctionalGroupClass.ACID: {1: "acid"},
    FunctionalGroupClass.ESTER: {1: "ester"},
    FunctionalGroupClass.KETONE: {1: "ketone", 2: "dione"},
    FunctionalGroupClass.ALDEHYDE: {1: "aldehyde"},
    FunctionalGroupClass.NITRILE: {1: "nitrile"},
    FunctionalGroupClass.AMIDE: {1: "amide"},
    FunctionalGroupClass.ALCOHOL: {1: "alcohol"},
    FunctionalGroupClass.AMINE: {1: "amine"},
    FunctionalGroupClass.NONE: {0: "alkane"}
}
_FIELDS = {
    FunctionalGroupClass.RADICAL: ("radical_c_idx", "radical_c_idxs"),
    FunctionalGroupClass.ACID: ("cooh_c_idx", "cooh_c_idxs"),
    FunctionalGroupClass.ESTER: ("ester_c_idx", "ester_c_idxs"),
    FunctionalGroupClass.KETONE: ("ketone_c_idx", "ketone_c_idxs"),
    FunctionalGroupClass.ALDEHYDE: ("aldehyde_c_idx", "aldehyde_c_idxs"),
    FunctionalGroupClass.NITRILE: ("nitrile_c_idx", "nitrile_c_idxs"),
    FunctionalGroupClass.AMIDE: ("amide_c_idx", "amide_c_idxs"),
    FunctionalGroupClass.ALCOHOL: ("oh_c_idx", "oh_c_idxs"),
    FunctionalGroupClass.AMINE: ("amine_c_idx", "amine_c_idxs"),
     FunctionalGroupClass.NONE:("none_c_idx","none_c_idxs"),
}


def _covered(selection: PrincipalGroupSelection, skeleton: ParentSkeleton):
    ids = skeleton.covered_principal_ids
    return tuple(o for o in selection.occurrences if o.id in ids)


def _chain_kind(group_class: FunctionalGroupClass, count: int) -> str | None:
    if group_class in (FunctionalGroupClass.ACID, FunctionalGroupClass.ALCOHOL,
                       FunctionalGroupClass.AMINE):
        # 任意主基团数(≥1)→ 基团名;count=0(仲/叔胺等未覆盖)不产 facts。
        return _CHAIN_KINDS[group_class][1] if count >= 1 else None
    kinds = _CHAIN_KINDS.get(group_class)
    if kinds is None:
        return None
    return kinds.get(count)


def _principal_fields(group_class: FunctionalGroupClass, anchors: list[int]) -> dict:
    single, plural = _FIELDS[group_class]
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


def _carbocycle_kind(mol, skeleton: ParentSkeleton) -> str:
    atoms = set(skeleton.atom_ids)
    unsaturated = sum(
        bond.GetBondTypeAsDouble() > 1
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms
    )
    return "cycloalkane" if not unsaturated else "cycloalkene" if unsaturated == 1 else "cyclopolyene"


def _ring_ketone_kind(scaffold_id: str | None, selection, skeleton: ParentSkeleton, count: int) -> str | None:
    if selection.group_class is not FunctionalGroupClass.KETONE:
        return None
    if scaffold_id is not None and scaffold_id not in {"cycloalkane", "cycloketone", "cycloalkanedione"}:
        return None
    anchors = set(_anchors(selection.occurrences))
    if not anchors or not anchors <= set(skeleton.atom_ids):
        return None
    return "cycloketone" if count == 1 else "cycloalkanedione" if count == 2 else None


def _generic_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    mol = info["mol"]
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in skeleton.atom_ids):
        return None
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids)
    return _carbocycle_kind(mol, skeleton) if all_carbon else None


def _resolved_ring_kind(scaffold, info: dict, skeleton: ParentSkeleton) -> str | None:
    return scaffold.id if scaffold and scaffold.id != "carbocycle" else _generic_ring_kind(info, skeleton)


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int, scaffold) -> str | None:
    ketone = _ring_ketone_kind(scaffold.id if scaffold else None, selection, skeleton, count)
    if ketone:
        return ketone
    if selection.group_class is FunctionalGroupClass.AMINE and count != 1:
        return None
    # 苯基取代基 radical 保持 'phenyl'（P-22.2.4），非保留名 parent；
    # 其余苯系保留名（benzoic/phenol/...）由 L5 typed_kinds 决定。
    if selection.group_class is FunctionalGroupClass.RADICAL and _is_benzene(info, skeleton):
        return "phenyl"
    return _resolved_ring_kind(scaffold, info, skeleton)


def _ring_fields(selection: PrincipalGroupSelection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    fields = _principal_fields(selection.group_class, anchors) if selection.group_class in _FIELDS else {}
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
    from namepredict.layer2.arene_carbonyl import _benzoate_alkoxy, _benzoate_side_fields
    e = info["esters"][0]
    side = _benzoate_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"]) or {}
    return {**fields, "o_idx": e["o_idx"], "alkoxy_c_idx": e["alkoxy_c_idx"],
            **_benzoate_side_fields(side)}


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
    fields = {**_ring_fact_fields(_ring_fields(selection, occurrences), facts),
              **_scaffold_fields(info, skeleton, facts, scaffold)}
    if facts.group_class is FunctionalGroupClass.ESTER and facts.multiplicity == 1:
        fields = _benzoate_ester_fields(info, fields)
    return _parent_dict(kind, skeleton, occurrences, fields, facts)


def _chain_fields(selection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    return {**_principal_fields(selection.group_class, anchors),
            **_expression_flags(selection, occurrences)}


def _unsat_bond_fields(dbs: list[dict], tbs: list[dict]) -> dict:
    """Unsaturation → double/triple-bond field dict (empty for mixed/none)."""
    if len(tbs) == 1 and not dbs:
        return {"triple_bond": (tbs[0]["c1"], tbs[0]["c2"])}
    if len(dbs) == 1 and not tbs:
        return {"double_bond": (dbs[0]["c1"], dbs[0]["c2"])}
    if len(dbs) >= 2 and not tbs:
        return {"double_bonds": [(d["c1"], d["c2"]) for d in dbs]}
    return {}


def _chain_unsat_fields(info: dict, skeleton: ParentSkeleton, fields: dict) -> dict:
    """Attach double/triple-bond locant fields for C=C/C≡C inside the skeleton."""
    dbs, tbs = _chain_polys(info, set(skeleton.atom_ids))
    return {**fields, **_unsat_bond_fields(dbs, tbs)}


def _chain_ester_fields(info: dict, occurrences, fields: dict) -> dict:
    """Ester alkoxy-side fields (o_idx/alkoxy_c_idx/alkoxy_*) for L5 ester naming."""
    if len(occurrences) != 1:
        return fields
    match = next((e for e in (info.get("esters") or [])
                  if e["c_idx"] in occurrences[0].characteristic_atoms), None)
    if match is None:
        return fields
    from namepredict.tools.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(info["mol"], match["o_idx"], match["alkoxy_c_idx"])
    return {**fields, "o_idx": match["o_idx"], "alkoxy_c_idx": match["alkoxy_c_idx"], **side}


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
    from namepredict.layer2.naphthalene import _naph_chains, _naph_parent_dict
    chains = _naph_chains(info)
    if chains and set(chains[0]) == atoms:
        return {**_naph_parent_dict(info, "naphthalene"), **_scaffold_fields(info, skeleton, None, scaffold)}
    return None


def _saturated_ring_parent(info: dict, skeleton: ParentSkeleton, atoms: set[int]) -> dict | None:
    """非芳香环：按环内不饱和度分配 cycloalkane/cycloalkene/cyclopolyene。"""
    dbs, tbs = _chain_polys(info, atoms)
    if tbs:
        return None  # 环炔暂不支持
    chain = _mono_ring_chain(info, atoms) or list(skeleton.atom_ids)
    bf = _unsat_bond_fields(dbs, tbs)  # tbs 已排除：只可能 double_bond/double_bonds/空
    if not dbs:
        kind = "cycloalkane"
    elif len(dbs) == 1:
        kind = "cycloalkene"
    else:
        kind = "cyclopolyene"
    return {"kind": kind, "chain": chain, "n_carbons": len(chain), **bf}


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
