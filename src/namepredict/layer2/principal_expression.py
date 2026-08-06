"""Annotate selected skeletons with principal-group expression facts."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology
from namepredict.layer2.principal_selection import PrincipalGroupSelection


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


_CHAIN_KINDS = {
    FunctionalGroupClass.ACID: {1: "acid", 2: "diacid", 3: "polycarboxylic"},
    FunctionalGroupClass.KETONE: {1: "ketone", 2: "dione"},
    FunctionalGroupClass.ALCOHOL: {1: "alcohol", 2: "diol", 3: "triol"},
    FunctionalGroupClass.AMINE: {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"},
}
_FIELDS = {
    FunctionalGroupClass.ACID: ("cooh_c_idx", "cooh_c_idxs"),
    FunctionalGroupClass.KETONE: ("ketone_c_idx", "ketone_c_idxs"),
    FunctionalGroupClass.ALCOHOL: ("oh_c_idx", "oh_c_idxs"),
    FunctionalGroupClass.AMINE: ("amine_c_idx", "amine_c_idxs"),
}


def _covered(selection: PrincipalGroupSelection, skeleton: ParentSkeleton):
    ids = skeleton.covered_principal_ids
    return tuple(o for o in selection.occurrences if o.id in ids)


def _anchors(occurrences) -> list[int]:
    return sorted({a for occurrence in occurrences for a in occurrence.parent_anchors})


def _chain_kind(group_class: FunctionalGroupClass, count: int) -> str | None:
    kinds = _CHAIN_KINDS.get(group_class)
    if kinds is None or count < 1:
        return None
    return kinds.get(count, "polycarboxylic" if group_class is FunctionalGroupClass.ACID else None)


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


_RETAINED_RING_KINDS = {
    FunctionalGroupClass.ACID: "benzoic",
    FunctionalGroupClass.ALDEHYDE: "benzaldehyde",
    FunctionalGroupClass.NITRILE: "benzonitrile",
    FunctionalGroupClass.AMIDE: "benzamide",
    FunctionalGroupClass.ALCOHOL: "phenol",
    FunctionalGroupClass.AMINE: "aniline",
}
_RING_FIELDS = {
    FunctionalGroupClass.ACID: ("cooh_c_idx", "cooh_c_idxs"),
    FunctionalGroupClass.KETONE: ("ketone_c_idx", "ketone_c_idxs"),
    FunctionalGroupClass.ALDEHYDE: ("aldehyde_c_idx", "aldehyde_c_idxs"),
    FunctionalGroupClass.NITRILE: ("nitrile_c_idx", "nitrile_c_idxs"),
    FunctionalGroupClass.AMIDE: ("amide_c_idx", "amide_c_idxs"),
    FunctionalGroupClass.ALCOHOL: ("oh_c_idx", "oh_c_idxs"),
    FunctionalGroupClass.AMINE: ("amine_c_idx", "amine_c_idxs"),
}


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


def _ring_ketone_kind(selection, skeleton: ParentSkeleton, count: int) -> str | None:
    if selection.group_class is not FunctionalGroupClass.KETONE:
        return None
    if skeleton.scaffold_id and skeleton.scaffold_id not in {"cycloalkane", "cycloketone", "cycloalkanedione"}:
        return None
    anchors = {i for o in selection.occurrences for i in o.parent_anchors}
    if not anchors or not anchors <= set(skeleton.atom_ids):
        return None
    return "cycloketone" if count == 1 else "cycloalkanedione" if count == 2 else None


def _generic_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    mol = info["mol"]
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in skeleton.atom_ids):
        return None
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids)
    return _carbocycle_kind(mol, skeleton) if all_carbon else None


def _resolved_ring_kind(info: dict, skeleton: ParentSkeleton) -> str | None:
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
    scaffold = resolve_ring_scaffold(info, skeleton)
    return scaffold.id if scaffold and scaffold.id != "carbocycle" else _generic_ring_kind(info, skeleton)


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int) -> str | None:
    ketone = _ring_ketone_kind(selection, skeleton, count)
    if ketone:
        return ketone
    if selection.group_class is FunctionalGroupClass.AMINE and count != 1:
        return None
    if _is_benzene(info, skeleton) and count == 1:
        return _RETAINED_RING_KINDS.get(selection.group_class)
    return _resolved_ring_kind(info, skeleton)


def _ring_fields(selection: PrincipalGroupSelection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    names = _RING_FIELDS.get(selection.group_class)
    fields = _principal_fields(selection.group_class, anchors) if selection.group_class in _FIELDS else {}
    if names and anchors:
        fields.update({names[0]: anchors[0], names[1]: anchors} if len(anchors) == 1 else {names[1]: anchors})
    characteristic = sorted({i for o in occurrences for i in o.characteristic_atoms})
    return {**fields, "principal_characteristic_atoms": characteristic,
            "principal_attachment_atoms": anchors,
            "principal_group_class": selection.group_class.value}


def _ring_fact_fields(fields: dict, facts: PrincipalExpressionFacts) -> dict:
    attachments = sorted(facts.attachment_atoms)
    extra = {"principal_attachment_atoms": attachments}
    if facts.group_class is FunctionalGroupClass.ACID and len(attachments) == 1:
        extra["ring_attach_idx"] = attachments[0]
    return {**fields, **extra}


def _scaffold_fields(info: dict, skeleton: ParentSkeleton, facts=None) -> dict:
    from namepredict.layer2.ring_expression_policy import supports_ring_expression
    from namepredict.layer2.ring_scaffold import resolve_ring_scaffold
    scaffold = resolve_ring_scaffold(info, skeleton)
    if not scaffold:
        return {}
    supported = supports_ring_expression(scaffold, facts) if facts else False
    return {"scaffold_id": scaffold.id, "scaffold_identity": scaffold,
            "typed_ring_expression_supported": supported}


def express_ring_principal(info: dict, selection: PrincipalGroupSelection,
                           skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    occurrences = _covered(selection, skeleton)
    kind = _ring_kind(info, selection, skeleton, len(occurrences))
    if kind is None:
        return None
    facts = _facts(selection, skeleton, occurrences, info["mol"])
    fields = {**_ring_fact_fields(_ring_fields(selection, occurrences), facts),
              **_scaffold_fields(info, skeleton, facts)}
    return _parent_dict(kind, skeleton, occurrences, fields, facts)


def _chain_fields(selection, occurrences) -> dict:
    anchors = _anchors(occurrences)
    return {**_principal_fields(selection.group_class, anchors),
            **_expression_flags(selection, occurrences)}


def express_chain_principal(selection: PrincipalGroupSelection, skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.ACYCLIC:
        return None
    occurrences = _covered(selection, skeleton)
    kind = _chain_kind(selection.group_class, len(occurrences))
    if kind is None:
        return None
    return _parent_dict(kind, skeleton, occurrences, _chain_fields(selection, occurrences),
                        _facts(selection, skeleton, occurrences))
