"""Annotate selected skeletons with principal-group expression facts."""
from __future__ import annotations

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass
from namepredict.layer2.parent_skeleton import ParentSkeleton, SkeletonTopology
from namepredict.layer2.principal_selection import PrincipalGroupSelection


_CHAIN_KINDS = {
    FunctionalGroupClass.ACID: {1: "acid", 2: "diacid", 3: "polycarboxylic"},
    FunctionalGroupClass.ALCOHOL: {1: "alcohol", 2: "diol", 3: "triol"},
    FunctionalGroupClass.AMINE: {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"},
}
_FIELDS = {
    FunctionalGroupClass.ACID: ("cooh_c_idx", "cooh_c_idxs"),
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


def _parent_dict(kind: str, skeleton: ParentSkeleton, occurrences, fields: dict) -> dict:
    return {"kind": kind, "chain": list(skeleton.atom_ids), "n_carbons": len(skeleton.atom_ids),
            "covered_principal_ids": tuple(o.id for o in occurrences),
            "principal_group_count": len(occurrences), **fields}


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


def _retained_scaffold_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    from namepredict.layer2.retained_registry import match_systems
    atoms = set(skeleton.atom_ids)
    return next((sid for sid, system, _ in match_systems(info)
                 if set(system.get("atom_ids") or ()) == atoms), None)


def _carbocycle_kind(mol, skeleton: ParentSkeleton) -> str:
    atoms = set(skeleton.atom_ids)
    unsaturated = sum(
        bond.GetBondTypeAsDouble() > 1
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms
    )
    return "cycloalkane" if not unsaturated else "cycloalkene" if unsaturated == 1 else "cyclopolyene"


def _ring_kind(info: dict, selection: PrincipalGroupSelection, skeleton: ParentSkeleton, count: int) -> str | None:
    if _is_benzene(info, skeleton) and count == 1:
        return _RETAINED_RING_KINDS.get(selection.group_class)
    mol = info["mol"]
    if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in skeleton.atom_ids):
        return _retained_scaffold_id(info, skeleton)
    if all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids):
        return _carbocycle_kind(mol, skeleton)
    return None


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


def express_ring_principal(info: dict, selection: PrincipalGroupSelection,
                           skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    occurrences = _covered(selection, skeleton)
    kind = _ring_kind(info, selection, skeleton, len(occurrences))
    if kind is None:
        return None
    return _parent_dict(kind, skeleton, occurrences, _ring_fields(selection, occurrences))


def express_chain_principal(selection: PrincipalGroupSelection, skeleton: ParentSkeleton) -> dict | None:
    if skeleton.topology is not SkeletonTopology.ACYCLIC:
        return None
    occurrences = _covered(selection, skeleton)
    kind = _chain_kind(selection.group_class, len(occurrences))
    if kind is None:
        return None
    anchors = _anchors(occurrences)
    fields = {**_principal_fields(selection.group_class, anchors),
              **_expression_flags(selection, occurrences)}
    return _parent_dict(kind, skeleton, occurrences, fields)
