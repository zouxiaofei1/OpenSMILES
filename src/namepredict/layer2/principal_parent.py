"""Rule-driven principal parent-skeleton selection entry point."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass, inventory_from_info
from namepredict.layer2.parent_skeleton import SkeletonSelection, SkeletonTopology, select_principal_skeletons
from namepredict.layer2.principal_expression import express_chain_principal, express_ring_principal
from namepredict.layer2.principal import PrincipalGroupSelection, select_principal_group

@dataclass(frozen=True)
class PrincipalParentSelection:
    principal: PrincipalGroupSelection | None
    skeletons: SkeletonSelection | None

def select_principal_parent_skeletons(info: dict) -> PrincipalParentSelection:
    principal = select_principal_group(inventory_from_info(info))
    occurrences = principal.occurrences if principal else ()
    skeletons = select_principal_skeletons(info, occurrences)
    return PrincipalParentSelection(principal, skeletons)

def _owned(parent: dict | None, selection: PrincipalParentSelection) -> dict | None:
    if parent is None:
        return None
    occurrences = selection.principal.occurrences
    return {**parent, "covered_principal_ids": tuple(o.id for o in occurrences),
            "principal_group_count": len(occurrences)}

def _special_expression(selection: PrincipalParentSelection, info: dict) -> dict | None:

    return None

def _unsupported_typed_ring(parent: dict, selection: PrincipalParentSelection) -> bool:
    return (selection.principal.group_class is FunctionalGroupClass.KETONE
            and parent.get("scaffold_identity") is not None
            and parent.get("typed_ring_expression_supported") is False)

def _express_selected(selection: PrincipalParentSelection, info: dict) -> list[dict]:
    parents = []
    for skeleton in selection.skeletons.candidates:
        parent = ((express_ring_principal(info, selection.principal, skeleton))
                  if skeleton.topology is SkeletonTopology.RING_SYSTEM
                  else express_chain_principal(info, selection.principal, skeleton))
        if parent is not None and not _unsupported_typed_ring(parent, selection):
            parents.append(parent)
    return parents

def _retained_ketone_skeleton(selection: PrincipalParentSelection, info: dict) -> bool:
    if selection.principal.group_class is not FunctionalGroupClass.KETONE:
        return False
    from namepredict.layer2.ring_scaffold import _producer_id
    anchors = {a for o in selection.principal.occurrences for a in o.parent_anchors}
    generic = {"cycloalkane", "cycloketone", "cycloalkanedione"}
    return any(anchors <= set(s.atom_ids) and (sid := _producer_id(info, s)) is not None
               and sid not in generic for s in selection.skeletons.candidates)

def _needs_special(selection: PrincipalParentSelection, parents: list[dict], info: dict) -> bool:
    if _retained_ketone_skeleton(selection, info):
        return False
    if selection.principal.group_class is not FunctionalGroupClass.KETONE:
        return True
    facts = [p.get("principal_expression_facts") for p in parents]
    return not parents or any(f and f.relation.value == "exocyclic" for f in facts)

def rule_driven_parent_candidates(info: dict) -> list[dict]:
    selection = select_principal_parent_skeletons(info)
    if selection.skeletons is None:
        return []
    if selection.principal is None:
        from namepredict.layer2.principal_expression import express_hydrocarbon_principal
        return [p for s in selection.skeletons.candidates
                if (p := express_hydrocarbon_principal(info, s)) is not None]
    parents = _express_selected(selection, info)
    special = _owned(_special_expression(selection, info), selection) if _needs_special(selection, parents, info) else None
    extra = [special] if special else []
    typed_first = selection.principal.group_class in {
        FunctionalGroupClass.KETONE, FunctionalGroupClass.AMINE,
    }
    return parents + extra if typed_first and parents else extra + parents
