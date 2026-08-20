"""规则驱动的母体骨架选择入口。"""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass, inventory_from_info
from namepredict.layer2.parent_skeleton import (
    SkeletonSelection,
    SkeletonTopology,
    select_principal_skeletons,
)
from namepredict.layer2.principal_expression import express_chain_principal, express_ring_principal
from namepredict.layer2.principal import PrincipalGroupSelection, select_principal_group

@dataclass(frozen=True)
class PrincipalParentSelection:
    principal: PrincipalGroupSelection | None
    skeletons: SkeletonSelection | None

def select_principal_parent_skeletons(info: dict) -> PrincipalParentSelection:
    """选出主官能团并枚举其骨架选择结果。"""
    principal = select_principal_group(inventory_from_info(info))
    occurrences = principal.occurrences if principal else ()
    skeletons = select_principal_skeletons(info, occurrences)
    return PrincipalParentSelection(principal, skeletons)

def _unsupported_typed_ring(parent: dict, selection: PrincipalParentSelection) -> bool:
    """判断环酮是否因未支持的 typed 表达被排除。"""
    return (selection.principal.group_class is FunctionalGroupClass.KETONE
            and parent.get("scaffold_identity") is not None
            and parent.get("typed_ring_expression_supported") is False)

def _express_selected(selection: PrincipalParentSelection, info: dict) -> list[dict]:
    """按骨架拓扑表达主基团，过滤不支持者。"""
    parents = []
    for skeleton in selection.skeletons.candidates:
        parent = ((express_ring_principal(info, selection.principal, skeleton))
                  if skeleton.topology is SkeletonTopology.RING_SYSTEM
                  else express_chain_principal(info, selection.principal, skeleton))
        if parent is not None and not _unsupported_typed_ring(parent, selection):
            parents.append(parent)
    return parents

def rule_driven_parent_candidates(info: dict) -> list[dict]:
    """规则驱动入口：返回最终母体候选（无主官能团时纯烃）。"""
    
    selection = select_principal_parent_skeletons(info)
    return _express_selected(selection, info)
