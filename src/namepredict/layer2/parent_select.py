"""Layer2 母体选择（P-44）：主官能团骨架表达、候选收集与并列最优返回。"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem import Mol

from namepredict.layer1.functional_group_inventory import FunctionalGroupClass, inventory_from_info
from namepredict.layer2.parent_skeleton import (
    SkeletonSelection,
    SkeletonTopology,
    select_principal_skeletons,
)
from namepredict.layer2.principal_expression import express_chain_principal, express_ring_principal
from namepredict.layer2.principal import PrincipalGroupSelection, select_principal_group
from namepredict.constants import OXO_CENTER_KINDS


@dataclass(frozen=True)
class PrincipalParentSelection:
    """L2 选择结果：主官能团选择与骨架选择（各自可能为 None）。"""
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


def _chain_atoms(parent: dict) -> set[int]:
    """取母体链原子集合。"""
    return set(parent.get("chain") or ())


_WHOLE_FG_ATOMS = frozenset({FunctionalGroupClass.OXOACID, FunctionalGroupClass.SULFONAMIDE})  # 特征原子全归主基团：中心非骨架成员，无可外借臂


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    """主官能团所有权原子：骨架内（或邻骨架）锚点及其特征原子。"""
    facts = parent.get("principal_expression_facts")
    occurrences = parent.get("principal_occurrences") or ()
    if facts is None:
        return set()
    chain = _chain_atoms(parent)
    anchors = {a for o in occurrences for a in o.parent_anchors}
    atoms = {a for o in occurrences for a in o.characteristic_atoms}
    seeds = anchors & chain
    if not seeds:  # 锚点全在骨架外：exocyclic 基团（苯甲酸的羧基），改取与骨架相邻的锚点
        linked = {n.GetIdx() for i in chain for n in mol.GetAtomWithIdx(i).GetNeighbors()}
        seeds = anchors & linked
    out = set(seeds)
    if facts.group_class in _WHOLE_FG_ATOMS:  # 含氧酸：中心的氧/卤素跨两跳，仅归本骨架覆盖的 occurrence
        covered = set(parent.get("covered_principal_ids") or ())
        return out | {i for o in occurrences if o.id in covered for i in o.characteristic_atoms}
    for i in tuple(out):
        out |= {n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
                if n.GetIdx() in atoms and n.GetIdx() not in anchors}
    return out


def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
    """一次性复制候选，生成不可变 owned_atoms frozenset。"""
    if isinstance(parent.get("owned_atoms"), frozenset):
        return parent
    fg_atoms = _kind_fg_atoms(parent, mol)
    if parent.get("kind") in OXO_CENTER_KINDS:  # 含氧酸中心母体：链仅供编号，臂一律退为取代基
        return {**parent, "owned_atoms": frozenset(fg_atoms)}
    owned = frozenset(_chain_atoms(parent) | fg_atoms)  # 链与主官能团特征原子的并集
    return {**parent, "owned_atoms": owned}


def _collect_candidates(info: dict) -> list[dict]:
    """收集并去重规则驱动的母体候选（剔除 None）。"""
    raw = rule_driven_parent_candidates(info)
    return [c for c in raw if c is not None]


def _p45_2_prefix_count(info: dict, parent: dict) -> int:
    """P-45.2.1 键：owned_atoms 边界外的 claim 个数。"""
    mol = info["mol"]
    from namepredict.layer3.claimable_block import iter_claims
    return len(iter_claims(mol, parent.get("owned_atoms") or frozenset()))


def _condensed_rank(info: dict, cand: dict) -> int:
    """候选所辖缩合含氧酸中心的链内桥氧数（P-67.2.1：多核磷酸/硫酸以链中中心为功能母体）；其余恒 0。"""
    mol = info.get("mol")
    if mol is None or cand.get("oxo_kind") not in ("phosphate", "sulfate"):
        return 0
    from namepredict.layer1.analyzer import _oxo_bridge_arms
    ids = set(cand.get("covered_principal_ids") or ())
    zs = [o.payload["oxo_z"] for o in cand.get("principal_occurrences") or ()
          if o.id in ids and o.payload.get("oxo_z") is not None]
    return max((_oxo_bridge_arms(mol, int(z)) for z in zs), default=0)


def _reorder_p45_2(info: dict, cands: list[dict], *, tied: bool = False) -> list[dict]:
    """P-45.2 流水线：按前缀取代基团数目（P-45.2.1）稳定重排候选。"""
    if len(cands) <= 1:
        return cands
    keyed = sorted(((_p45_2_prefix_count(info, c), _condensed_rank(info, c), i, c)
                    for i, c in enumerate(cands)), key=lambda t: (-t[0], -t[1], t[2]))
    if not tied:
        return [c for _, _, _, c in keyed]
    top = keyed[0][0]
    return [c for k, _, _, c in keyed if k == top]


def _finalize_ranked(info: dict, cands: list[dict]) -> list[dict]:
    """补齐契约/词干/编号并固化 owned_atoms。"""
    from namepredict.layer2.kind_registry import pack_parent_stem
    mol = info.get("mol")
    return [
        finalize_parent_ownership(
            pack_parent_stem(c, mol), mol,
        )
        for c in cands
    ]


def select_parent(info: dict) -> list[dict]:
    """P-44 收集->P-45.2->P-45.2.1 并列最优"""
    cands = _finalize_ranked(info, _collect_candidates(info))
    return _reorder_p45_2(info, cands, tied=True)
