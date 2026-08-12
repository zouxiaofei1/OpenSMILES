"""Resolve ring skeletons to topology-only scaffold identities."""
from __future__ import annotations

from namepredict.layer2.parent_skeleton import ParentSkeleton
from namepredict.layer2.retained_registry import match_systems
from namepredict.layer2.identity import ScaffoldIdentity
from namepredict.layer2.specs import get_identity


def _matched_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    atoms = set(skeleton.atom_ids)
    return next((sid for sid, system, _ in match_systems(info)
                 if set(system.get("atom_ids") or ()) == atoms), None)





def _producer_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    # 开链骨架不是环母体：直接短路。
    from namepredict.layer2.parent_skeleton import SkeletonTopology

    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    # 只查 core 表缓存。完整的 producer 回退曾是死代码（真实数据上从未命中，
    # 且会抢占羰基母环的专用 core 适配器），已删除。
    return None


def _generic_carbocycle(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids):
        return None
    return ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")


def resolve_ring_scaffold(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    direct = get_identity(skeleton.scaffold_id or "")
    if direct:
        return direct
    sid = _matched_id(info, skeleton)
    if sid:
        return get_identity(sid)
    sid = _producer_id(info, skeleton)
    return get_identity(sid) if sid else _generic_carbocycle(info, skeleton)
