"""Resolve ring skeletons to topology-only scaffold identities."""
from __future__ import annotations

from namepredict.layer2.parent_skeleton import ParentSkeleton
from namepredict.layer2.scaffold.retained_registry import match_systems
from namepredict.layer2.scaffold.identity import ScaffoldIdentity
from namepredict.layer2.scaffold.specs import get_identity


def _matched_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    atoms = set(skeleton.atom_ids)
    return next((sid for sid, system, _ in match_systems(info)
                 if set(system.get("atom_ids") or ()) == atoms), None)


def _producer_scaffold_ids(info: dict) -> dict[frozenset[int], str]:
    """Mother-ring scaffold identification via the core table, memoized per info.

    `ring_core` matches mother rings WITHOUT the substituent/outside/FG gates,
    so substituted mothers (e.g. 2,3,6-trimethylquinoline) are still tagged.
    The full ring producers are NOT consulted here: a fallback over
    `ring_try_fns` never hit on real data and was removed.
    """
    cache = info.get("_scaffold_ids")
    if cache is not None:
        return cache
    from namepredict.layer2.scaffold.ring_core import ring_core_fns
    found: dict[frozenset[int], str] = {}
    for fn in ring_core_fns():
        for atoms, sid in fn(info) or ():
            found.setdefault(frozenset(atoms), sid)
    info["_scaffold_ids"] = found
    return found


def _producer_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    # 开链骨架不是环母体：直接短路。
    from namepredict.layer2.parent_skeleton import SkeletonTopology

    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    # 只查 core 表缓存。完整的 producer 回退曾是死代码（真实数据上从未命中，
    # 且会抢占羰基母环的专用 core 适配器），已删除。
    return _producer_scaffold_ids(info).get(frozenset(skeleton.atom_ids))


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
