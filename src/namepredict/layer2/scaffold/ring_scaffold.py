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
    FG-variant retained names (benzofuranamine, ...) are NOT covered here —
    `_producer_id` falls back to the full producers lazily for those skeletons.
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
    # 开链骨架不是环母体：直接短路，避免为每个开链候选跑完整 producer 兜底。
    from namepredict.layer2.parent_skeleton import SkeletonTopology

    if skeleton.topology is not SkeletonTopology.RING_SYSTEM:
        return None
    key = frozenset(skeleton.atom_ids)
    if key in (found := _producer_scaffold_ids(info)):
        return found[key]
    # Core table missed this skeleton (FG-variant retained name): lazily ask the
    # full producers for exactly this ring-atom set.
    from namepredict.layer2 import kind_registry as registry
    for producer in registry.ring_try_fns():
        parent = producer(info)
        atoms = frozenset(parent.get("chain") or ()) if parent else frozenset()
        if atoms == key:
            return parent.get("scaffold_id") or parent.get("kind")
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
