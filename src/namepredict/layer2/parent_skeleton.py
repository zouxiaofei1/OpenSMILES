"""P-44 的拓扑优先母体骨架枚举。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.constants import Al, As, B, Bi, C, Ga, Ge, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl
from namepredict.layer1.functional_group_inventory import FunctionalGroupOccurrence
from namepredict.layer2.chain_walk import _chain_through, _chain_through_two, _longest_chain


class SkeletonTopology(str, Enum):
    ACYCLIC = "acyclic"
    RING_SYSTEM = "ring_system"


@dataclass(frozen=True)
class ParentSkeleton:
    topology: SkeletonTopology
    atom_ids: tuple[int, ...]
    covered_principal_ids: frozenset[str]
    scaffold_id: str | None = None


@dataclass(frozen=True)
class PrincipalSkeletons:
    candidates: tuple[ParentSkeleton, ...]
    unsupported_ids: frozenset[str]


@dataclass(frozen=True)
class SkeletonSelection:
    candidates: tuple[ParentSkeleton, ...]
    next_rule: str | None
    unsupported_ids: frozenset[str] = frozenset()


def _anchors(occurrences: tuple[FunctionalGroupOccurrence, ...]) -> list[int]:
    return sorted({a for occurrence in occurrences for a in occurrence.parent_anchors})


def _pair_chains(mol: Mol, anchors: list[int]) -> list[list[int]]:
    pairs = [(a, b) for i, a in enumerate(anchors) for b in anchors[i + 1:]]
    paths = [(a, b, _chain_through_two(mol, a, b)) for a, b in pairs]
    return [path for a, b, path in paths if a in path and b in path]


def _open_chains(mol: Mol, anchors: list[int]) -> list[list[int]]:
    open_anchors = [a for a in anchors if not mol.GetAtomWithIdx(a).IsInRing()]
    singles = [_chain_through({"mol": mol}, anchor) for anchor in open_anchors]
    out = singles + _pair_chains(mol, open_anchors)
    if out:
        return out
    # 无主官能团（纯烃）：最长链作为唯一开链骨架候选。
    chain = _longest_chain(mol)
    return [chain] if chain else []


def _chain_coverage(chain: list[int], occurrences) -> frozenset[str]:
    atoms = set(chain)
    return frozenset(o.id for o in occurrences if o.parent_anchors and o.parent_anchors <= atoms)


def _ring_attaches(mol: Mol, ring: set[int], occurrence: FunctionalGroupOccurrence) -> bool:
    if occurrence.parent_anchors & ring:
        return True
    return any(n.GetIdx() in ring for a in occurrence.parent_anchors for n in mol.GetAtomWithIdx(a).GetNeighbors())


def _ring_candidate(mol: Mol, system: dict, occurrences) -> ParentSkeleton | None:
    atoms = set(system.get("atom_ids") or ())
    covered = frozenset(o.id for o in occurrences if _ring_attaches(mol, atoms, o))
    # 有主官能团时要求环附着至少一个 occurrence；纯烃（无 occurrence）则全部枚举。
    if occurrences and not covered:
        return None
    return ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(sorted(atoms)), covered)


def _ring_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    # scaffold 身份 (scaffold_id) 不参与骨架选择，只对最终胜出的少数骨架有意义，延迟到表达阶段 resolve_ring_scaffold 再识别；此处不跑 producer，避免为每个环系统支付完整 parent 生成器成本。
    mol = info["mol"]
    basic = [c for system in info.get("ring_systems") or () if (c := _ring_candidate(mol, system, occurrences))]
    return [ParentSkeleton(c.topology, c.atom_ids, c.covered_principal_ids) for c in basic]


def _chain_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    paths = _open_chains(info["mol"], _anchors(occurrences))
    unique = {frozenset(path): path for path in paths if path}
    return [ParentSkeleton(SkeletonTopology.ACYCLIC, tuple(path), _chain_coverage(path, occurrences)) for path in unique.values()]


_SENIOR_ATOMS = (N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl, O, S, Se, Te, C)


def _senior_atom(mol: Mol, skeleton: ParentSkeleton) -> int:
    present = {mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids}
    return next((z for z in _SENIOR_ATOMS if z in present), 0)


def keep_senior_atom(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    seniority = {z: i for i, z in enumerate(reversed(_SENIOR_ATOMS), 1)}
    best = max((seniority.get(_senior_atom(mol, c), 0) for c in candidates), default=0)
    return tuple(c for c in candidates if seniority.get(_senior_atom(mol, c), 0) == best)


def keep_ring_over_chain(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    topologies = {c.topology for c in candidates}
    return tuple(c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM) if len(topologies) > 1 else candidates


def keep_p44_1_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    topologies = {c.topology for c in candidates}
    return keep_ring_over_chain(keep_senior_atom(mol, candidates)) if len(topologies) > 1 else candidates


def _hetero_count(mol: Mol, skeleton: ParentSkeleton) -> int:
    return sum(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in skeleton.atom_ids)


def _element_counts(mol: Mol, skeleton: ParentSkeleton) -> tuple[int, ...]:
    numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids]
    return tuple(numbers.count(z) for z in _SENIOR_ATOMS if z != 6)


def p44_3_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    return _hetero_count(mol, skeleton), len(skeleton.atom_ids), _element_counts(mol, skeleton)


def keep_p44_3(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    chains = tuple(c for c in candidates if c.topology is SkeletonTopology.ACYCLIC)
    best = max((p44_3_key(mol, c) for c in chains), default=())
    return tuple(c for c in chains if p44_3_key(mol, c) == best)


def _ring_count(mol: Mol, skeleton: ParentSkeleton) -> int:
    atoms = set(skeleton.atom_ids)
    return sum(set(ring) <= atoms for ring in mol.GetRingInfo().AtomRings())


def p44_2_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids]
    hetero = [z for z in numbers if z != C]
    has_n = N in hetero
    senior = next((len(_SENIOR_ATOMS) - i for i, z in enumerate(_SENIOR_ATOMS) if z in hetero), 0)
    return bool(hetero), numbers.count(N) if has_n else 0, senior if not has_n else 0, _ring_count(mol, skeleton), len(numbers), len(hetero)


def keep_p44_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    rings = tuple(c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM)
    best = max((p44_2_key(mol, c) for c in rings), default=())
    return tuple(c for c in rings if p44_2_key(mol, c) == best)


def keep_max_principal_coverage(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    maximum = max((len(c.covered_principal_ids) for c in candidates), default=0)
    return tuple(c for c in candidates if len(c.covered_principal_ids) == maximum)


def _principal_multiple_edges(mol: Mol, occurrences) -> set[frozenset[int]]:
    edges = set()
    for occurrence in occurrences:
        atoms = occurrence.characteristic_atoms
        edges |= {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()
                  if b.GetBondTypeAsDouble() > 1 and {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= atoms}
    return edges


def p44_4_unsaturation_key(mol: Mol, skeleton: ParentSkeleton, occurrences=()) -> tuple[int, int]:
    atoms, excluded = set(skeleton.atom_ids), _principal_multiple_edges(mol, occurrences)
    bonds = [b for b in mol.GetBonds() if not b.GetIsAromatic()
             and {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= atoms
             and frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) not in excluded]
    return sum(b.GetBondTypeAsDouble() > 1 for b in bonds), sum(b.GetBondTypeAsDouble() == 2 for b in bonds)


def keep_p44_4_unsaturation(mol: Mol, candidates: tuple[ParentSkeleton, ...], occurrences=()) -> tuple[ParentSkeleton, ...]:
    best = max((p44_4_unsaturation_key(mol, c, occurrences) for c in candidates), default=())
    return tuple(c for c in candidates if p44_4_unsaturation_key(mol, c, occurrences) == best)


def _finish_skeleton_selection(info: dict, candidates, occurrences, unsupported) -> SkeletonSelection:
    candidates = keep_p44_4_unsaturation(info["mol"], candidates, occurrences)
    next_rule = "L4:P-44.4.1.3+" if candidates else None
    return SkeletonSelection(candidates, next_rule, unsupported)


def select_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> SkeletonSelection:
    enumerated = enumerate_principal_skeletons(info, occurrences)
    candidates = keep_max_principal_coverage(enumerated.candidates)
    candidates = keep_p44_1_2(info["mol"], candidates)
    topologies = {candidate.topology for candidate in candidates}
    if topologies == {SkeletonTopology.ACYCLIC}:
        candidates = keep_p44_3(info["mol"], candidates)
    else:
        candidates = keep_p44_2(info["mol"], candidates)
    return _finish_skeleton_selection(info, candidates, occurrences, enumerated.unsupported_ids)


def enumerate_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> PrincipalSkeletons:
    candidates = tuple(_chain_candidates(info, occurrences) + _ring_candidates(info, occurrences))
    covered = frozenset(i for candidate in candidates for i in candidate.covered_principal_ids)
    return PrincipalSkeletons(candidates, frozenset(o.id for o in occurrences) - covered)
