"""Topology-first parent skeleton enumeration for P-44."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.constants import Al, As, B, Bi, C, Ga, Ge, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl
from namepredict.layer1.functional_group_inventory import FunctionalGroupOccurrence
from namepredict.layer2.chain_walk import _chain_through, _chain_through_two


class SkeletonTopology(str, Enum):
    ACYCLIC = "acyclic"
    RING_SYSTEM = "ring_system"


@dataclass(frozen=True)
class ParentSkeleton:
    topology: SkeletonTopology
    atom_ids: tuple[int, ...]
    covered_principal_ids: frozenset[str]


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
    return singles + _pair_chains(mol, open_anchors)


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
    return ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(sorted(atoms)), covered) if covered else None


def _ring_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    mol = info["mol"]
    return [c for system in info.get("ring_systems") or () if (c := _ring_candidate(mol, system, occurrences))]


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


def keep_max_ring_system_size(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    rings = tuple(c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM)
    if not rings:
        return candidates
    maximum = max(len(c.atom_ids) for c in rings)
    return tuple(c for c in rings if len(c.atom_ids) == maximum)


def keep_max_principal_coverage(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    maximum = max((len(c.covered_principal_ids) for c in candidates), default=0)
    return tuple(c for c in candidates if len(c.covered_principal_ids) == maximum)


def select_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> SkeletonSelection:
    enumerated = enumerate_principal_skeletons(info, occurrences)
    candidates = keep_max_principal_coverage(enumerated.candidates)
    candidates = keep_p44_1_2(info["mol"], candidates)
    topologies = {candidate.topology for candidate in candidates}
    if topologies == {SkeletonTopology.ACYCLIC}:
        return SkeletonSelection(keep_p44_3(info["mol"], candidates), "P-44.4", enumerated.unsupported_ids)
    candidates = keep_max_ring_system_size(candidates)
    return SkeletonSelection(candidates, "P-44.2-complete" if candidates else None, enumerated.unsupported_ids)


def enumerate_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> PrincipalSkeletons:
    candidates = tuple(_chain_candidates(info, occurrences) + _ring_candidates(info, occurrences))
    covered = frozenset(i for candidate in candidates for i in candidate.covered_principal_ids)
    return PrincipalSkeletons(candidates, frozenset(o.id for o in occurrences) - covered)
