"""P-44 的拓扑优先母体骨架枚举。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.constants import Al, As, B, Bi, C, Ga, Ge, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl
from namepredict.layer1.functional_group_inventory import (
    FunctionalGroupClass, FunctionalGroupOccurrence, inventory_from_info,
)
from namepredict.layer2.chain_walk import _all_chains_through, _chain_through_two, _longest_chain
from namepredict.layer1.ring_systems import sssr_rings
from namepredict.layer4.numbering_engine import narrow


class SkeletonTopology(str, Enum):
    """母体骨架拓扑类型（开链/环系）。"""
    ACYCLIC = "acyclic"
    RING_SYSTEM = "ring_system"


@dataclass(frozen=True)
class ParentSkeleton:
    """一个母体骨架候选：拓扑、原子集与覆盖的主基团 id。"""
    topology: SkeletonTopology
    atom_ids: tuple[int, ...]
    covered_principal_ids: frozenset[str]


@dataclass(frozen=True)
class SkeletonSelection:
    """骨架候选集合：枚举结果或筛选后的胜出候选。"""
    candidates: tuple[ParentSkeleton, ...]


def _anchors(occurrences: tuple[FunctionalGroupOccurrence, ...]) -> list[int]:
    """收集所有 occurrence 的主官能团锚点碳（去重排序）。"""
    return sorted({a for occurrence in occurrences for a in occurrence.parent_anchors})


def _demoted_leaf_carbons(info: dict) -> set[int]:
    """被压制为前缀叶的碳：demoted 条目中心碳不进主链（P-61.1.3）。"""
    mol: Mol = info["mol"]
    return {int(o.payload["center_idx"]) for o in inventory_from_info(info).demoted_entries()
            if o.payload.get("center_idx") is not None
            and mol.GetAtomWithIdx(int(o.payload["center_idx"])).GetAtomicNum() == C}


def _pair_chains(mol: Mol, anchors: list[int], banned: set[int] = frozenset()) -> list[list[int]]:
    """生成连接两锚点的链（穿过两点的最长链）。"""
    pairs = [(a, b) for i, a in enumerate(anchors) for b in anchors[i + 1:]]
    paths = [(a, b, _chain_through_two(mol, a, b, banned)) for a, b in pairs]
    return [path for a, b, path in paths if a in path and b in path]


def _open_chains(mol: Mol, anchors: list[int], banned: set[int] = frozenset()) -> list[list[int]]:
    """开环锚点的单链与两两链候选（无锚点时退化为最长链）。"""
    open_anchors = [a for a in anchors if not mol.GetAtomWithIdx(a).IsInRing()]
    singles = [chain for anchor in open_anchors for chain in _all_chains_through(mol, anchor, banned)]  # 穿过锚点的等长最长链全枚举（平局交 P-44.4/45.2 裁决）
    out = singles + _pair_chains(mol, open_anchors, banned)
    if out:
        return out
    chain = _longest_chain(mol, banned=banned)  # 无主官能团（纯烃）：最长链作为唯一开链骨架候选。
    return [chain] if chain else []


def _chain_coverage(chain: list[int], occurrences) -> frozenset[str]:
    """链覆盖的 occurrence id 集合（胺任一臂在链即覆盖，P-62.2）。"""
    atoms = set(chain)
    out: set[str] = set()
    for o in occurrences:
        if not o.parent_anchors:
            continue
        if o.group_class is FunctionalGroupClass.AMINE:
            if o.parent_anchors & atoms:
                out.add(o.id)
        elif o.parent_anchors <= atoms:
            out.add(o.id)
    return frozenset(out)


def _ring_attaches(mol: Mol, ring: set[int], occurrence: FunctionalGroupOccurrence) -> bool:
    """判断 occurrence 是否附着于环（胺/醇/酮/自由基只认直接附着）。"""
    if occurrence.parent_anchors & ring:
        return True
    if occurrence.group_class in (FunctionalGroupClass.AMINE, FunctionalGroupClass.ALCOHOL,FunctionalGroupClass.RADICAL,FunctionalGroupClass.KETONE):
        return False
    return any(n.GetIdx() in ring for a in occurrence.parent_anchors for n in mol.GetAtomWithIdx(a).GetNeighbors())


def _ring_candidate(mol: Mol, system: dict, occurrences) -> ParentSkeleton | None:
    """构造环骨架候选（有 occurrence 时要求附着）。"""
    atoms = set(system.get("atom_ids") or ())
    covered = frozenset(o.id for o in occurrences if _ring_attaches(mol, atoms, o))
  
    return ParentSkeleton(SkeletonTopology.RING_SYSTEM, tuple(sorted(atoms)), covered)


def _ring_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    """枚举全部环系统的骨架候选。"""
    mol = info["mol"]  # scaffold 身份延迟到表达阶段识别，此处不跑 producer
    basic = [c for system in info.get("ring_systems") or () if (c := _ring_candidate(mol, system, occurrences))]
    return basic


def _cation_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    """单核母体阳离子：只以阳离子原子本身作骨架候选（P-73.1.1）。"""
    return [ParentSkeleton(SkeletonTopology.ACYCLIC, (a,),
                           frozenset(o.id for o in occurrences if a in o.parent_anchors))
            for a in _anchors(occurrences)]


def _chain_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    """枚举去重后的开链骨架候选。"""
    paths = _open_chains(info["mol"], _anchors(occurrences), _demoted_leaf_carbons(info))
    unique = {frozenset(path): path for path in paths if path}
    return [ParentSkeleton(SkeletonTopology.ACYCLIC, tuple(path), _chain_coverage(path, occurrences)) for path in unique.values()]


_SENIOR_ATOMS = (N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl, O, S, Se, Te, C)
_SENIORITY = {z: i for i, z in enumerate(reversed(_SENIOR_ATOMS), 1)}  # 元素 → 优先序数（N 最高，C 最低）


def _senior_atom(mol: Mol, skeleton: ParentSkeleton) -> int:
    """骨架中存在的最优先元素原子序数。"""
    present = {mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids}
    return next((z for z in _SENIOR_ATOMS if z in present), 0)


def keep_senior_atom(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """保留含最优先元素（senior）的候选。"""
    return tuple(narrow(list(candidates), lambda c: _SENIORITY.get(_senior_atom(mol, c), 0),
                        reverse=True))


def keep_p44_1_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.1.2：混合拓扑时取 senior 元素。"""
    topologies = {c.topology for c in candidates}
    return keep_senior_atom(mol, candidates) if len(topologies) > 1 else candidates

def p44_3_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    """P-44.3 比较键：杂原子数、原子数、元素计数。"""
    return (sum(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in skeleton.atom_ids),  # 杂原子数
            len(skeleton.atom_ids), None)


def keep_p44_3(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.3：开链候选按 p44_3_key 取最大。"""
    chains = [c for c in candidates if c.topology is SkeletonTopology.ACYCLIC]
    return tuple(narrow(chains, lambda c: p44_3_key(mol, c), reverse=True))


def _ring_count(mol: Mol, skeleton: ParentSkeleton) -> int:
    """计算骨架覆盖的完整环数。"""
    atoms = set(skeleton.atom_ids)
    return sum(set(ring) <= atoms for ring in sssr_rings(mol))


def p44_2_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    """P-44.2 比较键：杂原子、N 计数、senior、环数等。"""
    numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids]
    hetero = [z for z in numbers if z != C]
    has_n = N in hetero
    senior = next((_SENIORITY[z] for z in _SENIOR_ATOMS if z in hetero), 0)
    return bool(hetero), numbers.count(N) if has_n else 0, senior if not has_n else 0, _ring_count(mol, skeleton), len(numbers), len(hetero)


def keep_p44_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.2：环候选按 p44_2_key 取最大。"""
    rings = [c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM]
    return tuple(narrow(rings, lambda c: p44_2_key(mol, c), reverse=True))


def keep_max_principal_coverage(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """保留覆盖主官能团最多的候选。"""
    return tuple(narrow(list(candidates), lambda c: len(c.covered_principal_ids), reverse=True))

def p44_4_unsaturation_key(mol: Mol, skeleton: ParentSkeleton, occurrences=()) -> tuple[int, int]:
    """P-44.4 不饱和度键：(多重键数, 双键数)。 """
    return memo.by_key(
        "p44_4", (id(mol), tuple(skeleton.atom_ids), tuple(id(o) for o in occurrences)),
        lambda: _p44_4_unsaturation_key_uncached(mol, skeleton, occurrences),
        mol, *occurrences,  # 键里用了 id()，须保活防复用
    )


def _p44_4_unsaturation_key_uncached(mol: Mol, skeleton: ParentSkeleton, occurrences=()) -> tuple[int, int]:
    """实际计算 P-44.4 不饱和度键（无记忆版本）。"""
    atoms = set(skeleton.atom_ids)
    non_arom: list = []
    n_arom = 0
    for b in mol.GetBonds():
        edge = frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        if edge <= atoms :
            if b.GetIsAromatic():
                n_arom += 1
            else:
                non_arom.append(b)
    n_multi = sum(b.GetBondTypeAsDouble() > 1 for b in non_arom) + n_arom // 2
    n_db = sum(b.GetBondTypeAsDouble() == 2 for b in non_arom) + n_arom // 2
    return (n_multi, n_db)


def keep_p44_4_unsaturation(mol: Mol, candidates: tuple[ParentSkeleton, ...], occurrences=()) -> tuple[ParentSkeleton, ...]:
    """P-44.4：按不饱和度键取最大候选。"""
    return tuple(narrow(list(candidates),
                        lambda c: p44_4_unsaturation_key(mol, c, occurrences), reverse=True))


def select_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> SkeletonSelection:
    """按 P-44 顺序筛选骨架（覆盖度→环/链→P-44.3/2）。"""
    enumerated = enumerate_principal_skeletons(info, occurrences)
    candidates = keep_max_principal_coverage(enumerated.candidates)
    candidates = keep_p44_1_2(info["mol"], candidates)
    topologies = {candidate.topology for candidate in candidates}
    if topologies == {SkeletonTopology.ACYCLIC}:
        candidates = keep_p44_3(info["mol"], candidates)
    else:
        candidates = keep_p44_2(info["mol"], candidates)
    return SkeletonSelection(keep_p44_4_unsaturation(info["mol"], candidates, occurrences))  # 末位应用不饱和度规则


def enumerate_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> SkeletonSelection:
    """枚举全部骨架候选；主基团为单核阳离子时只向阳离子原子收敛，不进链/环枚举。"""
    if occurrences and all(o.group_class is FunctionalGroupClass.CATION for o in occurrences):
        # P-73.7(c)：多阳离子中心时取优先元素（N > P > … > O > S），比 P-44 拓扑规则更专
        return SkeletonSelection(tuple(keep_senior_atom(info["mol"], tuple(_cation_candidates(info, occurrences)))))
    return SkeletonSelection(tuple(_chain_candidates(info, occurrences) + _ring_candidates(info, occurrences)))
