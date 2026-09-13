"""P-44 的拓扑优先母体骨架枚举。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.constants import Al, As, B, Bi, C, Ga, Ge, In, N, O, P, Pb, S, Sb, Se, Si, Sn, Te, Tl
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass, FunctionalGroupOccurrence
from namepredict.layer2.chain_walk import _all_chains_through, _chain_through_two, _longest_chain
from namepredict.layer1.ring_systems import sssr_rings


class SkeletonTopology(str, Enum):
    """母体骨架拓扑类型（开链/环系）。"""
    ACYCLIC = "acyclic"
    RING_SYSTEM = "ring_system"


@dataclass(frozen=True)
class ParentSkeleton:
    """一个母体骨架候选：拓扑、原子集、覆盖的主基团 id 与 scaffold 身份。"""
    topology: SkeletonTopology
    atom_ids: tuple[int, ...]
    covered_principal_ids: frozenset[str]
    scaffold_id: str | None = None


@dataclass(frozen=True)
class PrincipalSkeletons:
    """骨架枚举结果：全部候选与未被任何候选覆盖的 occurrence id。"""
    candidates: tuple[ParentSkeleton, ...]
    unsupported_ids: frozenset[str]


@dataclass(frozen=True)
class SkeletonSelection:
    """骨架筛选结果：胜出候选、交 L4 的下一规则标记与未覆盖 id。"""
    candidates: tuple[ParentSkeleton, ...]
    next_rule: str | None
    unsupported_ids: frozenset[str] = frozenset()


def _anchors(occurrences: tuple[FunctionalGroupOccurrence, ...]) -> list[int]:
    """收集所有 occurrence 的主官能团锚点碳（去重排序）。"""
    return sorted({a for occurrence in occurrences for a in occurrence.parent_anchors})


def _demoted_leaf_carbons(info: dict) -> set[int]:
    """被压制（降级）为前缀叶的碳集合：`demoted_*` 叶列表的中心碳（羧酸 carboxy 叶 P-61.1.3、腈 cyano 叶等）不得进入开链主链——否则词干链会把叶碳当饱和碳吞掉、其杂原子悬空误命名成 hydroxy/amino。"""
    mol: Mol = info["mol"]
    out: set[int] = set()
    for key, entries in info.items():
        if not key.startswith("demoted_"):
            continue
        for e in entries or []:
            if e.get("center_idx") is not None and mol.GetAtomWithIdx(int(e["center_idx"])).GetAtomicNum() == C:
                out.add(int(e["center_idx"]))
    return out


def _pair_chains(mol: Mol, anchors: list[int], banned: set[int] = frozenset()) -> list[list[int]]:
    """生成连接两锚点的链（穿过两点的最长链）。"""
    pairs = [(a, b) for i, a in enumerate(anchors) for b in anchors[i + 1:]]
    paths = [(a, b, _chain_through_two(mol, a, b, banned)) for a, b in pairs]
    return [path for a, b, path in paths if a in path and b in path]


def _open_chains(mol: Mol, anchors: list[int], banned: set[int] = frozenset()) -> list[list[int]]:
    """开环锚点的单链与两两链候选（无锚点时退化为最长链）。"""
    open_anchors = [a for a in anchors if not mol.GetAtomWithIdx(a).IsInRing()]
    singles = [chain for anchor in open_anchors for chain in _all_chains_through(mol, anchor, banned)]  # 穿过锚点的等长最长链全部枚举（平局候选让 P-44.4/P-45.2 裁决，如醛端连甲基 vs 羟甲基）。
    out = singles + _pair_chains(mol, open_anchors, banned)
    if out:
        return out
    chain = _longest_chain(mol, banned=banned)  # 无主官能团（纯烃）：最长链作为唯一开链骨架候选。
    return [chain] if chain else []


def _chain_coverage(chain: list[int], occurrences) -> frozenset[str]:
    """计算链覆盖的 occurrence id 集合：锚点全在链上（二级/三级胺 N 任一臂在链即覆盖，P-62.2 多臂选优，余臂作 N- 取代基留在所有权外）。"""
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
    """判断 occurrence 是否附着于环（锚点本身/邻居在环内；胺/醇/酮/自由基仅认锚点直接附着），避免苄基胺、苄醇等隔碳连芳环错选环母体。"""
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
    mol = info["mol"]  # scaffold 身份 (scaffold_id) 不参与骨架选择，只对最终胜出的少数骨架有意义，延迟到表达阶段 resolve_ring_scaffold 再识别；此处不跑 producer，避免为每个环系统支付完整 parent 生成器成本。
    basic = [c for system in info.get("ring_systems") or () if (c := _ring_candidate(mol, system, occurrences))]
    return [ParentSkeleton(c.topology, c.atom_ids, c.covered_principal_ids) for c in basic]


def _chain_candidates(info: dict, occurrences) -> list[ParentSkeleton]:
    """枚举去重后的开链骨架候选。"""
    paths = _open_chains(info["mol"], _anchors(occurrences), _demoted_leaf_carbons(info))
    unique = {frozenset(path): path for path in paths if path}
    return [ParentSkeleton(SkeletonTopology.ACYCLIC, tuple(path), _chain_coverage(path, occurrences)) for path in unique.values()]


_SENIOR_ATOMS = (N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl, O, S, Se, Te, C)


def _senior_atom(mol: Mol, skeleton: ParentSkeleton) -> int:
    """骨架中存在的最优先元素原子序数。"""
    present = {mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids}
    return next((z for z in _SENIOR_ATOMS if z in present), 0)


def keep_senior_atom(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """保留含最优先元素（senior）的候选。"""
    seniority = {z: i for i, z in enumerate(reversed(_SENIOR_ATOMS), 1)}
    best = max((seniority.get(_senior_atom(mol, c), 0) for c in candidates), default=0)
    return tuple(c for c in candidates if seniority.get(_senior_atom(mol, c), 0) == best)


def keep_ring_over_chain(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """混合拓扑时保留环骨架。"""
    topologies = {c.topology for c in candidates}
    return tuple(c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM) if len(topologies) > 1 else candidates


def keep_p44_1_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.1.2：混合拓扑时环优先且取 senior 元素。"""
    topologies = {c.topology for c in candidates}
    return keep_ring_over_chain(keep_senior_atom(mol, candidates)) if len(topologies) > 1 else candidates


def _hetero_count(mol: Mol, skeleton: ParentSkeleton) -> int:
    """统计骨架中非碳（杂）原子数。"""
    return sum(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in skeleton.atom_ids)


def _element_counts(mol: Mol, skeleton: ParentSkeleton) -> tuple[int, ...]:
    """统计骨架中各 senior 元素（除碳）的出现次数元组。"""
    numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids]
    return tuple(numbers.count(z) for z in _SENIOR_ATOMS if z != 6)


def p44_3_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    """P-44.3 比较键：杂原子数、原子数、元素计数。"""
    return _hetero_count(mol, skeleton), len(skeleton.atom_ids), _element_counts(mol, skeleton)


def keep_p44_3(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.3：开链候选按 p44_3_key 取最大。"""
    chains = tuple(c for c in candidates if c.topology is SkeletonTopology.ACYCLIC)
    best = max((p44_3_key(mol, c) for c in chains), default=())
    return tuple(c for c in chains if p44_3_key(mol, c) == best)


def _ring_count(mol: Mol, skeleton: ParentSkeleton) -> int:
    """计算骨架覆盖的完整环数。"""
    atoms = set(skeleton.atom_ids)
    return sum(set(ring) <= atoms for ring in sssr_rings(mol))


def p44_2_key(mol: Mol, skeleton: ParentSkeleton) -> tuple:
    """P-44.2 比较键：杂原子、N 计数、senior、环数等。"""
    numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in skeleton.atom_ids]
    hetero = [z for z in numbers if z != C]
    has_n = N in hetero
    senior = next((len(_SENIOR_ATOMS) - i for i, z in enumerate(_SENIOR_ATOMS) if z in hetero), 0)
    return bool(hetero), numbers.count(N) if has_n else 0, senior if not has_n else 0, _ring_count(mol, skeleton), len(numbers), len(hetero)


def keep_p44_2(mol: Mol, candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """P-44.2：环候选按 p44_2_key 取最大。"""
    rings = tuple(c for c in candidates if c.topology is SkeletonTopology.RING_SYSTEM)
    best = max((p44_2_key(mol, c) for c in rings), default=())
    return tuple(c for c in rings if p44_2_key(mol, c) == best)


def keep_max_principal_coverage(candidates: tuple[ParentSkeleton, ...]) -> tuple[ParentSkeleton, ...]:
    """保留覆盖主官能团最多的候选。"""
    maximum = max((len(c.covered_principal_ids) for c in candidates), default=0)
    return tuple(c for c in candidates if len(c.covered_principal_ids) == maximum)


def _principal_multiple_edges(mol: Mol, occurrences) -> set[frozenset[int]]:
    """收集主官能团特征原子间的多重键边（不计入不饱和）。"""
    edges = set()
    for occurrence in occurrences:
        atoms = occurrence.characteristic_atoms
        edges |= {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol.GetBonds()
                  if b.GetBondTypeAsDouble() > 1 and {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= atoms}
    return edges


def p44_4_unsaturation_key(mol: Mol, skeleton: ParentSkeleton, occurrences=()) -> tuple[int, int]:
    """P-44.4 不饱和度键：(多重键数, 双键数)。 """
    return memo.by_key(
        "p44_4", (id(mol), tuple(skeleton.atom_ids), tuple(id(o) for o in occurrences)),
        lambda: _p44_4_unsaturation_key_uncached(mol, skeleton, occurrences),
        mol, *occurrences,  # 键里用了 id()，须保活防复用
    )


def _p44_4_unsaturation_key_uncached(mol: Mol, skeleton: ParentSkeleton, occurrences=()) -> tuple[int, int]:
    """实际计算 P-44.4 不饱和度键（无记忆版本，见 p44_4_unsaturation_key）。"""
    atoms, excluded = set(skeleton.atom_ids), _principal_multiple_edges(mol, occurrences)
    non_arom: list = []
    n_arom = 0
    for b in mol.GetBonds():
        edge = frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        if edge <= atoms and edge not in excluded:
            if b.GetIsAromatic():
                n_arom += 1
            else:
                non_arom.append(b)
    n_multi = sum(b.GetBondTypeAsDouble() > 1 for b in non_arom) + n_arom // 2
    n_db = sum(b.GetBondTypeAsDouble() == 2 for b in non_arom) + n_arom // 2
    return (n_multi, n_db)


def keep_p44_4_unsaturation(mol: Mol, candidates: tuple[ParentSkeleton, ...], occurrences=()) -> tuple[ParentSkeleton, ...]:
    """P-44.4：按不饱和度键取最大候选。"""
    best = max((p44_4_unsaturation_key(mol, c, occurrences) for c in candidates), default=())
    return tuple(c for c in candidates if p44_4_unsaturation_key(mol, c, occurrences) == best)


def _finish_skeleton_selection(info: dict, candidates, occurrences, unsupported) -> SkeletonSelection:
    """末位应用不饱和度规则并封装选择结果。"""
    candidates = keep_p44_4_unsaturation(info["mol"], candidates, occurrences)
    next_rule = "L4:P-44.4.1.3+" if candidates else None
    return SkeletonSelection(candidates, next_rule, unsupported)


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
    return _finish_skeleton_selection(info, candidates, occurrences, enumerated.unsupported_ids)


def enumerate_principal_skeletons(info: dict, occurrences: tuple[FunctionalGroupOccurrence, ...]) -> PrincipalSkeletons:
    """枚举全部骨架候选，并计算未覆盖的 occurrence id。"""
    candidates = tuple(_chain_candidates(info, occurrences) + _ring_candidates(info, occurrences))
    covered = frozenset(i for candidate in candidates for i in candidate.covered_principal_ids)
    return PrincipalSkeletons(candidates, frozenset(o.id for o in occurrences) - covered)
