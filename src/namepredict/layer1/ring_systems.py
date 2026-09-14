"""由 SSSR 构建环系拓扑：稠合连通、螺环合并。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.constants import C

def _sssr(mol: Mol) -> list[tuple[int, ...]]:
    """返回分子的全部 SSSR 最小环原子序列。"""
    return memo.by_mol("sssr", lambda m: list(m.GetRingInfo().AtomRings()), mol)

def sssr_rings(mol: Mol) -> list[tuple[int, ...]]:
    """供各层统一调用的环访问器（与 _sssr 同一次记忆）。"""
    return _sssr(mol)

def _shared(a: tuple[int, ...], b: tuple[int, ...]) -> frozenset[int]:
    """返回两个环共享的原子集合。"""
    return frozenset(a) & frozenset(b)

def _fusion_edges(rings: list[tuple[int, ...]]) -> list[tuple[int, int, frozenset[int]]]:
    """共享 >=2 个原子的配对边 (i, j, shared_atoms)。"""
    out: list[tuple[int, int, frozenset[int]]] = []
    for i, ri in enumerate(rings):
        for j in range(i + 1, len(rings)):
            sh = _shared(ri, rings[j])
            if len(sh) >= 2:
                out.append((i, j, sh))
    return out

def _spiro_pairs(rings: list[tuple[int, ...]]) -> list[tuple[int, int, int]]:
    """恰好共享 1 个原子的配对：(i, j, atom)。"""
    out: list[tuple[int, int, int]] = []
    for i, ri in enumerate(rings):
        for j in range(i + 1, len(rings)):
            sh = _shared(ri, rings[j])
            if len(sh) == 1:
                out.append((i, j, next(iter(sh))))
    return out

def _uf_find(parent: list[int], x: int) -> int:
    """并查集查找根（含路径压缩）。"""
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def _uf_union(parent: list[int], a: int, b: int) -> None:
    """并查集合并两个根。"""
    ra, rb = _uf_find(parent, a), _uf_find(parent, b)
    if ra != rb:
        parent[rb] = ra

def _components(n: int, edges: list[tuple[int, int, frozenset[int]]]) -> list[list[int]]:
    """对环索引做并查集连通分量。"""
    parent = list(range(n))
    for i, j, _ in edges:
        _uf_union(parent, i, j)
    buckets: dict[int, list[int]] = {}
    for i in range(n):
        buckets.setdefault(_uf_find(parent, i), []).append(i)
    return list(buckets.values())

def _hetero_atoms(mol: Mol, atom_ids: set[int]) -> list[dict]:
    """返回环系中非碳原子（索引与原子序数）列表。"""
    out: list[dict] = []
    for i in sorted(atom_ids):
        z = mol.GetAtomWithIdx(i).GetAtomicNum()
        if z != C:
            out.append({"idx": i, "Z": z})
    return out

def _member_atoms(rings: list[tuple[int, ...]], members: list[int]) -> set[int]:
    """汇总分量内全部环成员的原子集合。"""
    atom_ids: set[int] = set()
    for i in members:
        atom_ids |= set(rings[i])
    return atom_ids

def _member_edges(
    members: list[int], fusion_edges: list[tuple[int, int, frozenset[int]]],
) -> list[tuple[int, int, list[int]]]:
    """筛选仅连接分量内环的稠合边。"""
    mset = set(members)
    return [
        (i, j, sorted(sh)) for i, j, sh in fusion_edges if i in mset and j in mset
    ]

def _system_dict(
    atom_ids: set[int], members: list[int], edges: list, mol: Mol,
) -> dict:
    """组装单个环系的事实 dict。"""
    return {
        "atom_ids": sorted(atom_ids),
        "sssr_indices": sorted(members),
        "fusion_edges": edges,
        "n_rings": len(members),
        "n_atoms": len(atom_ids),
        "hetero_atoms": _hetero_atoms(mol, atom_ids),
        "is_aromatic_mancude": None,
        "topology": None,
    }

def _system_entry(
    mol: Mol,
    rings: list[tuple[int, ...]],
    members: list[int],
    fusion_edges: list[tuple[int, int, frozenset[int]]],
) -> dict:
    """为连通分量构建环系条目。"""
    atoms = _member_atoms(rings, members)
    edges = _member_edges(members, fusion_edges)
    return _system_dict(atoms, members, edges, mol)

def _collect_merged_fields(
    systems: list[dict], indices: list[int],
) -> tuple[set[int], set[int], list[dict]]:
    """收集待合并系统的 sssr、原子与杂原子字段。"""
    sssr: set[int] = set()
    atoms: set[int] = set()
    hetero: list[dict] = []
    for idx in indices:
        s = systems[idx]
        sssr.update(s["sssr_indices"])
        atoms.update(s["atom_ids"])
        hetero.extend(s["hetero_atoms"])
    return sssr, atoms, hetero

def _merged_spiro_system(
    rings: list[tuple[int, ...]],
    sys_indices: list[int], systems: list[dict],
) -> dict:
    """将共享一个螺原子的多个单环系合并为一个螺环系。"""
    sssr, atoms, hetero = _collect_merged_fields(systems, sys_indices)
    return {
        "atom_ids": sorted(atoms),
        "sssr_indices": sorted(sssr),
        "fusion_edges": [],
        "n_rings": len(sssr),
        "n_atoms": len(atoms),
        "hetero_atoms": hetero,
        "is_aromatic_mancude": None,
        "topology": "spiro",
        "ring_sizes": sorted(len(rings[ri]) - 1 for ri in sssr),
    }

def _spiro_sys_indices(
    systems: list[dict], spiro_pairs: list[tuple[int, int, int]],
) -> list[int]:
    """通过螺原子连接的系统并查集父数组。"""
    n = len(systems)
    parent = list(range(n))
    for i, j, _ in spiro_pairs:
        si = _find_sys_for_ring(systems, i)
        sj = _find_sys_for_ring(systems, j)
        if si >= 0 and sj >= 0 and si != sj:
            _uf_union(parent, si, sj)
    return parent

def _find_sys_for_ring(systems: list[dict], ring_idx: int) -> int:
    """返回包含指定环索引的系统下标；无则返回 -1。"""
    for idx, s in enumerate(systems):
        if ring_idx in s["sssr_indices"]:
            return idx
    return -1

def _group_by_root(parent: list[int], n: int) -> dict[int, list[int]]:
    """按并查集根分组下标。"""
    buckets: dict[int, list[int]] = {}
    for idx in range(n):
        buckets.setdefault(_uf_find(parent, idx), []).append(idx)
    return buckets

def _merge_spiro(
    rings: list[tuple[int, ...]],
    systems: list[dict], spiro_pairs: list[tuple[int, int, int]],
) -> list[dict]:
    """将共享螺原子的单环系合并为螺环系。"""
    if not spiro_pairs or len(systems) <= 1:
        return systems
    parent = _spiro_sys_indices(systems, spiro_pairs)
    groups = _group_by_root(parent, len(systems))
    result: list[dict] = []
    for indices in groups.values():
        result.append(
            _merged_spiro_system(rings, indices, systems)
            if len(indices) > 1 else systems[indices[0]]
        )
    return result

def build_ring_systems(mol: Mol) -> list[dict]:
    """返回环系：稠合连通 + 螺环合并。"""
    rings = _sssr(mol)
    if not rings:
        return []
    fused = _fusion_edges(rings)
    comps = _components(len(rings), fused)
    spiro = _spiro_pairs(rings)
    systems = [
        _system_entry(mol, rings, m, fused) for m in comps
    ]
    return _merge_spiro(rings, systems, spiro)
