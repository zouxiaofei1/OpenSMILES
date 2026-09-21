"""环系拓扑 Layer1,4,5共用"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.constants import C

def kekulized(mol: Mol) -> Mol | None:
    """返回 Kekulize 并清芳香标志后的分子副本；失败返 None。"""
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        return None
    return kek

def _sssr(mol: Mol) -> list[tuple[int, ...]]:
    """返回分子的全部 SSSR 最小环原子序列。"""
    return memo.by_mol("sssr", lambda m: list(m.GetRingInfo().AtomRings()), mol)

def sssr_rings(mol: Mol) -> list[tuple[int, ...]]:
    """供各层统一调用的环访问器（与 _sssr 同一次记忆）。"""
    return _sssr(mol)

def _ring_pairs(rings: list[tuple[int, ...]]) -> tuple[list, list]:
    """单遍扫全部环对：共享 >=2 原子为稠合边、恰好 1 个为螺环对。"""
    sets = [frozenset(r) for r in rings]  # 各环原子集只构造一次，供全部 O(R²) 对交集复用
    fused: list[tuple[int, int, frozenset[int]]] = []
    spiro: list[tuple[int, int, int]] = []
    for i, si in enumerate(sets):
        for j in range(i + 1, len(sets)):
            sh = si & sets[j]
            if len(sh) >= 2:
                fused.append((i, j, sh))
            elif len(sh) == 1:
                spiro.append((i, j, next(iter(sh))))
    return fused, spiro

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

def _split_count(adj: dict[int, list[int]], removed: int) -> int:
    """去掉 removed 后，其环内邻居分属的连通分量数。"""
    seen: set[int] = set()
    count = 0
    for start in adj.get(removed) or ():
        if start in seen:
            continue
        count += 1
        stack = [start]
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            stack.extend(u for u in adj[v] if u != removed)
    return count


def _free_spiro_atoms(mol: Mol, rings: list[tuple[int, ...]], members: list[int],
                      spiro_edges: list[tuple[int, int, int]]) -> tuple[int, ...]:
    """分量的自由螺原子：去掉后环子图断成 >=2 分量（P-24.1 自由螺连接）。"""
    if not spiro_edges:
        return ()
    atoms: set[int] = set()
    for i in members:
        atoms |= set(rings[i])
    adj = {v: sorted(n.GetIdx() for n in mol.GetAtomWithIdx(v).GetNeighbors()
                     if n.GetIdx() in atoms) for v in atoms}
    out: list[int] = []
    for _, _, s in spiro_edges:
        if s not in out and _split_count(adj, s) >= 2:
            out.append(s)
    return tuple(sorted(out))


def _system_dict(
    atom_ids: set[int], members: list[int], edges: list, mol: Mol,
    spiro_edges: list | None = None, free_spiro: tuple[int, ...] = (),
) -> dict:
    """组装单个环系的事实 dict（螺环系附 spiro_edges 与自由螺原子）。"""
    return {
        "atom_ids": sorted(atom_ids),
        "sssr_indices": sorted(members),
        "fusion_edges": edges,
        "spiro_edges": list(spiro_edges or ()),
        "free_spiro_atoms": list(free_spiro),
        "n_rings": len(members),
        "n_atoms": len(atom_ids),
        "hetero_atoms": None,
        "is_aromatic_mancude": None,
        "topology": "spiro" if free_spiro else None,
    }

def _system_entry(
    mol: Mol,
    rings: list[tuple[int, ...]],
    members: list[int],
    fusion_edges: list[tuple[int, int, frozenset[int]]],
    spiro_pairs: list[tuple[int, int, int]] = (),
) -> dict:
    """为连通分量构建环系条目（含螺环合并后的自由螺原子判定）。"""
    atoms = _member_atoms(rings, members)
    edges = _member_edges(members, fusion_edges)
    mset = set(members)
    spiro_edges = [(i, j, s) for i, j, s in spiro_pairs if i in mset and j in mset]
    return _system_dict(atoms, members, edges, mol, spiro_edges,
                        _free_spiro_atoms(mol, rings, members, spiro_edges))



def build_ring_systems(mol: Mol) -> list[dict]:
    """返回环系：稠合连通 + 螺环合并。"""
    rings = _sssr(mol)
    if not rings:
        return []
    fused, spiro = _ring_pairs(rings)
    comps = _components(len(rings), fused + spiro)  # 螺环对并入连通分量
    return [ _system_entry(mol, rings, m, fused, spiro) for m in comps]
