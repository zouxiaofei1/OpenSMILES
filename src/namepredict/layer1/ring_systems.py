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



def build_ring_systems(mol: Mol) -> list[dict]:
    """返回环系：稠合连通 + 螺环合并。"""
    rings = _sssr(mol)
    if not rings:
        return []
    fused, spiro = _ring_pairs(rings)
    comps = _components(len(rings), fused)
    return [ _system_entry(mol, rings, m, fused) for m in comps]
