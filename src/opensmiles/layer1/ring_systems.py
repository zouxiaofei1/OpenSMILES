"""环系拓扑 Layer1,4,5共用"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from opensmiles.tools import memo
from opensmiles.constants import C

def _kekulize_uncached(mol: Mol) -> Mol | None:
    """Kekulize 并清芳香标志的分子副本；失败返 None。"""
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        return None
    return kek


def kekulized(mol: Mol) -> Mol | None:
    """分子的 Kekulize 副本（各层统一入口，按分子记忆）；失败返 None。"""
    return memo.by_mol("kekulized", _kekulize_uncached, mol)

def sssr_rings(mol: Mol) -> list[tuple[int, ...]]:
    """分子的全部 SSSR 最小环原子序列（各层统一入口，按分子记忆）。"""
    return memo.by_mol("sssr", lambda m: list(m.GetRingInfo().AtomRings()), mol)

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

def _components(n: int, edges: list[tuple[int, int, frozenset[int]]]) -> list[list[int]]:
    """环索引的连通分量，按最小成员升序、分量内成员升序（环数规模小，不需并查集）。"""
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for i, j, _ in edges:
        adj[i].append(j)
        adj[j].append(i)
    seen: set[int] = set()
    out: list[list[int]] = []
    for start in range(n):
        if start in seen:
            continue
        comp: list[int] = []
        stack = [start]
        while stack:
            v = stack.pop()
            if v in seen:
                continue
            seen.add(v)
            comp.append(v)
            stack.extend(adj[v])
        out.append(sorted(comp))
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


def _free_spiro_atoms(mol: Mol, atoms: set[int],
                      spiro_edges: list[tuple[int, int, int]]) -> tuple[int, ...]:
    """分量的自由螺原子：去掉后环子图断成 >=2 分量（P-24.1 自由螺连接）。"""
    if not spiro_edges:
        return ()
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
        "hetero_atoms": [{"idx": i, "Z": z} for i in sorted(atom_ids)
                         if (z := mol.GetAtomWithIdx(i).GetAtomicNum()) not in (1, C)],
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
                        _free_spiro_atoms(mol, atoms, spiro_edges))



def build_ring_systems(mol: Mol) -> list[dict]:
    """返回环系：稠合连通 + 螺环合并。"""
    rings = sssr_rings(mol)
    if not rings:
        return []
    fused, spiro = _ring_pairs(rings)
    comps = _components(len(rings), fused + spiro)  # 螺环对并入连通分量
    return [ _system_entry(mol, rings, m, fused, spiro) for m in comps]
