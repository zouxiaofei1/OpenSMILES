"""由 SSSR 构建环系拓扑（稠合图 + 杂原子汇总）：共享 >=2 个原子（邻位稠合）的环连成分量，共享 1 个原子（螺）的单环系合并为螺环系，共享 3+ 原子检测桥环（von Baeyer）系。
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.constants import C

def _sssr(mol: Mol) -> list[tuple[int, ...]]:
    """返回分子的全部 SSSR 最小环原子序列。

    `GetRingInfo().AtomRings()` 每次调用都要重建 Python 嵌套对象（实测 25~45 µs/环），
    而同一分子在一次命名内会被各层问十几次；环感知对分子恒定，故按 mol 记忆。
    """
    return memo.by_mol("sssr", lambda m: list(m.GetRingInfo().AtomRings()), mol)

def sssr_rings(mol: Mol) -> list[tuple[int, ...]]:
    """供各层统一调用的环访问器：与 `_sssr` 同一次记忆，避免各层各自重建 AtomRings。

    下游只遍历/求和，不得原地修改返回的列表（同一次命名内所有调用者共享同一对象）。
    """
    return _sssr(mol)

def _shared(a: tuple[int, ...], b: tuple[int, ...]) -> frozenset[int]:
    """返回两个环共享的原子集合。"""
    return frozenset(a) & frozenset(b)

def _ring_adjacent(ring: tuple[int, ...], a: int, b: int) -> bool:
    """若 a 与 b 在环中相邻（连续）则返回 True。"""
    n = len(ring)
    for i in range(n):
        if (ring[i] == a and ring[(i + 1) % n] == b) or \
           (ring[i] == b and ring[(i + 1) % n] == a):
            return True
    return False

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

def _is_arom_mancude(mol: Mol, atom_ids: set[int]) -> bool:
    """判断整个原子集合是否全为芳香原子。"""
    if not atom_ids:
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)

def _topology(n_rings: int, n_fusion: int, has_spiro: bool,
              is_bridged: bool = False) -> str:
    """按环数、稠合、螺与桥标志判定拓扑类型。"""
    if n_rings == 1:
        return "mono"
    if is_bridged:
        return "bridged"
    if n_fusion > 0:
        return "fused"
    if has_spiro:
        return "spiro"
    return "other"

def _non_adjacent_pairs(mrings, shared):
    """返回桥头原子：在至少一个环中不相邻的共享原子。"""
    sh_list, bh_set = sorted(shared), set()
    for i, a in enumerate(sh_list):
        for b in sh_list[i + 1:]:
            if not all(_ring_adjacent(r, a, b) for r in mrings):
                bh_set.add(a); bh_set.add(b)
    return bh_set

def _bridgeheads(rings: list[tuple[int, ...]], members: list[int],
                 shared: frozenset[int]) -> list[int]:
    """在桥环分量中寻找桥头原子。"""
    mrings = [rings[i] for i in members]
    return sorted(_non_adjacent_pairs(mrings, shared))

def _walk_path(ring, start, end, exclude, direction):
    """在环上从 start 沿 direction 走到 end，避开 exclude 收集路径原子。"""
    n, path, cur = len(ring), [], (start + direction) % len(ring)
    while cur != end and ring[cur] not in exclude:
        path.append(ring[cur]); cur = (cur + direction) % n
    return [] if cur != end else path

def _paths_between(ring: tuple[int, ...], a: int, b: int,
                   exclude: set[int]) -> list[list[int]]:
    """环中从 a 到 b 避开 exclude 集合的所有路径（不含 a、b）。"""
    try:
        ia, ib = ring.index(a), ring.index(b)
    except ValueError:
        return []
    return [p for d in (1, -1) if (p := _walk_path(ring, ia, ib, exclude, d))]

def _member_atoms(rings: list[tuple[int, ...]], members: list[int]) -> set[int]:
    """汇总分量内全部环成员的原子集合。"""
    atom_ids: set[int] = set()
    for i in members:
        atom_ids |= set(rings[i])
    return atom_ids

def _dedup_paths(paths):
    """按 frozenset 对路径去重，保持插入顺序。"""
    return list({frozenset(p): p for p in paths}.values())

def _bridge_paths(rings: list[tuple[int, ...]], members: list[int],
                  bridgeheads: list[int]) -> list[list[int]]:
    """计算桥路径（不含桥头原子的原子列表），按长度降序排列。"""
    if len(bridgeheads) < 2:
        return []
    bh_set, a, b = set(bridgeheads), bridgeheads[0], bridgeheads[1]
    mrings = [rings[i] for i in members]
    all_paths = [p for r in mrings for p in _paths_between(r, a, b, bh_set)]
    return sorted(_dedup_paths(all_paths), key=len, reverse=True)

def _bridge_info(rings: list[tuple[int, ...]], members: list[int],
                 bridgeheads: list[int]) -> dict | None:
    """计算桥环分量的桥信息。"""
    if len(bridgeheads) < 2:
        return None
    paths = _bridge_paths(rings, members, bridgeheads)
    return {"bridgeheads": bridgeheads,
            "bridge_lengths": [len(p) for p in paths],
            "bridge_paths": paths}

def _member_edges(
    members: list[int], fusion_edges: list[tuple[int, int, frozenset[int]]],
) -> list[tuple[int, int, list[int]]]:
    """筛选仅连接分量内环的稠合边。"""
    mset = set(members)
    return [
        (i, j, sorted(sh)) for i, j, sh in fusion_edges if i in mset and j in mset
    ]

def _system_dict(
    atom_ids: set[int], members: list[int], edges: list, mol: Mol, spiro_in: bool,
    is_bridged: bool = False, bridge_info: dict | None = None,
) -> dict:
    """组装单个环系的事实 dict。"""
    result = {
        "atom_ids": sorted(atom_ids),
        "sssr_indices": sorted(members),
        "fusion_edges": edges,
        "n_rings": len(members),
        "n_atoms": len(atom_ids),
        "hetero_atoms": _hetero_atoms(mol, atom_ids),
        "is_aromatic_mancude": _is_arom_mancude(mol, atom_ids),
        "topology": _topology(len(members), len(edges),
                              spiro_in and len(members) == 1, is_bridged),
    }
    if bridge_info:
        result["bridgeheads"] = bridge_info["bridgeheads"]
        result["bridge_lengths"] = bridge_info["bridge_lengths"]
        result["bridge_paths"] = bridge_info["bridge_paths"]
    return result

def _try_bridge_component(rings, members, fusion_edges):
    """尝试寻找桥信息：在双环分量中遍历 3+ 共享原子的稠合边。"""
    mset = set(members)
    for i, j, sh in fusion_edges:
        if i in mset and j in mset and len(sh) >= 3:
            bh = _bridgeheads(rings, members, sh)
            if len(bh) >= 2:
                return _bridge_info(rings, members, bh)
    return None

def _compute_bridged_info(
    rings: list[tuple[int, ...]], members: list[int],
    fusion_edges: list[tuple[int, int, frozenset[int]]],
) -> tuple[bool, dict | None]:
    """检查分量是否桥环并计算桥信息。"""
    if len(members) != 2:
        return False, None
    bi = _try_bridge_component(rings, members, fusion_edges)
    return (True, bi) if bi else (False, None)

def _system_entry(
    mol: Mol,
    rings: list[tuple[int, ...]],
    members: list[int],
    fusion_edges: list[tuple[int, int, frozenset[int]]],
    spiro_in: bool,
) -> dict:
    """为连通分量构建环系条目（含桥检测）。"""
    atoms = _member_atoms(rings, members)
    edges = _member_edges(members, fusion_edges)
    is_bridged, bridge_info = _compute_bridged_info(rings, members, fusion_edges)
    return _system_dict(atoms, members, edges, mol, spiro_in, is_bridged, bridge_info)

def _spiro_touching(members: list[int], spiro: list[tuple[int, int, int]]) -> bool:
    """判断分量是否与某条螺共享配对相连。"""
    mset = set(members)
    return any(i in mset or j in mset for i, j, _ in spiro)

def _collect_merged_fields(
    systems: list[dict], indices: list[int],
) -> tuple[set[int], set[int], list[dict]]:
    """收集待合并系统的 sssr_indices、atom_ids、hetero_atoms。"""
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
    mol: Mol, rings: list[tuple[int, ...]],
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
        "is_aromatic_mancude": _is_arom_mancude(mol, atoms),
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
    mol: Mol, rings: list[tuple[int, ...]],
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
            _merged_spiro_system(mol, rings, indices, systems)
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
        _system_entry(mol, rings, m, fused, _spiro_touching(m, spiro))
        for m in comps
    ]
    return _merge_spiro(mol, rings, systems, spiro)
