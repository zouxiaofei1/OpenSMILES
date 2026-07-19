"""Ring-system topology from SSSR (fusion graph + hetero summary).

Builds connected components where rings share >=2 atoms (ortho-fusion).
Mono-ring systems sharing exactly 1 atom (spiro) are merged into one spiro system.
"""
from __future__ import annotations

from rdkit.Chem import Mol


def _sssr(mol: Mol) -> list[tuple[int, ...]]:
    return list(mol.GetRingInfo().AtomRings())


def _shared(a: tuple[int, ...], b: tuple[int, ...]) -> frozenset[int]:
    return frozenset(a) & frozenset(b)


def _fusion_edges(rings: list[tuple[int, ...]]) -> list[tuple[int, int, frozenset[int]]]:
    """Edges (i, j, shared_atoms) for pairs sharing >=2 atoms."""
    out: list[tuple[int, int, frozenset[int]]] = []
    for i, ri in enumerate(rings):
        for j in range(i + 1, len(rings)):
            sh = _shared(ri, rings[j])
            if len(sh) >= 2:
                out.append((i, j, sh))
    return out


def _spiro_pairs(rings: list[tuple[int, ...]]) -> list[tuple[int, int, int]]:
    """Pairs sharing exactly 1 atom: (i, j, atom)."""
    out: list[tuple[int, int, int]] = []
    for i, ri in enumerate(rings):
        for j in range(i + 1, len(rings)):
            sh = _shared(ri, rings[j])
            if len(sh) == 1:
                out.append((i, j, next(iter(sh))))
    return out


def _uf_find(parent: list[int], x: int) -> int:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _uf_union(parent: list[int], a: int, b: int) -> None:
    ra, rb = _uf_find(parent, a), _uf_find(parent, b)
    if ra != rb:
        parent[rb] = ra


def _components(n: int, edges: list[tuple[int, int, frozenset[int]]]) -> list[list[int]]:
    """Union-find components over ring indices."""
    parent = list(range(n))
    for i, j, _ in edges:
        _uf_union(parent, i, j)
    buckets: dict[int, list[int]] = {}
    for i in range(n):
        buckets.setdefault(_uf_find(parent, i), []).append(i)
    return list(buckets.values())


def _hetero_atoms(mol: Mol, atom_ids: set[int]) -> list[dict]:
    out: list[dict] = []
    for i in sorted(atom_ids):
        z = mol.GetAtomWithIdx(i).GetAtomicNum()
        if z != 6:
            out.append({"idx": i, "Z": z})
    return out


def _is_arom_mancude(mol: Mol, atom_ids: set[int]) -> bool:
    if not atom_ids:
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)


def _topology(n_rings: int, n_fusion: int, has_spiro: bool) -> str:
    if n_rings == 1:
        return "mono"
    if n_fusion > 0:
        return "fused"
    if has_spiro:
        return "spiro"
    return "other"


def _member_atoms(rings: list[tuple[int, ...]], members: list[int]) -> set[int]:
    atom_ids: set[int] = set()
    for i in members:
        atom_ids |= set(rings[i])
    return atom_ids


def _member_edges(
    members: list[int], fusion_edges: list[tuple[int, int, frozenset[int]]],
) -> list[tuple[int, int, list[int]]]:
    mset = set(members)
    return [
        (i, j, sorted(sh)) for i, j, sh in fusion_edges if i in mset and j in mset
    ]


def _system_dict(
    atom_ids: set[int], members: list[int], edges: list, mol: Mol, spiro_in: bool,
) -> dict:
    return {
        "atom_ids": sorted(atom_ids),
        "sssr_indices": sorted(members),
        "fusion_edges": edges,
        "n_rings": len(members),
        "n_atoms": len(atom_ids),
        "hetero_atoms": _hetero_atoms(mol, atom_ids),
        "is_aromatic_mancude": _is_arom_mancude(mol, atom_ids),
        "topology": _topology(len(members), len(edges), spiro_in and len(members) == 1),
    }


def _system_entry(
    mol: Mol,
    rings: list[tuple[int, ...]],
    members: list[int],
    fusion_edges: list[tuple[int, int, frozenset[int]]],
    spiro_in: bool,
) -> dict:
    atoms = _member_atoms(rings, members)
    edges = _member_edges(members, fusion_edges)
    return _system_dict(atoms, members, edges, mol, spiro_in)


def _spiro_touching(members: list[int], spiro: list[tuple[int, int, int]]) -> bool:
    mset = set(members)
    return any(i in mset or j in mset for i, j, _ in spiro)


def _collect_merged_fields(
    systems: list[dict], indices: list[int],
) -> tuple[set[int], set[int], list[dict]]:
    """Collect sssr_indices, atom_ids, hetero_atoms from systems to merge."""
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
    """Merge several mono-ring systems sharing a spiro atom into one spiro system."""
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
    """Union-find parent array for systems connected by spiro atoms."""
    n = len(systems)
    parent = list(range(n))
    for i, j, _ in spiro_pairs:
        si = _find_sys_for_ring(systems, i)
        sj = _find_sys_for_ring(systems, j)
        if si >= 0 and sj >= 0 and si != sj:
            _uf_union(parent, si, sj)
    return parent


def _find_sys_for_ring(systems: list[dict], ring_idx: int) -> int:
    for idx, s in enumerate(systems):
        if ring_idx in s["sssr_indices"]:
            return idx
    return -1


def _group_by_root(parent: list[int], n: int) -> dict[int, list[int]]:
    buckets: dict[int, list[int]] = {}
    for idx in range(n):
        buckets.setdefault(_uf_find(parent, idx), []).append(idx)
    return buckets


def _merge_spiro(
    mol: Mol, rings: list[tuple[int, ...]],
    systems: list[dict], spiro_pairs: list[tuple[int, int, int]],
) -> list[dict]:
    """Merge mono-ring systems that share a spiro atom into spiro systems."""
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
    """Return ring systems: fusion-connected + spiro-merged."""
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
