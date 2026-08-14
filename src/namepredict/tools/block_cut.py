"""母体边界块切割：母体原子、侧链根、连通块。"""
from __future__ import annotations

from collections import deque

from rdkit.Chem import Mol


def parent_atom_set(parent: dict, mol: Mol) -> frozenset[int]:
    del mol  # 母体定稿后不再需要
    owned = parent.get("owned_atoms")
    if not isinstance(owned, frozenset):
        raise ValueError(f"unfinalized parent (missing owned_atoms): kind={parent.get('kind')!r}")
    return owned


def _is_heavy_out(atom, parent_atoms: frozenset[int]) -> bool:
    return atom.GetAtomicNum() != 1 and atom.GetIdx() not in parent_atoms


def _heavy_outside(mol: Mol, idx: int, parent_atoms: frozenset[int]) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if _is_heavy_out(n, parent_atoms)]


def side_roots(mol: Mol, parent_atoms: frozenset[int]) -> list[int]:
    """parent_atoms 的邻居原子中位于该集合之外者（块根）。"""
    roots: set[int] = set()
    for p in parent_atoms:
        roots.update(_heavy_outside(mol, p, parent_atoms))
    return sorted(roots)


def _visit(mol: Mol, cur: int, parent_atoms: frozenset[int], seen: set[int], q: deque) -> None:
    if cur in seen or cur in parent_atoms:
        return
    if mol.GetAtomWithIdx(cur).GetAtomicNum() == 1:
        return
    seen.add(cur)
    q.extend(_heavy_outside(mol, cur, parent_atoms))


def _bfs_block(mol: Mol, root: int, parent_atoms: frozenset[int]) -> set[int]:
    seen: set[int] = set()
    q: deque[int] = deque([root])
    while q:
        _visit(mol, q.popleft(), parent_atoms, seen, q)
    return seen


def cut_block(mol: Mol, root: int, parent_atoms: frozenset[int]) -> frozenset[int] | None:
    """从 root 出发不进入 parent_atoms 的连通分量；为空/无效时返回 None。"""
    if root in parent_atoms or mol.GetAtomWithIdx(root).GetAtomicNum() == 1:
        return None
    seen = _bfs_block(mol, root, parent_atoms)
    return frozenset(seen) if seen else None


def side_atoms(
    mol: Mol, owned_atoms: frozenset[int], attach_idx: int,
    seed_atoms: frozenset[int],
) -> frozenset[int]:
    roots = {
        n.GetIdx()
        for n in mol.GetAtomWithIdx(attach_idx).GetNeighbors()
        if n.GetIdx() not in owned_atoms and n.GetIdx() in seed_atoms
    }
    comp: set[int] = set()
    for r in roots:
        block = cut_block(mol, r, owned_atoms)
        if block:
            comp |= block
    return frozenset(comp)



