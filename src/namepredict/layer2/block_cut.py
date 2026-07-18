"""Parent-boundary block cutting: parent atoms, side roots, connected blocks."""
from __future__ import annotations

from collections import deque

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _dbl_o_idx


def _chain_atoms(parent: dict) -> set[int]:
    return set(parent.get("chain") or [])


def _amide_n_from_c(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _add_opt(out: set[int], idx: int | None) -> set[int]:
    if idx is not None:
        out.add(idx)
    return out


def _amide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    c_idx = parent.get("amide_c_idx")
    if c_idx is None:
        return set()
    out = {c_idx}
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _amide_n_from_c(mol, c_idx))


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    kind = parent.get("kind")
    if kind in ("benzamide", "amide"):
        return _amide_fg_atoms(mol, parent)
    return set()


def parent_atom_set(parent: dict, mol: Mol) -> frozenset[int]:
    """Union chain + kind-specific FG atoms (amide C/N/O for benzamide/amide)."""
    return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))


def _is_heavy_out(atom, parent_atoms: frozenset[int]) -> bool:
    return atom.GetAtomicNum() != 1 and atom.GetIdx() not in parent_atoms


def _heavy_outside(mol: Mol, idx: int, parent_atoms: frozenset[int]) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if _is_heavy_out(n, parent_atoms)]


def side_roots(mol: Mol, parent_atoms: frozenset[int]) -> list[int]:
    """Neighbor atoms of parent_atoms that are outside the set (block roots)."""
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
    """Connected component from root not entering parent_atoms; None if empty/invalid."""
    if root in parent_atoms or mol.GetAtomWithIdx(root).GetAtomicNum() == 1:
        return None
    seen = _bfs_block(mol, root, parent_atoms)
    return frozenset(seen) if seen else None
