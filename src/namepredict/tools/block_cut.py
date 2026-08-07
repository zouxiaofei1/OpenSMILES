"""Parent-boundary block cutting: parent atoms, side roots, connected blocks."""
from __future__ import annotations

from collections import deque

from rdkit.Chem import Mol


def parent_atom_set(parent: dict, mol: Mol) -> frozenset[int]:
    """Parent boundary atoms; caller must supply a finalized parent.

    L2 parent_selector always finalizes (owned_atoms = frozenset) before L3,
    so this adapter no longer falls back to L2 parent_ownership finalization
    (keeps block_cut free of pipeline-layer dependencies).
    """
    del mol  # not needed once parent is finalized
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


def side_atoms(
    mol: Mol, owned_atoms: frozenset[int], attach_idx: int,
    seed_atoms: frozenset[int],
) -> frozenset[int]:
    """Full non-parent connected components of a substituent.

    Root the cut at the seed atoms (the substituent's already-known atoms)
    adjacent to the parent attach atom, then union their cut_block components.
    attach_idx is the parent-side atom; seed_atoms are the substituent-side
    atoms (must contain at least one atom bonded to attach_idx, otherwise the
    result is empty).
    """
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



