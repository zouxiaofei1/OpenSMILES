"""L2 topology for pure saturated-carbon rooted side trees."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from rdkit.Chem import Mol


@dataclass(frozen=True)
class RootedAlkylTree:
    root: int
    atoms: frozenset[int]
    children: dict[int, tuple[int, ...]]  # parent -> sorted child idxs
    depth: dict[int, int]  # atom -> depth from root


def _is_pure_sat_c(mol: Mol, idx: int, *, atoms: frozenset[int]) -> bool:
    """True when idx is non-ring C and all *in-claim* heavy edges are single C–C."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetIsAromatic():
        return False
    for n in atom.GetNeighbors():
        j = n.GetIdx()
        if n.GetAtomicNum() == 1 or j not in atoms:
            continue  # H or outside claim (parent attach) ignored
        if n.GetAtomicNum() != 6 or n.GetIsAromatic() or n.IsInRing():
            return False
        bond = mol.GetBondBetweenAtoms(idx, j)
        if bond is None or bond.GetBondTypeAsDouble() != 1.0:
            return False
    return True


def _c_neighbors_in(mol: Mol, idx: int, allowed: frozenset[int]) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetIdx() in allowed and n.GetAtomicNum() == 6
    ]


def _bfs_tree(
    mol: Mol, root: int, atoms: frozenset[int],
) -> tuple[dict[int, tuple[int, ...]], dict[int, int]] | None:
    children: dict[int, list[int]] = {a: [] for a in atoms}
    depth = {root: 0}
    q: deque[int] = deque([root])
    seen = {root}
    while q:
        cur = q.popleft()
        for n in _c_neighbors_in(mol, cur, atoms):
            if n in seen:
                continue
            seen.add(n)
            depth[n] = depth[cur] + 1
            children[cur].append(n)
            q.append(n)
    if seen != set(atoms):
        return None
    frozen = {k: tuple(sorted(v)) for k, v in children.items()}
    return frozen, depth


def build_rooted_alkyl_tree(
    mol: Mol,
    *,
    root: int,
    atoms: frozenset[int],
    max_atoms: int = 12,
    max_depth: int = 3,
) -> RootedAlkylTree | None:
    """Accept pure saturated non-ring carbon tree rooted at claim.root."""
    if root not in atoms or not atoms or len(atoms) > max_atoms:
        return None
    if any(not _is_pure_sat_c(mol, i, atoms=atoms) for i in atoms):
        return None
    built = _bfs_tree(mol, root, atoms)
    if built is None:
        return None
    children, depth = built
    if max(depth.values()) > max_depth:
        return None
    return RootedAlkylTree(root=root, atoms=atoms, children=children, depth=depth)
