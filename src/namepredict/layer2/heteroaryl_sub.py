"""Heteroaryl side chains: unsubstituted pyridin-n-yl (P-29).

Parent is open chain; pyridine N=1; attachment locant via ring walk.
"""
from __future__ import annotations

from rdkit.Chem import Mol


def _nb_out(mol: Mol, i: int, ring: set[int]) -> list:
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ring
    ]


def _is_pyridine_ring(mol: Mol, r: tuple) -> bool:
    if len(r) != 6:
        return False
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in r]
    if zs.count(7) != 1 or any(z not in (6, 7) for z in zs):
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)


def _pyridine_rings(mol: Mol) -> list[set[int]]:
    return [
        set(r) for r in mol.GetRingInfo().AtomRings() if _is_pyridine_ring(mol, r)
    ]


def _n_idx(mol: Mol, ring: set[int]) -> int:
    return next(i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)


def _unsub_ok(mol: Mol, ring: set[int], attach: int, parent: int) -> bool:
    """No exo heavy except parent link; attach is ring C not N."""
    if mol.GetAtomWithIdx(attach).GetAtomicNum() != 6:
        return False
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if not (i == attach and nb.GetIdx() == parent):
                return False
    return True


def _pyridine_at(mol: Mol, attach: int, parent: int) -> set[int] | None:
    hits = [r for r in _pyridine_rings(mol) if attach in r and parent not in r]
    if len(hits) != 1:
        return None
    ring = hits[0]
    return ring if _unsub_ok(mol, ring, attach, parent) else None


def _ring_nbrs(mol: Mol, i: int, ring: set[int]) -> list[int]:
    return sorted(
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in ring
    )


def _walk(mol: Mol, ring: set[int], start: int, nxt: int) -> list[int]:
    path, prev, cur = [start], start, nxt
    while cur != start:
        path.append(cur)
        opts = [j for j in _ring_nbrs(mol, cur, ring) if j != prev]
        if not opts:
            break
        prev, cur = cur, opts[0]
    return path


def _locant(mol: Mol, ring: set[int], attach: int) -> int:
    n1 = _n_idx(mol, ring)
    nbrs = _ring_nbrs(mol, n1, ring)
    if len(nbrs) != 2:
        return 1
    a = _walk(mol, ring, n1, nbrs[0])
    b = _walk(mol, ring, n1, nbrs[1])
    la = a.index(attach) + 1 if attach in a else 99
    lb = b.index(attach) + 1 if attach in b else 99
    return min(la, lb)


def _parent_link(mol: Mol, start: int, parent: set[int]) -> int | None:
    links = [
        n.GetIdx() for n in mol.GetAtomWithIdx(start).GetNeighbors()
        if n.GetIdx() in parent
    ]
    return links[0] if len(links) == 1 else None


def _make_pyridinyl(att: int, start: int, ring: set[int], loc: int) -> dict:
    return {
        "attach": att, "outer_c": start, "ring": ring, "atoms": list(ring),
        "en": f"pyridin-{loc}-yl", "zh": f"吡啶-{loc}-基", "paren": True,
    }


def _from_start(mol: Mol, start: int, parent: set[int]) -> dict | None:
    if start in parent:
        return None
    att = _parent_link(mol, start, parent)
    if att is None:
        return None
    ring = _pyridine_at(mol, start, att)
    if ring is None:
        return None
    return _make_pyridinyl(att, start, ring, _locant(mol, ring, start))


def _side_starts(mol: Mol, parent: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for r in parent
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() not in parent
    ]


def ring_pyridinyls(mol: Mol, parent: set[int]) -> list[dict]:
    out: list[dict] = []
    for s in _side_starts(mol, parent):
        one = _from_start(mol, s, parent)
        if one is not None:
            out.append(one)
    return out


def pyridinyl_atoms(mol: Mol, parent: set[int]) -> set[int]:
    out: set[int] = set()
    for p in ring_pyridinyls(mol, parent):
        out |= set(p["atoms"])
    return out
