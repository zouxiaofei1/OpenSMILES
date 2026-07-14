"""Depth-1 unsubstituted (or mono-halo) phenyl / phenoxy topology (P-29.3).

Used by L2 parent selection (allowed atoms / ring pick) and L3 extractors.
Nested aryl (Ph-on-Ph) is intentionally rejected this round (depth=1).
"""
from __future__ import annotations

from rdkit.Chem import Mol

_HALO = frozenset({9, 17, 35, 53})
_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}


def _is_arom_c6(mol: Mol, atoms) -> bool:
    if len(atoms) != 6:
        return False
    return all(
        mol.GetAtomWithIdx(i).GetIsAromatic()
        and mol.GetAtomWithIdx(i).GetAtomicNum() == 6
        for i in atoms
    )


def _c6_rings_at(mol: Mol, c_idx: int) -> list[set[int]]:
    return [
        set(r) for r in mol.GetRingInfo().AtomRings()
        if c_idx in r and _is_arom_c6(mol, r)
    ]


def _nb_outside(mol: Mol, i: int, ring: set[int]):
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ring
    ]


def _is_terminal_halo(atom) -> bool:
    if atom.GetAtomicNum() not in _HALO:
        return False
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() != 1) == 1


def _count_side_halos(mol: Mol, ring: set[int], attach: int, parent: int) -> int | None:
    n = 0
    for i in ring:
        for nb in _nb_outside(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if not _is_terminal_halo(nb):
                return None
            n += 1
    return n


def _side_ok(mol: Mol, ring: set[int], attach: int, parent: int) -> bool:
    n = _count_side_halos(mol, ring, attach, parent)
    return n is not None and n <= 1


def _phenyl_at(mol: Mol, attach: int, parent: int) -> set[int] | None:
    hits = _c6_rings_at(mol, attach)
    if len(hits) != 1 or parent in hits[0]:
        return None
    ring = hits[0]
    return ring if _side_ok(mol, ring, attach, parent) else None


def _halo_pair(mol: Mol, ring: set[int]) -> tuple[int, int] | None:
    found = [
        (i, nb.GetAtomicNum())
        for i in ring
        for nb in _nb_outside(mol, i, ring)
        if _is_terminal_halo(nb)
    ]
    return found[0] if len(found) == 1 else None


def _ring_nbrs(mol: Mol, i: int, ring: set[int]) -> list[int]:
    return sorted(
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetIdx() in ring
    )


def _walk_ring(mol: Mol, ring: set[int], start: int, nxt: int) -> list[int]:
    path, prev, cur = [start], start, nxt
    while cur != start:
        path.append(cur)
        opts = [j for j in _ring_nbrs(mol, cur, ring) if j != prev]
        if not opts:
            break
        prev, cur = cur, opts[0]
    return path


def _ring_order(mol: Mol, ring: set[int], start: int) -> list[int]:
    nbrs = _ring_nbrs(mol, start, ring)
    if len(nbrs) != 2:
        return [start]
    a = _walk_ring(mol, ring, start, nbrs[0])
    b = _walk_ring(mol, ring, start, nbrs[1])
    return a if a <= b else b


def _ph_locant(mol: Mol, ring: set[int], attach: int, site: int) -> int:
    order = _ring_order(mol, ring, attach)
    return order.index(site) + 1 if site in order else 1


def _phenoxy_from_ether(mol: Mol, e: dict, parent: set[int]) -> dict | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in parent) == (c2 in parent):
        return None
    ring_c, outer = (c1, c2) if c1 in parent else (c2, c1)
    ph = _phenyl_at(mol, outer, o)
    if ph is None:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "ph": ph,
            "atoms": [o, *ph]}


def _parent_link(mol: Mol, start: int, parent: set[int]) -> int | None:
    links = [
        n.GetIdx() for n in mol.GetAtomWithIdx(start).GetNeighbors()
        if n.GetIdx() in parent
    ]
    return links[0] if len(links) == 1 else None


def _phenyl_from_start(mol: Mol, start: int, parent: set[int]) -> dict | None:
    if start in parent:
        return None
    att = _parent_link(mol, start, parent)
    if att is None:
        return None
    ph = _phenyl_at(mol, start, att)
    return None if ph is None else {
        "attach": att, "outer_c": start, "ph": ph, "atoms": list(ph),
    }


def _ring_phenoxys(info: dict, parent: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _phenoxy_from_ether(mol, e, parent)
        if one is not None:
            out.append(one)
    return out


def _side_c_starts(mol: Mol, parent: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for r in parent
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() not in parent
    ]


def _ring_phenyls(mol: Mol, parent: set[int]) -> list[dict]:
    out: list[dict] = []
    for s in _side_c_starts(mol, parent):
        one = _phenyl_from_start(mol, s, parent)
        if one is not None:
            out.append(one)
    return out


def _halo_atoms_on(mol: Mol, ph: set[int]) -> set[int]:
    return {
        nb.GetIdx()
        for i in ph
        for nb in _nb_outside(mol, i, ph)
        if _is_terminal_halo(nb)
    }


def _aryl_atoms(info: dict, parent: set[int]) -> set[int]:
    mol: Mol = info["mol"]
    out: set[int] = set()
    for p in _ring_phenoxys(info, parent):
        out.update(p["atoms"])
        out.update(_halo_atoms_on(mol, p["ph"]))
    for p in _ring_phenyls(mol, parent):
        out.update(p["atoms"])
        out.update(_halo_atoms_on(mol, p["ph"]))
    return out


def _aryl_sub_n(info: dict, parent: set[int]) -> int:
    mol: Mol = info["mol"]
    return len(_ring_phenyls(mol, parent)) + len(_ring_phenoxys(info, parent))


def _halo_ph_names(mol: Mol, ph: set[int], attach: int, stem_en: str, stem_zh: str):
    hp = _halo_pair(mol, ph)
    if hp is None:
        return stem_en, stem_zh, False
    site, z = hp
    loc = _ph_locant(mol, ph, attach, site)
    return f"{loc}-{_HALO_EN[z]}{stem_en}", f"{loc}-{_HALO_ZH[z]}{stem_zh}", True


def _phenyl_name(mol: Mol, ph: set[int], attach: int) -> tuple[str, str, bool]:
    return _halo_ph_names(mol, ph, attach, "phenyl", "苯基")


def _phenoxy_name(mol: Mol, ph: set[int], outer: int) -> tuple[str, str, bool]:
    return _halo_ph_names(mol, ph, outer, "phenoxy", "苯氧基")


def _arom_c6_rings(mol: Mol) -> list[set[int]]:
    return [set(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]


def _is_unfused_benzene_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _phenyl_starts_set(mol: Mol, parent: set[int]) -> set[int]:
    return {p["outer_c"] for p in _ring_phenyls(mol, parent)}


def _unsub_phenyl_at(mol: Mol, c_idx: int, n_idx: int) -> bool:
    """True if c_idx is sole N-attachment of unsubstituted phenyl (amide N-Ph)."""
    if mol.GetAtomWithIdx(n_idx).GetAtomicNum() != 7:
        return False
    ph = _phenyl_at(mol, c_idx, n_idx)
    return ph is not None and _halo_pair(mol, ph) is None


def _arom_c6_ring_lists(mol: Mol) -> list[list[int]]:
    return [list(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]
