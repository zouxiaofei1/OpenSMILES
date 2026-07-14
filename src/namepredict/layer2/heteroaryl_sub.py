"""Heteroaryl side chains: pyridin-n-yl with simple leaves (P-29).

Parent is open chain; pyridine N=1; attachment locant via ring walk.
Leaves on ring C: halo / Me / methoxy / ethoxy (≤3).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import _alkoxy_n
from namepredict.layer2.leaves.topo import is_halo, is_me_leaf


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


def _leaf_ok(mol: Mol, nb, ring_i: int) -> bool:
    if is_halo(nb):
        return True
    if nb.GetAtomicNum() == 6 and is_me_leaf(mol, nb.GetIdx(), ring_i):
        return True
    return _alkoxy_n(mol, nb, ring_i) in (1, 2)


def _count_leaves(mol: Mol, ring: set[int], attach: int, parent: int) -> int | None:
    n = 0
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 7 or not _leaf_ok(mol, nb, i):
                return None
            n += 1
    return n


def _side_ok(mol: Mol, ring: set[int], attach: int, parent: int) -> bool:
    if mol.GetAtomWithIdx(attach).GetAtomicNum() != 6:
        return False
    n = _count_leaves(mol, ring, attach, parent)
    return n is not None and n <= 3


def _pyridine_at(mol: Mol, attach: int, parent: int) -> set[int] | None:
    hits = [r for r in _pyridine_rings(mol) if attach in r and parent not in r]
    if len(hits) != 1:
        return None
    ring = hits[0]
    return ring if _side_ok(mol, ring, attach, parent) else None


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


def _ring_order(mol: Mol, ring: set[int], attach: int) -> list[int]:
    """N=1, walk so attach gets lowest locant (2 if ortho)."""
    n1 = _n_idx(mol, ring)
    nbrs = _ring_nbrs(mol, n1, ring)
    if len(nbrs) != 2:
        return [n1]
    a = _walk(mol, ring, n1, nbrs[0])
    b = _walk(mol, ring, n1, nbrs[1])
    la = a.index(attach) + 1 if attach in a else 99
    lb = b.index(attach) + 1 if attach in b else 99
    return a if la <= lb else b


def _locant_in(order: list[int], site: int) -> int:
    return order.index(site) + 1 if site in order else 1


def _locant(mol: Mol, ring: set[int], attach: int) -> int:
    return _locant_in(_ring_order(mol, ring, attach), attach)


def _parent_link(mol: Mol, start: int, parent: set[int]) -> int | None:
    links = [
        n.GetIdx() for n in mol.GetAtomWithIdx(start).GetNeighbors()
        if n.GetIdx() in parent
    ]
    return links[0] if len(links) == 1 else None


_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}


def _leaf_name(mol: Mol, nb, ring_i: int) -> tuple[str, str] | None:
    if is_halo(nb):
        z = nb.GetAtomicNum()
        return _HALO_EN[z], _HALO_ZH[z]
    if nb.GetAtomicNum() == 6 and is_me_leaf(mol, nb.GetIdx(), ring_i):
        return "methyl", "甲基"
    return _alkoxy_leaf(mol, nb, ring_i)


def _alkoxy_leaf(mol: Mol, nb, ring_i: int) -> tuple[str, str] | None:
    n = _alkoxy_n(mol, nb, ring_i)
    if n == 1:
        return "methoxy", "甲氧基"
    if n == 2:
        return "ethoxy", "乙氧基"
    return None


def _leaf_items(
    mol: Mol, ring: set[int], attach: int, parent: int, order: list[int],
) -> list[tuple[int, str, str]]:
    items: list[tuple[int, str, str]] = []
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            got = _leaf_name(mol, nb, i)
            if got:
                items.append((_locant_in(order, i), got[0], got[1]))
    return items


def _group_pref(en: str, locs_zh: list[tuple[int, str]]) -> tuple[str, str]:
    locs = sorted(l for l, _ in locs_zh)
    zh = locs_zh[0][1]
    ls = ",".join(str(l) for l in locs)
    m = {1: "", 2: "di", 3: "tri"}.get(len(locs), "")
    mz = {1: "", 2: "二", 3: "三"}.get(len(locs), "")
    return f"{ls}-{m}{en}", f"{ls}-{mz}{zh}"


def _pref_from_items(items: list[tuple[int, str, str]]) -> tuple[str, str]:
    if not items:
        return "", ""
    from collections import defaultdict
    g: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for loc, en, zh in items:
        g[en].append((loc, zh))
    parts = [_group_pref(en, g[en]) for en in sorted(g)]
    return "-".join(p[0] for p in parts), "-".join(p[1] for p in parts)


def _make_pyridinyl(
    mol: Mol, att: int, start: int, ring: set[int], parent: int,
) -> dict:
    order = _ring_order(mol, ring, start)
    loc = _locant_in(order, start)
    pe, pz = _pref_from_items(_leaf_items(mol, ring, start, att, order))
    en = f"{pe}pyridin-{loc}-yl" if pe else f"pyridin-{loc}-yl"
    zh = f"{pz}吡啶-{loc}-基" if pz else f"吡啶-{loc}-基"
    return {
        "attach": att, "outer_c": start, "ring": ring, "atoms": list(ring),
        "en": en, "zh": zh, "paren": True,
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
    return _make_pyridinyl(mol, att, start, ring, att)


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


def _add_alkoxy_atoms(mol: Mol, o_idx: int, ring_i: int, out: set[int]) -> None:
    out.add(o_idx)
    for n2 in mol.GetAtomWithIdx(o_idx).GetNeighbors():
        if n2.GetAtomicNum() == 6 and n2.GetIdx() != ring_i:
            out.add(n2.GetIdx())


def _leaf_atoms_on(mol: Mol, ring: set[int], attach: int, parent: int) -> set[int]:
    out: set[int] = set()
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if nb.GetAtomicNum() == 8:
                _add_alkoxy_atoms(mol, nb.GetIdx(), i, out)
            else:
                out.add(nb.GetIdx())
    return out


def pyridinyl_atoms(mol: Mol, parent: set[int]) -> set[int]:
    out: set[int] = set()
    for p in ring_pyridinyls(mol, parent):
        out |= set(p["atoms"])
        out |= _leaf_atoms_on(mol, p["ring"], p["outer_c"], p["attach"])
    return out
