"""Ph ring namer: number + assemble leaves via LeafRegistry."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import _join_pref
from namepredict.layer2.leaves.protocol import Match
from namepredict.layer2.leaves.registry import match_leaf, name_leaf
from namepredict.layer2.leaves.topo import nb_out

_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}
_MULT_EN = {1: "", 2: "di", 3: "tri"}
_MULT_ZH = {1: "", 2: "二", 3: "三"}


def _ring_nbrs(mol: Mol, i: int, ring: set[int]) -> list[int]:
    return sorted(
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in ring
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


def _collect_leaves(
    mol: Mol, ring: set[int], attach: int, parent: int, depth: int,
) -> list[tuple[Match, object]]:
    out: list[tuple[Match, object]] = []
    for i in ring:
        for nb in nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            got = match_leaf(mol, nb, i, depth)
            if got is not None:
                out.append(got)
    return out


def _sites(leaves: list[tuple[Match, object]]) -> list[int]:
    return [m["site"] for m, _ in leaves]


def _loc_tuple(order: list[int], sites: list[int]) -> tuple[int, ...]:
    return tuple(sorted(order.index(s) + 1 for s in sites if s in order))


def _ring_order(mol: Mol, ring: set[int], attach: int, sites: list[int]) -> list[int]:
    nbrs = _ring_nbrs(mol, attach, ring)
    if len(nbrs) != 2:
        return [attach]
    a = _walk_ring(mol, ring, attach, nbrs[0])
    b = _walk_ring(mol, ring, attach, nbrs[1])
    ta, tb = _loc_tuple(a, sites), _loc_tuple(b, sites)
    if ta != tb:
        return a if ta < tb else b
    return a if a <= b else b


def _loc(order: list[int], site: int) -> int:
    return order.index(site) + 1 if site in order else 1


def _wrap(stem: str) -> str:
    if not stem or stem in ("phenyl", "phenoxy"):
        return stem
    if stem[0].isdigit() or "(" in stem:
        return f"[{stem}]" if "(" in stem else f"({stem})"
    return stem


def _mult_pref(locs: list[int], en: str, zh: str) -> tuple[str, str]:
    n, ls = len(locs), ",".join(str(l) for l in locs)
    return (
        f"{ls}-{_MULT_EN.get(n, '')}{en}",
        f"{ls}-{_MULT_ZH.get(n, '')}{zh}",
    )


def _name_all(
    mol: Mol, leaves: list[tuple[Match, object]], order: list[int], depth: int,
) -> list[tuple[int, str, str, str, bool, set[int]]]:
    out = []
    for m, h in leaves:
        loc = _loc(order, m["site"])
        en, zh, at = name_leaf(mol, m, h, depth)
        out.append((loc, en, zh, h.alpha_key, h.complex, at))
    return out


def _stem_key(en: str) -> str:
    return en[1:] if en.startswith("(") else en


def _bucket_named(named: list) -> tuple[dict, list]:
    from collections import defaultdict

    simple: dict[str, list] = defaultdict(list)
    complex_items: list = []
    for item in named:
        if item[4]:
            complex_items.append(item)
        else:
            simple[item[1]].append(item)
    return simple, complex_items


def _alpha_merge(named: list) -> tuple[list[str], list[str]]:
    """Assemble prefix parts: simple grouped by stem en; complex by stem."""
    return _emit_alpha(*_bucket_named(named))


def _emit_simple(simple: dict) -> tuple[list[str], list[str]]:
    parts_e, parts_z = [], []
    for en in sorted(simple, key=_stem_key):
        items = simple[en]
        locs = sorted(l for l, *_ in items)
        pe, pz = _mult_pref(locs, items[0][1], items[0][2])
        parts_e.append(pe)
        parts_z.append(pz)
    return parts_e, parts_z


def _emit_complex(complex_items: list) -> tuple[list[str], list[str]]:
    parts_e, parts_z = [], []
    for loc, en, zh, _, _, _ in sorted(
        complex_items, key=lambda x: (_stem_key(x[1]), x[0]),
    ):
        parts_e.append(f"{loc}-{_wrap(en)}")
        parts_z.append(f"{loc}-{_wrap(zh)}")
    return parts_e, parts_z


def _emit_alpha(simple: dict, complex_items: list) -> tuple[list[str], list[str]]:
    se, sz = _emit_simple(simple)
    ce, cz = _emit_complex(complex_items)
    return se + ce, sz + cz


def name_ph_ring(
    mol: Mol, ph: set[int], attach: int, parent: int, depth: int,
) -> tuple[str, str, set[int]]:
    """Recursive Ph name via registry; returns en, zh, atoms."""
    leaves = _collect_leaves(mol, ph, attach, parent, depth)
    order = _ring_order(mol, ph, attach, _sites(leaves))
    named = _name_all(mol, leaves, order, depth)
    pe, pz = _alpha_merge(named)
    atoms = set().union(*(at for *_, at in named)) if named else set()
    return _join_pref(pe) + "phenyl", _join_pref(pz) + "苯基", set(ph) | atoms


def recursive_ph_name(
    mol: Mol, ph: set[int], attach: int, parent: int,
) -> tuple[str, str, bool, set[int]]:
    en, zh, atoms = name_ph_ring(mol, ph, attach, parent, depth=1)
    return en, zh, en != "phenyl", atoms
