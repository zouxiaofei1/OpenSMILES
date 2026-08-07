"""Ph ring namer: number + assemble leaves via LeafRegistry."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.tools.side_facts import (
    ArylLeafFact, ArylLeafKind, ring_leaf, ring_outside,
)

Match = ArylLeafFact


def _join_pref(parts: list[str]) -> str:
    return "-".join(part for part in parts if part)


def nb_out(mol: Mol, atom: int, ring: set[int]) -> list:
    return ring_outside(mol, atom, ring)


def match_leaf(mol: Mol, atom: object, ring_atom: int, depth: int) -> Match | None:
    return ring_leaf(mol, atom, ring_atom, depth)


_SIMPLE_NAMES = {
    ArylLeafKind.AMINO: ("amino", "氨基"),
    ArylLeafKind.CYANO: ("cyano", "氰基"),
    ArylLeafKind.HYDROXY: ("hydroxy", "羟基"),
    ArylLeafKind.METHYL: ("methyl", "甲基"),
    ArylLeafKind.METHYLSULFANYL: ("methylsulfanyl", "甲硫基"),
    ArylLeafKind.NITRO: ("nitro", "硝基"),
    ArylLeafKind.TRIFLUOROMETHYL: ("(trifluoromethyl)", "三氟甲基"),
}
_ALKOXY = {1: ("methoxy", "甲氧基"), 2: ("ethoxy", "乙氧基"),
            3: ("propoxy", "丙氧基"), 4: ("butoxy", "丁氧基")}
_ALKYL = {2: ("ethyl", "乙基"), 3: ("propyl", "丙基"), 4: ("butyl", "丁基")}
_COMPLEX = {ArylLeafKind.PHENYL: "phenyl", ArylLeafKind.PHENOXY: "phenoxy",
            ArylLeafKind.BENZYL: "benzyl"}


def _simple_name(fact: Match) -> tuple[str, str]:
    if fact.kind == ArylLeafKind.HALOGEN:
        return _HALO_EN[fact.value], _HALO_ZH[fact.value]
    if fact.kind == ArylLeafKind.ALKOXY:
        return _ALKOXY[fact.value]
    if fact.kind == ArylLeafKind.N_ALKYL:
        return _ALKYL[fact.value]
    return _SIMPLE_NAMES[fact.kind]


def _complex_name(mol: Mol, fact: Match, depth: int) -> tuple[str, str, set[int]]:
    en, zh, atoms = name_ph_ring(mol, set(fact.child_ring), fact.child_attach,
                                 fact.child_parent, depth + 1)
    stem = _COMPLEX[fact.kind]
    en, zh = en[:-6] + stem, zh[:-2] + {"phenyl": "苯基", "phenoxy": "苯氧基",
                                      "benzyl": "苄基"}[stem]
    return en, zh, atoms | ({fact.extra_atom} if fact.extra_atom >= 0 else set())


def name_leaf(mol: Mol, fact: Match, depth: int) -> tuple[str, str, set[int]]:
    if fact.kind in _COMPLEX:
        return _complex_name(mol, fact, depth)
    en, zh = _simple_name(fact)
    return en, zh, set(fact.atoms)

from namepredict.constants import HALO_EN as _HALO_EN, HALO_ZH as _HALO_ZH, MULT_EN as _MULT_EN, MULT_ZH as _MULT_ZH


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
) -> list[Match]:
    out: list[Match] = []
    for i in ring:
        for nb in nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            got = match_leaf(mol, nb, i, depth)
            if got is not None:
                out.append(got)
    return out


def _sites(leaves: list[Match]) -> list[int]:
    return [fact.site for fact in leaves]


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
    mol: Mol, leaves: list[Match], order: list[int], depth: int,
) -> list[tuple[int, str, str, str, bool, set[int]]]:
    out = []
    for fact in leaves:
        loc = _loc(order, fact.site)
        en, zh, atoms = name_leaf(mol, fact, depth)
        complex_leaf = fact.kind in _COMPLEX
        out.append((loc, en, zh, en, complex_leaf, atoms))
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
