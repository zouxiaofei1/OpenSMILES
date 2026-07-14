"""Recursive aryl substituent naming (Ph / OPh arms, depth can exceed 2).

Names a phenyl ring with simple leaves and nested Ph/OPh, each nested
ring named the same way until max_depth. Topology gates + prefix assembly.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import (
    _ALKOXY_EN,
    _ALKOXY_ZH,
    _alkoxy_n,
    _depth2_hetero,
    _heavies,
    _is_arom_c6_set,
    _is_cf3_leaf,
    _is_terminal_nh2,
    _is_terminal_nitro,
    _is_terminal_oh,
    _join_pref,
    _leaf_atoms_alkoxy,
    _leaf_atoms_cf3,
    _leaf_atoms_nitro,
    _nested_c6_at,
    _same_stem_prefix,
)
from namepredict.layer2.side_alkyl import _is_cf3_carbon

_MAX_DEPTH = 3
_HALO = frozenset({9, 17, 35, 53})
_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}
_MULT_EN = {1: "", 2: "di", 3: "tri"}
_MULT_ZH = {1: "", 2: "二", 3: "三"}


def _nb_out(mol: Mol, i: int, ring: set[int]) -> list:
    return [
        n for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ring
    ]


def _is_halo(nb) -> bool:
    if nb.GetAtomicNum() not in _HALO:
        return False
    return sum(1 for n in nb.GetNeighbors() if n.GetAtomicNum() != 1) == 1


def _is_me_leaf(mol: Mol, c_idx: int, ring_c: int) -> bool:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetIsAromatic():
        return False
    return all(n.GetIdx() == ring_c or n.GetAtomicNum() == 1 for n in atom.GetNeighbors())


def _simple_kind(mol: Mol, nb, ring_i: int) -> str | None:
    if _is_halo(nb):
        return "halo"
    if nb.GetAtomicNum() == 6 and _is_me_leaf(mol, nb.GetIdx(), ring_i):
        return "me"
    if nb.GetAtomicNum() == 6 and _is_cf3_leaf(mol, nb.GetIdx(), ring_i):
        return "cf3"
    return _depth2_hetero(mol, nb, ring_i)


def _ether_outer_c(o_atom, ring_i: int) -> int | None:
    if o_atom.GetAtomicNum() != 8 or o_atom.IsInRing():
        return None
    nbs = _heavies(o_atom)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def _phenoxy_outer(mol: Mol, o_atom, ring_i: int) -> int | None:
    """Ether O on ring → outer aromatic C of OPh, else None."""
    outer = _ether_outer_c(o_atom, ring_i)
    if outer is None:
        return None
    return outer if mol.GetAtomWithIdx(outer).GetIsAromatic() else None


def _nested_ph_ring(mol: Mol, nb, ring_i: int) -> set[int] | None:
    if nb.GetAtomicNum() != 6 or not nb.GetIsAromatic():
        return None
    ph = _nested_c6_at(mol, nb.GetIdx(), {ring_i})
    return None if ph is None or ring_i in ph else ph


def _ph_leaf_count(mol: Mol, ph: set[int], attach: int, parent: int, depth: int) -> int | None:
    """Count leaves if all valid, else None."""
    n = 0
    for i in ph:
        for nb in _nb_out(mol, i, ph):
            if i == attach and nb.GetIdx() == parent:
                continue
            if not _leaf_ok(mol, nb, i, depth):
                return None
            n += 1
    return n


def _ph_ok(mol: Mol, ph: set[int], attach: int, parent: int, depth: int) -> bool:
    """All outsides of ph (except parent link) are valid leaves; ≤3 leaves."""
    if depth > _MAX_DEPTH:
        return False
    n = _ph_leaf_count(mol, ph, attach, parent, depth)
    return n is not None and n <= 3


def _leaf_ok(mol: Mol, nb, ring_i: int, depth: int) -> bool:
    if _simple_kind(mol, nb, ring_i) is not None:
        return True
    ph = _nested_ph_ring(mol, nb, ring_i)
    if ph is not None:
        return _ph_ok(mol, ph, nb.GetIdx(), ring_i, depth + 1)
    outer = _phenoxy_outer(mol, nb, ring_i)
    if outer is None:
        return False
    ph2 = _nested_c6_at(mol, outer, {ring_i, nb.GetIdx()})
    return ph2 is not None and _ph_ok(mol, ph2, outer, nb.GetIdx(), depth + 1)


def _is_nested_ph_leaf(mol: Mol, nb, ring_i: int, depth: int = 1) -> bool:
    ph = _nested_ph_ring(mol, nb, ring_i)
    return ph is not None and _ph_ok(mol, ph, nb.GetIdx(), ring_i, depth + 1)


def _is_nested_oph_leaf(mol: Mol, nb, ring_i: int, depth: int = 1) -> bool:
    outer = _phenoxy_outer(mol, nb, ring_i)
    if outer is None:
        return False
    ph = _nested_c6_at(mol, outer, {ring_i, nb.GetIdx()})
    return ph is not None and _ph_ok(mol, ph, outer, nb.GetIdx(), depth + 1)


def _recurse_leaf_kind(mol: Mol, nb, ring_i: int) -> str | None:
    """Leaf kind for side_ok / counting (depth-1 gate uses depth=1)."""
    sk = _simple_kind(mol, nb, ring_i)
    if sk is not None:
        return sk
    if _is_nested_ph_leaf(mol, nb, ring_i, 1):
        return "phenyl"
    if _is_nested_oph_leaf(mol, nb, ring_i, 1):
        return "phenoxy"
    return None


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


def _loc_tuple(order: list[int], sites: list[int]) -> tuple[int, ...]:
    return tuple(sorted(order.index(s) + 1 for s in sites if s in order))


def _sub_sites_ring(mol: Mol, ring: set[int], attach: int, parent: int) -> list[int]:
    sites: list[int] = []
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if _recurse_leaf_kind(mol, nb, i) is not None:
                sites.append(i)
    return sites


def _ring_order(mol: Mol, ring: set[int], attach: int, parent: int) -> list[int]:
    nbrs = _ring_nbrs(mol, attach, ring)
    if len(nbrs) != 2:
        return [attach]
    a = _walk_ring(mol, ring, attach, nbrs[0])
    b = _walk_ring(mol, ring, attach, nbrs[1])
    sites = _sub_sites_ring(mol, ring, attach, parent)
    ta, tb = _loc_tuple(a, sites), _loc_tuple(b, sites)
    if ta != tb:
        return a if ta < tb else b
    return a if a <= b else b


def _loc(mol: Mol, ring: set[int], attach: int, parent: int, site: int) -> int:
    order = _ring_order(mol, ring, attach, parent)
    return order.index(site) + 1 if site in order else 1


def _halo_items(mol: Mol, ring: set[int], attach: int, parent: int) -> list[tuple[int, int]]:
    return sorted(
        (_loc(mol, ring, attach, parent, i), nb.GetAtomicNum())
        for i in ring
        for nb in _nb_out(mol, i, ring)
        if not (i == attach and nb.GetIdx() == parent) and _is_halo(nb)
    )


def _me_locs(mol: Mol, ring: set[int], attach: int, parent: int) -> list[int]:
    return [
        _loc(mol, ring, attach, parent, i)
        for i in ring
        for nb in _nb_out(mol, i, ring)
        if not (i == attach and nb.GetIdx() == parent)
        and nb.GetAtomicNum() == 6 and _is_me_leaf(mol, nb.GetIdx(), i)
    ]


def _same_z_pref(items: list[tuple[int, int]]) -> tuple[str, str]:
    z, n = items[0][1], len(items)
    locs = ",".join(str(l) for l, _ in items)
    return (
        f"{locs}-{_MULT_EN.get(n, '')}{_HALO_EN[z]}",
        f"{locs}-{_MULT_ZH.get(n, '')}{_HALO_ZH[z]}",
    )


def _mixed_z_pref(items: list[tuple[int, int]]) -> tuple[str, str]:
    by_en = sorted(items, key=lambda lz: _HALO_EN[lz[1]])
    return (
        "-".join(f"{l}-{_HALO_EN[z]}" for l, z in by_en),
        "-".join(f"{l}-{_HALO_ZH[z]}" for l, z in by_en),
    )


def _halo_pref(items: list[tuple[int, int]]) -> tuple[str, str]:
    if not items:
        return "", ""
    if len({z for _, z in items}) == 1:
        return _same_z_pref(items)
    return _mixed_z_pref(items)


def _me_pref(locs: list[int]) -> tuple[str, str]:
    if not locs:
        return "", ""
    n, ls = len(locs), ",".join(str(l) for l in sorted(locs))
    return f"{ls}-{_MULT_EN.get(n, '')}methyl", f"{ls}-{_MULT_ZH.get(n, '')}甲基"


def _wrap_complex(stem: str) -> str:
    """Paren compound nested leaf only (not bare phenyl/phenoxy)."""
    if not stem or stem in ("phenyl", "phenoxy"):
        return stem
    if stem[0].isdigit() or "(" in stem:
        return f"[{stem}]" if "(" in stem else f"({stem})"
    return stem


def _name_nested_ph(mol: Mol, ph: set[int], attach: int, parent: int, depth: int) -> tuple[str, str, set[int]]:
    en, zh, atoms = _name_ph(mol, ph, attach, parent, depth)
    return en, zh, atoms


def _to_phenoxy_stem(en: str, zh: str) -> tuple[str, str]:
    if en.endswith("phenyl"):
        en = en[: -len("phenyl")] + "phenoxy"
    if zh.endswith("苯基"):
        zh = zh[: -len("苯基")] + "苯氧基"
    return en, zh


def _name_nested_oph(mol: Mol, o_idx: int, ring_i: int, depth: int) -> tuple[str, str, set[int]]:
    outer = _phenoxy_outer(mol, mol.GetAtomWithIdx(o_idx), ring_i)
    if outer is None:
        return "phenoxy", "苯氧基", {o_idx}
    ph = _nested_c6_at(mol, outer, {ring_i, o_idx})
    if ph is None:
        return "phenoxy", "苯氧基", {o_idx}
    en, zh, atoms = _name_ph(mol, ph, outer, o_idx, depth)
    en, zh = _to_phenoxy_stem(en, zh)
    return en, zh, atoms | {o_idx}


def _collect_complex_items(
    mol: Mol, ring: set[int], attach: int, parent: int, depth: int,
) -> list[tuple[int, str, str, set[int]]]:
    items: list[tuple[int, str, str, set[int]]] = []
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            got = _try_complex_leaf(mol, nb, i, depth)
            if got is not None:
                loc = _loc(mol, ring, attach, parent, i)
                items.append((loc, got[0], got[1], got[2]))
    return items


def _complex_leaf_prefs(
    mol: Mol, ring: set[int], attach: int, parent: int, depth: int,
) -> tuple[list[str], list[str], set[int]]:
    """Named nested Ph/OPh leaves with locants; returns en/zh parts + atoms."""
    parts_e, parts_z, atoms = [], [], set()
    items = _collect_complex_items(mol, ring, attach, parent, depth)
    for loc, en, zh, at in sorted(items, key=lambda x: (x[1], x[0])):
        parts_e.append(f"{loc}-{_wrap_complex(en)}")
        parts_z.append(f"{loc}-{_wrap_complex(zh)}")
        atoms |= at
    return parts_e, parts_z, atoms


def _try_complex_leaf(mol: Mol, nb, ring_i: int, depth: int) -> tuple[str, str, set[int]] | None:
    ph = _nested_ph_ring(mol, nb, ring_i)
    if ph is not None and _ph_ok(mol, ph, nb.GetIdx(), ring_i, depth + 1):
        return _name_nested_ph(mol, ph, nb.GetIdx(), ring_i, depth + 1)
    if _is_nested_oph_leaf(mol, nb, ring_i, depth):
        return _name_nested_oph(mol, nb.GetIdx(), ring_i, depth + 1)
    return None


def _simple_pref_parts(
    mol: Mol, ring: set[int], attach: int, parent: int,
) -> tuple[list[str], list[str], set[int]]:
    """Halo/Me/alkoxy/nitro/OH/NH2/CF3 prefixes + atoms."""
    pe, pz = _halo_pref(_halo_items(mol, ring, attach, parent))
    me_e, me_z = _me_pref(_me_locs(mol, ring, attach, parent))
    de, dz, atoms = _d2_simple_prefs(mol, ring, attach, parent)
    en = [x for x in [de[0], pe, de[1], me_e, de[2], de[3], de[4]] if x]
    zh = [x for x in [dz[0], pz, dz[1], me_z, dz[2], dz[3], dz[4]] if x]
    atoms |= _simple_atoms(mol, ring, attach, parent)
    return en, zh, atoms


def _acc_hetero(mol, nb, i, loc, buckets, atoms) -> None:
    am, oh, no2, _cf3, alk = buckets
    k = _depth2_hetero(mol, nb, i)
    if k == "amino":
        am.append(loc); atoms.add(nb.GetIdx())
    elif k == "hydroxy":
        oh.append(loc); atoms.add(nb.GetIdx())
    elif k == "nitro":
        no2.append(loc); atoms |= _leaf_atoms_nitro(nb)
    elif k in _ALKOXY_EN.values() and (n := _alkoxy_n(mol, nb, i)):
        alk.append((loc, n)); atoms |= _leaf_atoms_alkoxy(mol, nb.GetIdx(), i)


def _acc_simple(mol, nb, i, loc, buckets, atoms) -> None:
    _acc_hetero(mol, nb, i, loc, buckets, atoms)
    if nb.GetAtomicNum() == 6 and _is_cf3_leaf(mol, nb.GetIdx(), i):
        buckets[3].append(loc)
        atoms |= _leaf_atoms_cf3(mol, nb.GetIdx())


def _d2_simple_prefs(
    mol: Mol, ring: set[int], attach: int, parent: int,
) -> tuple[list[str], list[str], set[int]]:
    buckets: tuple = ([], [], [], [], [])
    atoms: set[int] = set()
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            _acc_simple(mol, nb, i, _loc(mol, ring, attach, parent, i), buckets, atoms)
    am, oh, no2, cf3, alk = buckets
    return _pack_simple(am, alk, oh, no2, cf3), _pack_simple_zh(am, alk, oh, no2, cf3), atoms


def _pack_simple(am, alk, oh, no2, cf3) -> list[str]:
    from namepredict.layer2.aryl_depth2 import (
        _alkoxy_prefix, _amino_prefix, _cf3_prefix, _hydroxy_prefix, _nitro_prefix,
    )
    ae, _ = _amino_prefix(am)
    oe, _ = _alkoxy_prefix(alk)
    he, _ = _hydroxy_prefix(oh)
    ne, _ = _nitro_prefix(no2)
    ce, _ = _cf3_prefix(cf3)
    return [ae, oe, he, ne, ce]


def _pack_simple_zh(am, alk, oh, no2, cf3) -> list[str]:
    from namepredict.layer2.aryl_depth2 import (
        _alkoxy_prefix, _amino_prefix, _cf3_prefix, _hydroxy_prefix, _nitro_prefix,
    )
    _, az = _amino_prefix(am)
    _, oz = _alkoxy_prefix(alk)
    _, hz = _hydroxy_prefix(oh)
    _, nz = _nitro_prefix(no2)
    _, cz = _cf3_prefix(cf3)
    return [az, oz, hz, nz, cz]


def _simple_atoms(mol: Mol, ring: set[int], attach: int, parent: int) -> set[int]:
    out: set[int] = set()
    for i in ring:
        for nb in _nb_out(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if _is_halo(nb) or (
                nb.GetAtomicNum() == 6 and _is_me_leaf(mol, nb.GetIdx(), i)
            ):
                out.add(nb.GetIdx())
    return out


def _name_ph(
    mol: Mol, ph: set[int], attach: int, parent: int, depth: int,
) -> tuple[str, str, set[int]]:
    """Recursive Ph name: prefixes + phenyl; atoms include leaves + nested."""
    se, sz, sa = _simple_pref_parts(mol, ph, attach, parent)
    ce, cz, ca = _complex_leaf_prefs(mol, ph, attach, parent, depth)
    # alpha-ish merge: simple parts already ordered; complex by stem then loc
    en = _join_pref(se + ce) + "phenyl"
    zh = _join_pref(sz + cz) + "苯基"
    return en, zh, set(ph) | sa | ca


def _recursive_ph_name(
    mol: Mol, ph: set[int], attach: int, parent: int,
) -> tuple[str, str, bool, set[int]]:
    """Public: name arm Ph with recursion. paren if any leaf present."""
    en, zh, atoms = _name_ph(mol, ph, attach, parent, depth=1)
    has_leaf = en != "phenyl"
    return en, zh, has_leaf, atoms
