"""Aryl arms: Ph / OPh / CH2Ph / OCH2Ph with recursive nested leaves (P-29.3).

Ph may carry ≤3 leaves (simple + nested Ph/OPh up to max depth).
Locants from attach=1 (lowest set).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_depth2 import (
    _d2_pref_parts,
    _d2_sites,
    _depth2_atoms_on,
    _join_pref,
)
from namepredict.layer2.aryl_recurse import (
    _recurse_leaf_kind,
    _recursive_ph_name,
)

_HALO = frozenset({9, 17, 35, 53})
_HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
_HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}
_MULT_EN = {1: "", 2: "di", 3: "tri"}
_MULT_ZH = {1: "", 2: "二", 3: "三"}


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


def _is_terminal_me_leaf(mol: Mol, c_idx: int, ring_c: int) -> bool:
    """Methyl leaf on Ph: C only bonded to ring carbon (+H)."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetIsAromatic():
        return False
    return all(n.GetIdx() == ring_c or n.GetAtomicNum() == 1 for n in atom.GetNeighbors())


def _side_leaf_kind(mol: Mol, nb, ring_i: int) -> str | None:
    if _is_terminal_halo(nb):
        return "halo"
    if nb.GetAtomicNum() == 6 and _is_terminal_me_leaf(mol, nb.GetIdx(), ring_i):
        return "me"
    return _recurse_leaf_kind(mol, nb, ring_i)


def _count_side_leaves(mol: Mol, ring: set[int], attach: int, parent: int) -> int | None:
    """Count terminal halo + methyl leaves on Ph (exclude parent link)."""
    n = 0
    for i in ring:
        for nb in _nb_outside(mol, i, ring):
            if i == attach and nb.GetIdx() == parent:
                continue
            if _side_leaf_kind(mol, nb, i) is None:
                return None
            n += 1
    return n


def _side_ok(mol: Mol, ring: set[int], attach: int, parent: int) -> bool:
    n = _count_side_leaves(mol, ring, attach, parent)
    return n is not None and n <= 3


def _phenyl_at(mol: Mol, attach: int, parent: int) -> set[int] | None:
    hits = _c6_rings_at(mol, attach)
    if len(hits) != 1 or parent in hits[0]:
        return None
    ring = hits[0]
    return ring if _side_ok(mol, ring, attach, parent) else None


def _halo_list(mol: Mol, ring: set[int]) -> list[tuple[int, int]]:
    return [
        (i, nb.GetAtomicNum())
        for i in ring
        for nb in _nb_outside(mol, i, ring)
        if _is_terminal_halo(nb)
    ]


def _me_sites(mol: Mol, ring: set[int]) -> list[int]:
    sites: list[int] = []
    for i in ring:
        for nb in _nb_outside(mol, i, ring):
            if nb.GetAtomicNum() == 6 and _is_terminal_me_leaf(mol, nb.GetIdx(), i):
                sites.append(i)
    return sites


def _halo_pair(mol: Mol, ring: set[int]) -> tuple[int, int] | None:
    found = _halo_list(mol, ring)
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


def _loc_tuple(order: list[int], sites: list[int]) -> tuple[int, ...]:
    return tuple(sorted(order.index(s) + 1 for s in sites if s in order))


def _sub_sites(mol: Mol, ring: set[int], attach: int | None = None) -> list[int]:
    return (
        [s for s, _ in _halo_list(mol, ring)]
        + _me_sites(mol, ring)
        + _d2_sites(mol, ring, _nb_outside, attach)
    )


def _ring_order(mol: Mol, ring: set[int], start: int) -> list[int]:
    nbrs = _ring_nbrs(mol, start, ring)
    if len(nbrs) != 2:
        return [start]
    a = _walk_ring(mol, ring, start, nbrs[0])
    b = _walk_ring(mol, ring, start, nbrs[1])
    sites = _sub_sites(mol, ring, start)
    ta, tb = _loc_tuple(a, sites), _loc_tuple(b, sites)
    if ta != tb:
        return a if ta < tb else b
    return a if a <= b else b


def _ph_locant(mol: Mol, ring: set[int], attach: int, site: int) -> int:
    order = _ring_order(mol, ring, attach)
    return order.index(site) + 1 if site in order else 1


def _same_z_prefix(items: list[tuple[int, int]]) -> tuple[str, str]:
    z, n = items[0][1], len(items)
    locs = ",".join(str(l) for l, _ in items)
    return (
        f"{locs}-{_MULT_EN.get(n, '')}{_HALO_EN[z]}",
        f"{locs}-{_MULT_ZH.get(n, '')}{_HALO_ZH[z]}",
    )


def _mixed_z_prefix(items: list[tuple[int, int]]) -> tuple[str, str]:
    by_en = sorted(items, key=lambda lz: _HALO_EN[lz[1]])
    en = "-".join(f"{l}-{_HALO_EN[z]}" for l, z in by_en)
    zh = "-".join(f"{l}-{_HALO_ZH[z]}" for l, z in by_en)
    return en, zh


def _halo_prefix(items: list[tuple[int, int]]) -> tuple[str, str]:
    if len({z for _, z in items}) == 1:
        return _same_z_prefix(items)
    return _mixed_z_prefix(items)


def _me_prefix(locs: list[int]) -> tuple[str, str]:
    if not locs:
        return "", ""
    n = len(locs)
    ls = ",".join(str(l) for l in sorted(locs))
    return (
        f"{ls}-{_MULT_EN.get(n, '')}methyl",
        f"{ls}-{_MULT_ZH.get(n, '')}甲基",
    )


def _halo_me_prefs(mol: Mol, ph: set[int], attach: int) -> tuple[str, str, str, str]:
    items = sorted((_ph_locant(mol, ph, attach, s), z) for s, z in _halo_list(mol, ph))
    pe, pz = _halo_prefix(items) if items else ("", "")
    me = [_ph_locant(mol, ph, attach, s) for s in _me_sites(mol, ph)]
    me_e, me_z = _me_prefix(me)
    return pe, pz, me_e, me_z


def _ph_leaf_prefs(mol: Mol, ph: set[int], attach: int) -> tuple[str, str]:
    pe, pz, me_e, me_z = _halo_me_prefs(mol, ph, attach)
    de, dz = _d2_pref_parts(mol, ph, attach, _ph_locant, _nb_outside)
    # amino, alkoxy, halo, hydroxy, methyl, nitro, phenyl, CF3
    en = _join_pref([de[0], de[1], pe, de[2], me_e, de[3], de[4], de[5]])
    zh = _join_pref([dz[0], dz[1], pz, dz[2], me_z, dz[3], dz[4], dz[5]])
    return en, zh


def _ph_bridge_attach(mol: Mol, ph: set[int], bridge: int) -> int:
    for i in ph:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetIdx() == bridge:
                return i
    return min(ph)


def _stem_swap(en: str, zh: str, from_en: str, to_en: str, from_zh: str, to_zh: str):
    if en.endswith(from_en):
        en = en[: -len(from_en)] + to_en
    if zh.endswith(from_zh):
        zh = zh[: -len(from_zh)] + to_zh
    return en, zh


def _phenyl_name(mol: Mol, ph: set[int], attach: int) -> tuple[str, str, bool]:
    # parent of arm Ph is outside attach; use any non-ring heavy as parent proxy
    parent = _arm_parent_of(mol, ph, attach)
    en, zh, paren, _ = _recursive_ph_name(mol, ph, attach, parent)
    return en, zh, paren


def _arm_parent_of(mol: Mol, ph: set[int], attach: int) -> int:
    for n in mol.GetAtomWithIdx(attach).GetNeighbors():
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ph:
            return n.GetIdx()
    return -1


def _phenoxy_name(mol: Mol, ph: set[int], outer: int) -> tuple[str, str, bool]:
    parent = _arm_parent_of(mol, ph, outer)
    en, zh, paren, _ = _recursive_ph_name(mol, ph, outer, parent)
    en, zh = _stem_swap(en, zh, "phenyl", "phenoxy", "苯基", "苯氧基")
    return en, zh, paren


def _benzyl_name(mol: Mol, ph: set[int], ch2: int) -> tuple[str, str, bool]:
    att = _ph_bridge_attach(mol, ph, ch2)
    en, zh, paren, _ = _recursive_ph_name(mol, ph, att, ch2)
    en, zh = _stem_swap(en, zh, "phenyl", "benzyl", "苯基", "苄基")
    return en, zh, paren


def _benzyloxy_name(mol: Mol, ph: set[int], ch2: int) -> tuple[str, str, bool]:
    att = _ph_bridge_attach(mol, ph, ch2)
    en, zh, paren, _ = _recursive_ph_name(mol, ph, att, ch2)
    en, zh = _stem_swap(en, zh, "phenyl", "benzyloxy", "苯基", "苄氧基")
    return en, zh, paren


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _is_open_ch2(atom) -> bool:
    return atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing()


def _other_heavy(atom, parent: int) -> int | None:
    heavies = _heavies(atom)
    if len(heavies) != 2:
        return None
    ids = {h.GetIdx() for h in heavies}
    return None if parent not in ids else (ids - {parent}).pop()


def _ch2_ph_at(mol: Mol, ch2: int, parent: int) -> set[int] | None:
    atom = mol.GetAtomWithIdx(ch2)
    if not _is_open_ch2(atom):
        return None
    other = _other_heavy(atom, parent)
    return None if other is None else _phenyl_at(mol, other, ch2)


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


def _benzyloxy_from_ether(mol: Mol, e: dict, parent: set[int]) -> dict | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in parent) == (c2 in parent):
        return None
    ring_c, outer = (c1, c2) if c1 in parent else (c2, c1)
    ph = _ch2_ph_at(mol, outer, o)
    if ph is None:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "ch2": outer,
            "ph": ph, "atoms": [o, outer, *ph]}


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


def _benzyl_from_start(mol: Mol, start: int, parent: set[int]) -> dict | None:
    if start in parent:
        return None
    att = _parent_link(mol, start, parent)
    if att is None:
        return None
    ph = _ch2_ph_at(mol, start, att)
    if ph is None:
        return None
    return {"attach": att, "outer_c": start, "ch2": start, "ph": ph,
            "atoms": [start, *ph]}


def _ring_phenoxys(info: dict, parent: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _phenoxy_from_ether(mol, e, parent)
        if one is not None:
            out.append(one)
    return out


def _ring_benzyloxys(info: dict, parent: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _benzyloxy_from_ether(mol, e, parent)
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


def _ring_benzyls(mol: Mol, parent: set[int]) -> list[dict]:
    out: list[dict] = []
    for s in _side_c_starts(mol, parent):
        one = _benzyl_from_start(mol, s, parent)
        if one is not None:
            out.append(one)
    return out


def _simple_leaf_atoms(mol: Mol, ph: set[int]) -> set[int]:
    out: set[int] = set()
    for i in ph:
        for nb in _nb_outside(mol, i, ph):
            if _is_terminal_halo(nb):
                out.add(nb.GetIdx())
            elif nb.GetAtomicNum() == 6 and _is_terminal_me_leaf(mol, nb.GetIdx(), i):
                out.add(nb.GetIdx())
    return out


def _halo_atoms_on(mol: Mol, ph: set[int], attach: int | None = None) -> set[int]:
    """Halo + methyl + depth-2 leaf atom indices on Ph (legacy name)."""
    return _simple_leaf_atoms(mol, ph) | _depth2_atoms_on(mol, ph, _nb_outside, attach)


def _arm_attach_parent(mol: Mol, p: dict) -> tuple[int, int]:
    if "ch2" in p:
        return _ph_bridge_attach(mol, p["ph"], p["ch2"]), p["ch2"]
    att = p.get("outer_c", p.get("attach"))
    return att, _arm_parent_of(mol, p["ph"], att)


def _one_arm_atoms(mol: Mol, p: dict) -> set[int]:
    att, parent = _arm_attach_parent(mol, p)
    _, _, _, nested = _recursive_ph_name(mol, p["ph"], att, parent)
    return set(p["atoms"]) | nested


def _aryl_atoms(info: dict, parent: set[int]) -> set[int]:
    mol: Mol = info["mol"]
    out: set[int] = set()
    arms = (
        _ring_phenoxys(info, parent) + _ring_benzyloxys(info, parent)
        + _ring_phenyls(mol, parent) + _ring_benzyls(mol, parent)
    )
    for p in arms:
        out |= _one_arm_atoms(mol, p)
    return out


def _aryl_sub_n(info: dict, parent: set[int]) -> int:
    mol: Mol = info["mol"]
    return (
        len(_ring_phenyls(mol, parent)) + len(_ring_phenoxys(info, parent))
        + len(_ring_benzyls(mol, parent)) + len(_ring_benzyloxys(info, parent))
    )


def _arom_c6_rings(mol: Mol) -> list[set[int]]:
    return [set(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]


def _is_unfused_benzene_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _phenyl_starts_set(mol: Mol, parent: set[int]) -> set[int]:
    ps = {p["outer_c"] for p in _ring_phenyls(mol, parent)}
    return ps | {p["outer_c"] for p in _ring_benzyls(mol, parent)}


def _unsub_phenyl_at(mol: Mol, c_idx: int, n_idx: int) -> bool:
    """True if c_idx is sole N-attachment of unsubstituted phenyl (amide N-Ph)."""
    if mol.GetAtomWithIdx(n_idx).GetAtomicNum() != 7:
        return False
    ph = _phenyl_at(mol, c_idx, n_idx)
    return ph is not None and not _halo_list(mol, ph)


def _arom_c6_ring_lists(mol: Mol) -> list[list[int]]:
    return [list(r) for r in mol.GetRingInfo().AtomRings() if _is_arom_c6(mol, r)]


def _unfused_c6_at(mol: Mol, c_idx: int) -> set[int] | None:
    """Sole unfused aromatic C6 containing c_idx, else None."""
    hits = [
        set(r) for r in _arom_c6_ring_lists(mol)
        if c_idx in r and _is_unfused_benzene_ring(mol, set(r))
    ]
    return hits[0] if len(hits) == 1 else None


def _exocyclic_fg_ring(mol: Mol, fg_c: int) -> set[int] | None:
    """Unfused C6 attached to exocyclic FG carbon (COOH/CHO/CN/...)."""
    nbs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(fg_c).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIsAromatic()
    ]
    return _unfused_c6_at(mol, nbs[0]) if len(nbs) == 1 else None


def _aryl_exclude(info: dict, parent: set[int]) -> set[int]:
    mol: Mol = info["mol"]
    return _aryl_atoms(info, parent) | _phenyl_starts_set(mol, parent)
