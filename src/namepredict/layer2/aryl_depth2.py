"""Depth-2 simple leaves on Ph arms: alkoxy / nitro / CF3 (P-29.3).

Topology + prefix only; ring walk / attach locants stay in aryl_sub.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_alkoxy import _linear_alkyl_atoms
from namepredict.layer2.side_alkyl import _is_cf3_carbon

_ALKOXY_EN = {1: "methoxy", 2: "ethoxy"}
_ALKOXY_ZH = {1: "甲氧基", 2: "乙氧基"}
_MULT_EN = {1: "", 2: "di", 3: "tri"}
_MULT_ZH = {1: "", 2: "二", 3: "三"}


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _nitro_os_ok(nbs: list) -> bool:
    os = [n for n in nbs if n.GetAtomicNum() == 8]
    cs = [n for n in nbs if n.GetAtomicNum() == 6]
    if len(os) != 2 or len(cs) != 1:
        return False
    return all(len(_heavies(o)) == 1 for o in os)


def _is_terminal_nitro(atom) -> bool:
    """Ring-NO2: N with exactly two O neighbors (no further heavy)."""
    if atom.GetAtomicNum() != 7:
        return False
    nbs = _heavies(atom)
    return len(nbs) == 3 and _nitro_os_ok(nbs)


def _ether_outer(o_atom, ring_i: int) -> int | None:
    if o_atom.GetAtomicNum() != 8 or o_atom.IsInRing():
        return None
    nbs = _heavies(o_atom)
    if len(nbs) != 2:
        return None
    ids = {n.GetIdx() for n in nbs}
    return None if ring_i not in ids else (ids - {ring_i}).pop()


def _alkoxy_n(mol: Mol, o_atom, ring_i: int) -> int | None:
    """Ether O on ring → linear C1–C2 alkyl length, else None."""
    outer = _ether_outer(o_atom, ring_i)
    if outer is None:
        return None
    path = _linear_alkyl_atoms(mol, outer, o_atom.GetIdx())
    return len(path) if path is not None and len(path) in (1, 2) else None


def _is_cf3_leaf(mol: Mol, c_idx: int, ring_i: int) -> bool:
    if not _is_cf3_carbon(mol, c_idx):
        return False
    atom = mol.GetAtomWithIdx(c_idx)
    cs = [n for n in _heavies(atom) if n.GetAtomicNum() == 6]
    return len(cs) == 1 and cs[0].GetIdx() == ring_i


def _alkoxy_kind(n: int | None) -> str | None:
    return {1: "methoxy", 2: "ethoxy"}.get(n) if n else None


def _depth2_kind(mol: Mol, nb, ring_i: int) -> str | None:
    """Return 'nitro' | 'cf3' | 'methoxy' | 'ethoxy' | None for one outside nb."""
    if _is_terminal_nitro(nb):
        return "nitro"
    if nb.GetAtomicNum() == 6 and _is_cf3_leaf(mol, nb.GetIdx(), ring_i):
        return "cf3"
    return _alkoxy_kind(_alkoxy_n(mol, nb, ring_i))


def _nitro_sites(mol: Mol, ring: set[int], nb_outside) -> list[int]:
    return [
        i for i in ring
        for nb in nb_outside(mol, i, ring)
        if _is_terminal_nitro(nb)
    ]


def _cf3_sites(mol: Mol, ring: set[int], nb_outside) -> list[int]:
    return [
        i for i in ring
        for nb in nb_outside(mol, i, ring)
        if nb.GetAtomicNum() == 6 and _is_cf3_leaf(mol, nb.GetIdx(), i)
    ]


def _alkoxy_sites(mol: Mol, ring: set[int], nb_outside) -> list[tuple[int, int]]:
    """List of (ring_site, alkoxy_n) for methoxy/ethoxy leaves."""
    out: list[tuple[int, int]] = []
    for i in ring:
        for nb in nb_outside(mol, i, ring):
            n = _alkoxy_n(mol, nb, i)
            if n is not None:
                out.append((i, n))
    return out


def _same_stem_prefix(locs: list[int], en: str, zh: str) -> tuple[str, str]:
    if not locs:
        return "", ""
    n = len(locs)
    ls = ",".join(str(l) for l in sorted(locs))
    return (
        f"{ls}-{_MULT_EN.get(n, '')}{en}",
        f"{ls}-{_MULT_ZH.get(n, '')}{zh}",
    )


def _cf3_prefix(locs: list[int]) -> tuple[str, str]:
    """EN stem keeps parentheses: 4-(trifluoromethyl)."""
    if not locs:
        return "", ""
    n = len(locs)
    ls = ",".join(str(l) for l in sorted(locs))
    mult = _MULT_EN.get(n, "")
    return (
        f"{ls}-{mult}(trifluoromethyl)",
        f"{ls}-{_MULT_ZH.get(n, '')}三氟甲基",
    )


def _alkoxy_prefix(items: list[tuple[int, int]]) -> tuple[str, str]:
    """items: (locant, n). Group by n; alpha: ethoxy before methoxy."""
    if not items:
        return "", ""
    parts_e, parts_z = [], []
    for n in sorted({k for _, k in items}, key=lambda k: _ALKOXY_EN[k]):
        locs = [l for l, k in items if k == n]
        pe, pz = _same_stem_prefix(locs, _ALKOXY_EN[n], _ALKOXY_ZH[n])
        parts_e.append(pe)
        parts_z.append(pz)
    return "-".join(parts_e), "-".join(parts_z)


def _nitro_prefix(locs: list[int]) -> tuple[str, str]:
    return _same_stem_prefix(locs, "nitro", "硝基")


def _leaf_atoms_nitro(nb) -> set[int]:
    return {nb.GetIdx()} | {o.GetIdx() for o in _heavies(nb) if o.GetAtomicNum() == 8}


def _leaf_atoms_cf3(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {c_idx} | {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 9}


def _leaf_atoms_alkoxy(mol: Mol, o_idx: int, ring_i: int) -> set[int]:
    o = mol.GetAtomWithIdx(o_idx)
    nbs = _heavies(o)
    outer = next(n.GetIdx() for n in nbs if n.GetIdx() != ring_i)
    path = _linear_alkyl_atoms(mol, outer, o_idx) or []
    return {o_idx, *path}


def _leaf_atoms_one(mol: Mol, nb, ring_i: int) -> set[int]:
    kind = _depth2_kind(mol, nb, ring_i)
    if kind == "nitro":
        return _leaf_atoms_nitro(nb)
    if kind == "cf3":
        return _leaf_atoms_cf3(mol, nb.GetIdx())
    if kind in ("methoxy", "ethoxy"):
        return _leaf_atoms_alkoxy(mol, nb.GetIdx(), ring_i)
    return set()


def _depth2_atoms_on(mol: Mol, ph: set[int], nb_outside) -> set[int]:
    out: set[int] = set()
    for i in ph:
        for nb in nb_outside(mol, i, ph):
            out |= _leaf_atoms_one(mol, nb, i)
    return out
