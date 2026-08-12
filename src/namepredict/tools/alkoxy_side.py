"""Shared ester/carbamate alkoxy-side topology (P-65.6)."""
from __future__ import annotations


def _heavy_nbs(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _is_methyl_c(atom, center_idx: int) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    heavy = _heavy_nbs(atom)
    return len(heavy) == 1 and heavy[0].GetIdx() == center_idx


def _methyl_count(c, o_idx: int, c_idx: int) -> int:
    return sum(
        1 for n in _heavy_nbs(c)
        if n.GetIdx() != o_idx and _is_methyl_c(n, c_idx)
    )


def _is_branched_alkyl_c(mol, c_idx: int, o_idx: int, *, total_nbs: int, n_methyl: int) -> bool:
    """True if `c_idx` is C with `total_nbs` heavy neighbours, `n_methyl` of them methyl."""
    c = mol.GetAtomWithIdx(c_idx)
    nbs = _heavy_nbs(c)
    return c.GetAtomicNum() == 6 and len(nbs) == total_nbs and _methyl_count(c, o_idx, c_idx) == n_methyl


def _is_tert_butyl_c(mol, c_idx: int, o_idx: int) -> bool:
    """C attached to O is C(CH3)3."""
    return _is_branched_alkyl_c(mol, c_idx, o_idx, total_nbs=4, n_methyl=3)


def _is_isopropyl_c(mol, c_idx: int, o_idx: int) -> bool:
    """C attached to O is CH(CH3)2."""
    return _is_branched_alkyl_c(mol, c_idx, o_idx, total_nbs=3, n_methyl=2)


def _sole_c6_at(mol, attach: int, parent: int):
    from namepredict.layer3.aryl_sub import _c6_rings_at, _is_unfused_benzene_ring
    hits = _c6_rings_at(mol, attach)
    if len(hits) != 1 or parent in hits[0]:
        return None
    ring = hits[0]
    return ring if _is_unfused_benzene_ring(mol, ring) else None


def _ring_only_parent_link(mol, ring, attach: int, parent: int) -> bool:
    from namepredict.layer3.aryl_sub import _nb_outside
    for i in ring:
        outs = [n.GetIdx() for n in _nb_outside(mol, i, ring)]
        if i == attach:
            if outs != [parent]:
                return False
        elif outs:
            return False
    return True


def _bare_ph_ring(mol, attach: int, parent: int):
    """Unfused unsubstituted Ph: sole outside bond is attach→parent."""
    ring = _sole_c6_at(mol, attach, parent)
    return ring if ring is not None and _ring_only_parent_link(mol, ring, attach, parent) else None


def _is_phenyl_c(mol, c_idx: int, o_idx: int) -> bool:
    """O–Ph bare phenyl only (not heteroaryl / substituted / fused)."""
    return _bare_ph_ring(mol, c_idx, o_idx) is not None


def _benzyl_other(nbs, o_idx: int):
    other = [n for n in nbs if n.GetIdx() != o_idx]
    return other[0] if len(other) == 1 else None


def _is_open_ch2(atom) -> bool:
    if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing():
        return False
    return len(_heavy_nbs(atom)) == 2


def _benzyl_ph_attach(mol, ch2_idx: int, o_idx: int) -> int | None:
    """If ch2 is O–CH2–Ph (bare Ph), return Ph attach carbon."""
    c = mol.GetAtomWithIdx(ch2_idx)
    if not _is_open_ch2(c):
        return None
    a = _benzyl_other(_heavy_nbs(c), o_idx)
    if a is None or a.GetAtomicNum() != 6:
        return None
    return a.GetIdx() if _bare_ph_ring(mol, a.GetIdx(), ch2_idx) else None


def _linear_n(mol, start: int) -> int:
    from namepredict.tools.chain import _longest_from
    chain = _longest_from(mol, start)
    return len(chain) if chain else 1


def _special_alkoxy(mol, o_idx: int, alkoxy_c: int) -> dict | None:
    if _is_tert_butyl_c(mol, alkoxy_c, o_idx):
        return {"alkoxy_en": "tert-butyl", "alkoxy_zh": "叔丁", "alkoxy_n": None}
    if _is_isopropyl_c(mol, alkoxy_c, o_idx):
        return {"alkoxy_en": "propan-2-yl", "alkoxy_zh": "丙-2-基", "alkoxy_n": None}
    if _is_phenyl_c(mol, alkoxy_c, o_idx):
        return {"alkoxy_en": "phenyl", "alkoxy_zh": "苯", "alkoxy_n": None}
    return None


def _benzyl_or_linear(mol, o_idx: int, alkoxy_c: int) -> dict:
    if _benzyl_ph_attach(mol, alkoxy_c, o_idx) is not None:
        return {"alkoxy_en": "benzyl", "alkoxy_zh": "苄", "alkoxy_n": None}
    return {"alkoxy_en": "", "alkoxy_zh": "", "alkoxy_n": _linear_n(mol, alkoxy_c)}


def classify_alkoxy(mol, o_idx: int, alkoxy_c: int) -> dict:
    """Return labels for ester/carbamate O-side.

    Keys: alkoxy_en, alkoxy_zh, alkoxy_n (int|None for non-linear specials).
    """
    return _special_alkoxy(mol, o_idx, alkoxy_c) or _benzyl_or_linear(mol, o_idx, alkoxy_c)
