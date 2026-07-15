"""L1 detection of urea R2N–C(=O)–NR2 / enol tautomer (P-66.1.6.1.1)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol


def _bond(a, b):
    return a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())


def _bt(a, b):
    b0 = _bond(a, b)
    return b0.GetBondType() if b0 is not None else None


def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _n_rest(n, carbon):
    return [x for x in _heavies(n) if x.GetIdx() != carbon.GetIdx()]


def _n_ok(n, carbon) -> bool:
    if n.GetAtomicNum() != 7:
        return False
    rest = _n_rest(n, carbon)
    return len(rest) <= 2 and all(x.GetAtomicNum() == 6 for x in rest)


def _dbl_o(carbon):
    return next(
        (n for n in carbon.GetNeighbors()
         if n.GetAtomicNum() == 8 and _bt(carbon, n) == BondType.DOUBLE),
        None,
    )


def _sgl_oh(carbon):
    return next(
        (n for n in carbon.GetNeighbors()
         if n.GetAtomicNum() == 8 and n.GetTotalNumHs() >= 1
         and _bt(carbon, n) == BondType.SINGLE),
        None,
    )


def _n_pair(carbon, want_dbl: bool) -> list | None:
    """Two N neighbors: keto both single; enol one double + one single."""
    ns = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7]
    if len(ns) != 2:
        return None
    bts = {_bt(carbon, n) for n in ns}
    if want_dbl:
        ok = bts == {BondType.DOUBLE, BondType.SINGLE}
    else:
        ok = bts == {BondType.SINGLE}
    return ns if ok and all(_n_ok(n, carbon) for n in ns) else None


def _not_in_ring_urea(carbon, ns) -> bool:
    """Skip cyclic ureas (imidazolidinone etc.): C or both N in same ring."""
    if carbon.IsInRing():
        return False
    return not (ns[0].IsInRing() and ns[1].IsInRing())


def _pack(carbon, o, ns: list) -> dict:
    n0, n1 = ns[0], ns[1]
    return {
        "c_idx": carbon.GetIdx(),
        "o_idx": o.GetIdx(),
        "n1_idx": n0.GetIdx(),
        "n2_idx": n1.GetIdx(),
        "n1_c_idxs": [x.GetIdx() for x in _n_rest(n0, carbon)],
        "n2_c_idxs": [x.GetIdx() for x in _n_rest(n1, carbon)],
        "enol": o.GetTotalNumHs() >= 1,
    }


def _keto_entry(carbon) -> dict | None:
    o = _dbl_o(carbon)
    ns = _n_pair(carbon, want_dbl=False)
    if o is None or ns is None or not _not_in_ring_urea(carbon, ns):
        return None
    return _pack(carbon, o, ns)


def _enol_entry(carbon) -> dict | None:
    o = _sgl_oh(carbon)
    ns = _n_pair(carbon, want_dbl=True)
    if o is None or ns is None or not _not_in_ring_urea(carbon, ns):
        return None
    return _pack(carbon, o, ns)


def _is_urea_carbon(atom) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    return _keto_entry(atom) is not None or _enol_entry(atom) is not None


def _entry(atom) -> dict | None:
    return _keto_entry(atom) or _enol_entry(atom)


def urea_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry(a)) is not None]


def is_urea_oh(atom) -> bool:
    """True if atom is enol-OH oxygen of a urea carbon."""
    if atom.GetAtomicNum() != 8 or atom.GetTotalNumHs() < 1:
        return False
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() == 6 and _enol_entry(n) is not None:
            return True
    return False
