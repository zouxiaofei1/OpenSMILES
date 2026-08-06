"""L1 detection of guanidine H2N–C(=NH)–NH2 / tautomers (P-66.4.1.2.1)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, H, N, S

def _bond(a, b):
    return a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())

def _bt(a, b):
    b0 = _bond(a, b)
    return b0.GetBondType() if b0 is not None else None

def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]

def _n_rest(n, carbon):
    """Heavy neighbors of N excluding the guanidine carbon."""
    return [x for x in _heavies(n) if x.GetIdx() != carbon.GetIdx()]

def _rest_atom_ok(x) -> bool:
    """Allow C or S (for Ar–SO2–N–guanidine)."""
    return x.GetAtomicNum() in (C, S)

def _n_ok(n, carbon) -> bool:
    if n.GetAtomicNum() != N or n.GetIsAromatic() or n.IsInRing():
        return False
    rest = _n_rest(n, carbon)
    return len(rest) <= 2 and all(_rest_atom_ok(x) for x in rest)

def _n_triple(carbon) -> list | None:
    ns = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == N]
    return ns if len(ns) == 3 else None

def _bond_pattern_ok(carbon, ns) -> bool:
    """Standard tautomer: one C=N + two C–N, or three single (charged)."""
    bts = [_bt(carbon, n) for n in ns]
    if any(b is None for b in bts):
        return False
    n_dbl = sum(1 for b in bts if b == BondType.DOUBLE)
    n_sgl = sum(1 for b in bts if b == BondType.SINGLE)
    return (n_dbl == 1 and n_sgl == 2) or (n_dbl == 0 and n_sgl == 3)

def _not_in_ring(carbon, ns) -> bool:
    if carbon.IsInRing() or carbon.GetIsAromatic():
        return False
    return not any(n.IsInRing() for n in ns)

def _pack_n(n, carbon) -> dict:
    rest = _n_rest(n, carbon)
    return {
        "n_idx": n.GetIdx(),
        "rest_idxs": [x.GetIdx() for x in rest],
        "rest_zs": [x.GetAtomicNum() for x in rest],
    }

def _pack(carbon, ns: list) -> dict:
    sides = [_pack_n(n, carbon) for n in ns]
    return {
        "c_idx": carbon.GetIdx(),
        "n_idxs": [s["n_idx"] for s in sides],
        "n_sides": sides,
    }

def _entry(atom) -> dict | None:
    if atom.GetAtomicNum() != C:
        return None
    ns = _n_triple(atom)
    if ns is None or not _bond_pattern_ok(atom, ns):
        return None
    if not all(_n_ok(n, atom) for n in ns) or not _not_in_ring(atom, ns):
        return None
    return _pack(atom, ns)

def guanidine_entries(mol: Mol) -> list[dict]:
    return [e for a in mol.GetAtoms() if (e := _entry(a)) is not None]

def is_guanidine_n(atom) -> bool:
    """True if atom is one of the three N on a guanidine carbon."""
    if atom.GetAtomicNum() != N:
        return False
    return any(
        n.GetAtomicNum() == C and _entry(n) is not None
        for n in atom.GetNeighbors()
    )
