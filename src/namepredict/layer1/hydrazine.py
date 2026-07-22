"""L1 detection of open-chain hydrazine N–N (P-68.3.1.2)."""
from __future__ import annotations

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, H, N, O


def _bond(a, b):
    return a.GetOwningMol().GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())


def _is_sgl(a, b) -> bool:
    b0 = _bond(a, b)
    return b0 is not None and b0.GetBondType() == BondType.SINGLE


def _heavies(atom):
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]


def _rest_cs(n, other_n):
    """Carbon substituents on N excluding the hydrazine partner N."""
    return [
        x for x in _heavies(n)
        if x.GetIdx() != other_n.GetIdx() and x.GetAtomicNum() == C
    ]


def _rest_ok(n, other_n) -> bool:
    """N may only carry H and C (no O/N etc. beyond partner)."""
    rest = [x for x in _heavies(n) if x.GetIdx() != other_n.GetIdx()]
    return all(x.GetAtomicNum() == C for x in rest)


def _n_ok(n, other) -> bool:
    if n.GetAtomicNum() != N or n.GetIsAromatic() or n.IsInRing():
        return False
    return _rest_ok(n, other)


def _has_dbl_o(carbon) -> bool:
    for o in carbon.GetNeighbors():
        if o.GetAtomicNum() != O:
            continue
        b = _bond(carbon, o)
        if b is not None and b.GetBondType() == BondType.DOUBLE:
            return True
    return False


def _is_amide_like(n) -> bool:
    """True if N is attached to C(=O) (amide / urea / hydrazide)."""
    return any(
        c.GetAtomicNum() == C and _has_dbl_o(c) for c in n.GetNeighbors()
    )


def _pair_ok(n1, n2) -> bool:
    if not (_n_ok(n1, n2) and _n_ok(n2, n1) and _is_sgl(n1, n2)):
        return False
    return not (_is_amide_like(n1) or _is_amide_like(n2))


def _pack(n1, n2) -> dict:
    return {
        "n1_idx": n1.GetIdx(),
        "n2_idx": n2.GetIdx(),
        "n1_c_idxs": [c.GetIdx() for c in _rest_cs(n1, n2)],
        "n2_c_idxs": [c.GetIdx() for c in _rest_cs(n2, n1)],
    }


def _entry_for_bond(bond) -> dict | None:
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    if a.GetAtomicNum() != N or b.GetAtomicNum() != N:
        return None
    if not _pair_ok(a, b):
        return None
    # Canonical order by index so each N–N is listed once.
    n1, n2 = (a, b) if a.GetIdx() < b.GetIdx() else (b, a)
    return _pack(n1, n2)


def hydrazine_entries(mol: Mol) -> list[dict]:
    return [e for b in mol.GetBonds() if (e := _entry_for_bond(b)) is not None]
