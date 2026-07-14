"""Aliphatic-only FG attachment filters for chain poly-FG parents (L2)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _is_arom_c(mol: Mol, c_idx: int) -> bool:
    a = mol.GetAtomWithIdx(c_idx)
    return a.GetAtomicNum() == 6 and a.GetIsAromatic()


def _aliphatic_entries(info: dict, ekey: str) -> list:
    """FG entries whose attachment carbon is non-aromatic."""
    mol: Mol = info["mol"]
    return [
        e for e in (info.get(ekey) or [])
        if "c_idx" in e and not _is_arom_c(mol, e["c_idx"])
    ]


def _c_idxs(entries, n: int) -> list[int] | None:
    xs = [e["c_idx"] for e in entries or [] if "c_idx" in e]
    return xs if len(xs) == n and len(set(xs)) == n else None


def _aliph_c_idxs(info: dict, ekey: str, n: int) -> list[int] | None:
    return _c_idxs(_aliphatic_entries(info, ekey), n)
