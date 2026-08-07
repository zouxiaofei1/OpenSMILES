"""Amide/amine N-side substituent dicts with atom provenance for coverage."""
from __future__ import annotations

from rdkit.Chem import Mol

_N_ALKYL_EN = {1: "N-methyl", 2: "N-ethyl", 3: "N-propyl", 4: "N-butyl"}
_N_ALKYL_ZH = {1: "N-甲基", 2: "N-乙基", 3: "N-丙基", 4: "N-丁基"}
_N_STEM_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_STEM_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
_AMIDE_KINDS = frozenset({"amide", "benzamide"})


def _n_alkyl_sub(en: str, zh: str, attach: int, n: int, atoms: list[int] | None = None) -> dict:
    return {
        "kind": "n_alkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": list(atoms or []), "en": en, "zh": zh,
    }


def _n_phenyl_sub(attach: int, atoms: list[int] | None = None) -> dict:
    return {
        "kind": "n_phenyl", "n_carbons": 6, "attach_idx": attach,
        "atoms": list(atoms or []), "en": "N-phenyl", "zh": "N-苯基",
    }


def _tert_n_prefix(ns: list[int]) -> tuple[str, str] | None:
    if len(ns) != 2 or any(n not in _N_STEM_EN for n in ns):
        return None
    a, b = ns
    if a == b:
        return f"N,N-di{_N_STEM_EN[a]}", f"N,N-二{_N_STEM_ZH[a]}"
    x, y = sorted(ns, key=lambda n: _N_STEM_EN[n])
    return f"N-{_N_STEM_EN[x]}-N-{_N_STEM_EN[y]}", f"N-{_N_STEM_ZH[x]}-N-{_N_STEM_ZH[y]}"


def _n_alkyl_prefix(parent: dict) -> tuple[str, str, int] | None:
    kind, n = parent.get("kind"), parent.get("n_alkyl_n")
    if kind in ("sec_amine", "amide", "benzamide") and n in _N_ALKYL_EN:
        return _N_ALKYL_EN[n], _N_ALKYL_ZH[n], n
    if kind in ("tert_amine", "amide", "benzamide"):
        pref = _tert_n_prefix(list(parent.get("n_alkyl_ns") or []))
        return (*pref, 0) if pref else None
    return None


def _amide_n_idx(parent: dict, mol: Mol | None) -> int | None:
    c = parent.get("amide_c_idx")
    if c is None or mol is None:
        return None
    for n in mol.GetAtomWithIdx(int(c)).GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _n_side_roots(mol: Mol, n_idx: int, owned: frozenset[int]) -> list[int]:
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(n_idx).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in owned
    ]


def _cut_n_sides(mol: Mol, n_idx: int, owned: frozenset[int]) -> list[int]:
    from namepredict.tools.block_cut import cut_block

    atoms: set[int] = set()
    for root in _n_side_roots(mol, n_idx, owned):
        block = cut_block(mol, root, owned)
        if block:
            atoms |= set(block)
    return sorted(atoms)


def _n_side_atoms(mol: Mol | None, parent: dict) -> list[int]:
    owned = parent.get("owned_atoms")
    n_idx = _amide_n_idx(parent, mol) if mol is not None else None
    if mol is None or not isinstance(owned, frozenset) or n_idx is None or n_idx not in owned:
        return []
    return _cut_n_sides(mol, n_idx, owned)


def _n_attach(parent: dict) -> int | None:
    kind = parent.get("kind")
    key = "amide_c_idx" if kind in _AMIDE_KINDS else "amine_c_idx"
    attach = parent.get(key)
    return parent.get("ring_attach_idx", attach) if kind == "benzamide" else attach


def extract_n_alkyl(parent: dict, mol: Mol | None = None) -> list[dict]:
    attach, pref = _n_attach(parent), _n_alkyl_prefix(parent)
    if attach is None or not pref:
        return []
    kind = parent.get("kind")
    atoms = _n_side_atoms(mol, parent) if kind in _AMIDE_KINDS else []
    return [_n_alkyl_sub(*pref[:2], attach, pref[2], atoms)]


def extract_n_phenyl(parent: dict, mol: Mol | None = None) -> list[dict]:
    if parent.get("kind") not in _AMIDE_KINDS or not parent.get("n_phenyl"):
        return []
    attach = _n_attach(parent)
    return [] if attach is None else [_n_phenyl_sub(attach, _n_side_atoms(mol, parent))]
