"""Secondary amino substituents on non-amine parents (e.g. N-benzyl on ethanol)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _make_amino(attach: int, n_idx: int, en: str = "amino", zh: str = "氨基",
                atoms: list[int] | None = None, paren: bool = False) -> dict:
    return {
        "kind": "amino", "attach_idx": attach,
        "atoms": [n_idx] if atoms is None else list(atoms),
        "en": en, "zh": zh, "paren": paren,
    }


def _n_aryl_named(mol: Mol, start: int, n_idx: int, is_ch2: bool):
    from namepredict.layer2.aryl_sub import (
        _benzyl_name, _ch2_ph_at, _halo_atoms_on, _phenyl_at, _phenyl_name,
    )
    ph = _ch2_ph_at(mol, start, n_idx) if is_ch2 else _phenyl_at(mol, start, n_idx)
    if ph is None:
        return None
    en0, zh0, _ = (_benzyl_name if is_ch2 else _phenyl_name)(mol, ph, start)
    base = [n_idx, start, *ph] if is_ch2 else [n_idx, *ph]
    return f"({en0})amino", f"({zh0})氨基", base + list(_halo_atoms_on(mol, ph))


def _sec_n_side_name(mol: Mol, start: int, n_idx: int):
    return _n_aryl_named(mol, start, n_idx, True) or _n_aryl_named(mol, start, n_idx, False)


def _sec_on_off(a: dict, chain_set: set[int]) -> tuple[int, int] | None:
    cs = a.get("c_idxs") or []
    if a.get("degree") != 2 or len(cs) != 2:
        return None
    on = [c for c in cs if c in chain_set]
    off = [c for c in cs if c not in chain_set]
    return (on[0], off[0]) if len(on) == 1 and len(off) == 1 else None


def _sec_amino_off_chain(mol: Mol, a: dict, chain_set: set[int]) -> dict | None:
    pair = _sec_on_off(a, chain_set)
    if pair is None:
        return None
    named = _sec_n_side_name(mol, pair[1], a["n_idx"])
    if named is None:
        return None
    en, zh, atoms = named
    return _make_amino(pair[0], a["n_idx"], en, zh, atoms, True)


def _one_amino(info: dict, a: dict, chain_set: set[int]) -> dict | None:
    if "c_idx" in a and a["c_idx"] in chain_set:
        return _make_amino(a["c_idx"], a["n_idx"])
    return _sec_amino_off_chain(info["mol"], a, chain_set)


def _extract_aminos(info: dict, parent: dict, parent_nh2_kinds: set) -> list[dict]:
    if parent.get("kind") in parent_nh2_kinds:
        return []
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for a in info.get("amines") or []:
        one = _one_amino(info, a, chain_set)
        if one is not None:
            out.append(one)
    return out
