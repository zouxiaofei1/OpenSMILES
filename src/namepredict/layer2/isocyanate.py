"""L2 open-chain mono isocyanate / isothiocyanate parents (P-61.9)."""
from __future__ import annotations

from namepredict.layer2.chain_walk import _chain_through
from namepredict.layer2.parent_selector import (
    _is_mono_fg, _is_open_sat, _no_fgs, _parent_dict,
)

_ISO_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
    "has_amine", "has_alcohol", "has_thiol",
)


def _r_open(mol, r_c: int) -> bool:
    return not mol.GetAtomWithIdx(r_c).IsInRing()


def _r_not_arom(mol, r_c: int) -> bool:
    return not mol.GetAtomWithIdx(r_c).GetIsAromatic()


def _mono_iso_ok(info: dict, flag: str, key: str) -> bool:
    if not _is_open_sat(info) or not _is_mono_fg(info, flag, key):
        return False
    if not _no_fgs(info, _ISO_BAD):
        return False
    e, mol = info[key][0], info["mol"]
    r = e["r_c_idx"]
    return _r_open(mol, r) and _r_not_arom(mol, r)


def _iso_parent(info: dict, key: str, kind: str) -> dict:
    e = info[key][0]
    r = e["r_c_idx"]
    return _parent_dict(
        _chain_through(info, r), kind,
        r_c_idx=r, n_idx=e["n_idx"], c_idx=e["c_idx"], x_idx=e["x_idx"],
    )


def _try_isocyanate(info: dict) -> dict | None:
    if not _mono_iso_ok(info, "has_isocyanate", "isocyanates"):
        return None
    return _iso_parent(info, "isocyanates", "isocyanate")


def _try_isothiocyanate(info: dict) -> dict | None:
    if not _mono_iso_ok(info, "has_isothiocyanate", "isothiocyanates"):
        return None
    return _iso_parent(info, "isothiocyanates", "isothiocyanate")
