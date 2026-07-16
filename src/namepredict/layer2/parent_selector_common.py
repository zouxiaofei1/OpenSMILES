"""Shared functional-group eligibility predicates for parent selection."""
from __future__ import annotations

def _no_fgs(info: dict, keys: tuple) -> bool:
    return not any(info.get(k) for k in keys)
_CORE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
)
_DIOL_BAD = _CORE_BAD + ("has_amine",)
_DIACID_BAD = (  # P-65.1.2: hydroxy/amino/oxo are prefixes, not competing
    "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_anhydride",
)
_DIAMINE_BAD = _CORE_BAD + ("has_alcohol",)
_DIONE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_amine", "has_alcohol", "has_acyl_chloride", "has_anhydride",
)
_ANHYDRIDE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_acyl_chloride",
)
