"""L2 dialkyl sulfone functional parent (P-65.3.1.2)."""
from __future__ import annotations


_SULFONE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride",
    "has_sulfonate", "has_sulfonyl_chloride", "has_sulfonic_acid",
    "has_sulfonamide",
)


def _arm_n(mol, c_idx: int, forbid: int) -> int:
    from namepredict.layer2.chain_walk import _longest_from
    return len(_longest_from(mol, c_idx, {forbid}) or [c_idx])


def _arm_ok(mol, c_idx: int, forbid: int) -> bool:
    from namepredict.layer2.chain_walk import _longest_from
    from namepredict.layer2.parent_core import _arm_ok as _ok
    arm = _longest_from(mol, c_idx, {forbid}) or [c_idx]
    return _ok(mol, arm, forbid)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs
    return len(info.get("sulfones") or []) == 1 and _no_fgs(info, _SULFONE_BAD)


def _s_chain(mol, c_idx: int, s_idx: int) -> list[int]:
    from namepredict.layer2.chain_walk import _longest_from
    return _longest_from(mol, c_idx, {s_idx}) or [c_idx]


def _validate_sulfone_arms(mol, c1: int, c2: int, s_idx: int):
    if not _arm_ok(mol, c1, s_idx) or not _arm_ok(mol, c2, s_idx):
        return None
    a1 = _s_chain(mol, c1, s_idx)
    a2 = _s_chain(mol, c2, s_idx)
    n1, n2 = len(a1), len(a2)
    if n1 < 1 or n2 < 1 or n1 > 4 or n2 > 4:
        return None
    return n1, n2, a1, a2


def _sulfone_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e = info["sulfones"][0]
    v = _validate_sulfone_arms(info["mol"], e["c1"], e["c2"], e["s_idx"])
    if v is None:
        return None
    n1, n2, a1, a2 = v
    from namepredict.layer2.parent_core import _parent_dict
    return _parent_dict(a1 if n1 >= n2 else a2, "sulfone", s_idx=e["s_idx"], alkyl_ns=tuple(sorted((n1, n2), reverse=True)))
