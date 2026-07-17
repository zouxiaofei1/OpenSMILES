"""Open-chain simple dialkyl sulfoxide parent (P-63.3)."""
from __future__ import annotations


_SULFOXIDE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
    "has_amine", "has_alcohol", "has_thiol", "has_ether", "has_sulfide",
)


def _simple_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _is_open_sat, _no_fgs

    if info.get("has_ring"):
        return False
    return _is_open_sat(info) and _no_fgs(info, _SULFOXIDE_BAD)


def _arm_pair(info: dict, e: dict) -> tuple[list[int], list[int]] | None:
    from namepredict.layer2.chain_walk import _longest_from
    from namepredict.layer2.parent_core import _arm_ok, _hetero_open_chain

    mol, s_idx = info["mol"], e["s_idx"]
    if not _hetero_open_chain(mol, s_idx):
        return None
    arms = [_longest_from(mol, c, set()) for c in (e["c1"], e["c2"])]
    if not all(_arm_ok(mol, a, s_idx) for a in arms):
        return None
    return arms[0], arms[1]


def _make_parent(a1: list[int], a2: list[int], s_idx: int) -> dict:
    from namepredict.layer2.parent_core import _parent_dict

    parent = a1 if len(a1) >= len(a2) else a2
    return _parent_dict(parent, "sulfoxide", s_idx=s_idx, alkyl_ns=(len(a1), len(a2)))


def _sulfoxide_parent(info: dict) -> dict | None:
    sxs = info.get("sulfoxides") or []
    if len(sxs) != 1 or not _simple_ok(info):
        return None
    got = _arm_pair(info, sxs[0])
    return None if got is None else _make_parent(*got, sxs[0]["s_idx"])
