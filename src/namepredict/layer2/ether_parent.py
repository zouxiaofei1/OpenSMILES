"""Open-chain dialkyl ether parents (P-63.2.2).

Linear C1–C4 arms: alkoxyalkane / retained sym ether (existing).
Branched specials: isopropyl / hexafluoroisopropyl → functional-class fields.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.parent_core import (
    _arm_ok,
    _hetero_open_chain,
    _is_open_sat,
    _longest_from,
    _no_fgs,
    _parent_dict,
)

_ETHER_BAD = (
    "has_acid", "has_ester", "has_amide", "has_aldehyde", "has_ketone",
    "has_nitrile", "has_acyl_chloride", "has_anhydride",
    "has_amine", "has_alcohol", "has_thiol", "has_sulfide",
)


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _is_cf3(mol: Mol, c_idx: int, parent: int) -> bool:
    """C is CF3 attached only to parent (+3F)."""
    nbs = _heavies(mol.GetAtomWithIdx(c_idx))
    fs = [n for n in nbs if n.GetAtomicNum() == 9]
    rest = [n.GetIdx() for n in nbs if n.GetAtomicNum() != 9]
    return len(fs) == 3 and rest == [parent]


def _is_hfip(mol: Mol, c_idx: int, o_idx: int) -> bool:
    """O–CH(CF3)2 (1,1,1,3,3,3-hexafluoropropan-2-yl)."""
    nbs = _heavies(mol.GetAtomWithIdx(c_idx))
    if len(nbs) != 3:
        return False
    ids = {n.GetIdx() for n in nbs}
    if o_idx not in ids:
        return False
    cs = [n.GetIdx() for n in nbs if n.GetAtomicNum() == 6]
    return len(cs) == 2 and all(_is_cf3(mol, c, c_idx) for c in cs)


def _hfip_chain(mol: Mol, c_idx: int, o_idx: int) -> list[int]:
    cs = [n.GetIdx() for n in _heavies(mol.GetAtomWithIdx(c_idx)) if n.GetAtomicNum() == 6]
    return [c_idx, *cs]


def _is_ipr(mol: Mol, c_idx: int, o_idx: int) -> bool:
    from namepredict.layer2.alkoxy_side import _is_isopropyl_c
    return _is_isopropyl_c(mol, c_idx, o_idx)


def _ipr_chain(mol: Mol, c_idx: int, o_idx: int) -> list[int]:
    cs = [n.GetIdx() for n in _heavies(mol.GetAtomWithIdx(c_idx)) if n.GetAtomicNum() == 6]
    return [c_idx, *cs]


def _linear_arm(mol: Mol, c_idx: int, o_idx: int) -> dict | None:
    arm = _longest_from(mol, c_idx, set())
    if not _arm_ok(mol, arm, o_idx):
        return None
    return {"tag": "n", "n": len(arm), "chain": arm}


def _special_arm(mol: Mol, c_idx: int, o_idx: int) -> dict | None:
    if _is_hfip(mol, c_idx, o_idx):
        return {"tag": "hfip", "n": 0, "chain": _hfip_chain(mol, c_idx, o_idx)}
    if _is_ipr(mol, c_idx, o_idx):
        return {"tag": "ipr", "n": 0, "chain": _ipr_chain(mol, c_idx, o_idx)}
    return None


def _classify_arm(mol: Mol, c_idx: int, o_idx: int) -> dict | None:
    return _special_arm(mol, c_idx, o_idx) or _linear_arm(mol, c_idx, o_idx)


def _linear_ether_parent(a1: dict, a2: dict, e: dict) -> dict:
    parent, short = (a1, a2) if len(a1["chain"]) >= len(a2["chain"]) else (a2, a1)
    return _parent_dict(
        parent["chain"], "ether",
        ether_c_idx=parent["chain"][0], o_idx=e["o_idx"], alkoxy_n=short["n"],
    )


def _func_ether_parent(a1: dict, a2: dict, e: dict) -> dict:
    """Functional-class ether: store both arm tags for L5."""
    chain = a1["chain"] if len(a1["chain"]) >= len(a2["chain"]) else a2["chain"]
    arms = ((a1["tag"], a1["n"]), (a2["tag"], a2["n"]))
    return _parent_dict(
        chain, "ether",
        ether_c_idx=chain[0], o_idx=e["o_idx"], ether_arms=arms, alkoxy_n=None,
    )


def _ether_entry(info: dict) -> tuple[dict, Mol] | None:
    ets = info.get("ethers") or []
    if len(ets) != 1 or not _is_open_sat(info) or not _no_fgs(info, _ETHER_BAD):
        return None
    e, mol = ets[0], info["mol"]
    return (e, mol) if _hetero_open_chain(mol, e["o_idx"]) else None


def _ether_from_arms(a1: dict, a2: dict, e: dict) -> dict:
    if a1["tag"] == "n" and a2["tag"] == "n":
        return _linear_ether_parent(a1, a2, e)
    return _func_ether_parent(a1, a2, e)


def _ether_parent(info: dict) -> dict | None:
    got = _ether_entry(info)
    if got is None:
        return None
    e, mol = got
    a1 = _classify_arm(mol, e["c1"], e["o_idx"])
    a2 = _classify_arm(mol, e["c2"], e["o_idx"])
    return None if a1 is None or a2 is None else _ether_from_arms(a1, a2, e)
