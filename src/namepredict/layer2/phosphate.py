"""Simple monoalkyl phosphate and alkylphosphonic acid parents (P-67)."""
from __future__ import annotations


_P_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_anhydride",
    "has_thiol", "has_ether", "has_sulfide", "has_nitro",
)


def _simple_p_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs

    if info.get("has_ring") or info.get("has_alkene") or info.get("has_alkyne"):
        return False
    return _no_fgs(info, _P_BAD)


def _alkyl_chain(info: dict, c_idx: int) -> list[int]:
    from namepredict.layer2.parent_core import _longest_from

    return _longest_from(info["mol"], c_idx)


def _is_simple_phosphate(info: dict) -> bool:
    ps = info.get("phosphates") or []
    return len(ps) == 1 and not (info.get("phosphonics") or []) and _simple_p_ok(info)


def _is_simple_phosphonic(info: dict) -> bool:
    ps = info.get("phosphonics") or []
    return len(ps) == 1 and not (info.get("phosphates") or []) and _simple_p_ok(info)


def _chain_or_none(info: dict, c_idx: int) -> list[int] | None:
    chain = _alkyl_chain(info, c_idx)
    return chain if chain and c_idx in chain else None


def _phos_dict(chain, kind, **kw):
    from namepredict.layer2.parent_core import _parent_dict

    return _parent_dict(chain, kind, **kw)


def _phosphate_parent(info: dict) -> dict | None:
    if not _is_simple_phosphate(info):
        return None
    e = info["phosphates"][0]
    chain = _chain_or_none(info, e["alkoxy_c_idx"])
    if chain is None:
        return None
    return _phos_dict(
        chain, "phosphate", p_idx=e["p_idx"], alkoxy_c_idx=e["alkoxy_c_idx"],
        o_idx=e["o_idx"], alkoxy_n=len(chain),
    )


def _phosphonic_parent(info: dict) -> dict | None:
    if not _is_simple_phosphonic(info):
        return None
    e = info["phosphonics"][0]
    chain = _chain_or_none(info, e["c_idx"])
    if chain is None:
        return None
    return _phos_dict(chain, "phosphonic", p_idx=e["p_idx"], c_idx=e["c_idx"])
