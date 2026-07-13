"""Open-chain monounsaturated dicarboxylic acids (alkenedioic)."""
from __future__ import annotations

_DIACID_BAD = (
    "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_anhydride",
)


def _is_open_unsat_diacid(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _is_mono_alkene, _no_fgs

    if info.get("has_ring") or info.get("has_alkyne"):
        return False
    return _is_mono_alkene(info) and _no_fgs(info, _DIACID_BAD)


def _is_simple_alkenedioic(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _c_idxs

    if not _is_open_unsat_diacid(info):
        return False
    return _c_idxs(info.get("carboxyls") or [], 2) is not None


def _alkenedioic_atoms(info: dict) -> list[int] | None:
    from namepredict.layer2.parent_selector import _c_idxs

    cs = _c_idxs(info.get("carboxyls") or [], 2)
    if cs is None:
        return None
    db = (info.get("double_bonds") or [None])[0]
    if not db:
        return None
    return cs + [db["c1"], db["c2"]]


def _alkenedioic_meta(info: dict) -> dict:
    from namepredict.layer2.parent_selector import _c_idxs

    db = info["double_bonds"][0]
    cs = _c_idxs(info.get("carboxyls") or [], 2) or []
    return dict(cooh_c_idxs=cs, double_bond=(db["c1"], db["c2"]), mol=info["mol"])


def _alkenedioic_parent(info: dict) -> dict:
    from namepredict.layer2.parent_selector import _best_cover_pair, _parent_dict

    atoms = _alkenedioic_atoms(info) or []
    chain = _best_cover_pair(info["mol"], atoms)
    return _parent_dict(chain, "alkenedioic", **_alkenedioic_meta(info))
