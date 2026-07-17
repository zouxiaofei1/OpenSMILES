"""Open-chain mono-/polyunsaturated dicarboxylic acids (alkenedioic)."""
from __future__ import annotations


_DIACID_BAD = (
    "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_anhydride",
)


def _is_open_unsat_diacid(info: dict) -> bool:
    """Open-chain diacid with ≥1 C=C and no competing higher FG."""
    from namepredict.layer2.parent_core import _no_fgs

    if info.get("has_ring") or info.get("has_alkyne"):
        return False
    dbs = info.get("double_bonds") or []
    return len(dbs) >= 1 and _no_fgs(info, _DIACID_BAD)


def _is_simple_alkenedioic(info: dict) -> bool:
    from namepredict.layer2.parent_core import _c_idxs

    if not _is_open_unsat_diacid(info):
        return False
    return _c_idxs(info.get("carboxyls") or [], 2) is not None


def _alkenedioic_atoms(info: dict) -> list[int] | None:
    from namepredict.layer2.parent_core import _c_idxs

    cs = _c_idxs(info.get("carboxyls") or [], 2)
    if cs is None:
        return None
    atoms = list(cs)
    for db in info.get("double_bonds") or []:
        atoms.extend([db["c1"], db["c2"]])
    return atoms


def _db_meta(info: dict) -> dict:
    """Single C=C → double_bond; multi → double_bonds (polyene style)."""
    from namepredict.layer2.parent_core import _db_pairs

    pairs = _db_pairs(info)
    if len(pairs) >= 2:
        return {"double_bonds": pairs}
    if len(pairs) == 1:
        return {"double_bond": pairs[0]}
    return {}


def _alkenedioic_meta(info: dict) -> dict:
    from namepredict.layer2.parent_core import _c_idxs

    cs = _c_idxs(info.get("carboxyls") or [], 2) or []
    return dict(cooh_c_idxs=cs, mol=info["mol"], **_db_meta(info))


def _alkenedioic_parent(info: dict) -> dict:
    from namepredict.layer2.parent_core import _best_cover_pair, _parent_dict

    atoms = _alkenedioic_atoms(info) or []
    chain = _best_cover_pair(info["mol"], atoms)
    return _parent_dict(chain, "diacid", **_alkenedioic_meta(info))
