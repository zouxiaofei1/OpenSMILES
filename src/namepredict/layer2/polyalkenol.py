"""Open-chain polyunsaturated monoalcohols (polyalkenols, P-63.1.1 / P-31.1)."""
from __future__ import annotations

from namepredict.layer2.aliph_fg import _aliphatic_entries


def _polyalkenol_ok(info: dict, bad: tuple) -> bool:
    """≥2 C=C, one aliphatic OH, no ring/alkyne/higher FG."""
    if info.get("has_ring") or info.get("triple_bonds"):
        return False
    if len(_aliphatic_entries(info, "hydroxyls")) != 1:
        return False
    if len(info.get("double_bonds") or []) < 2:
        return False
    return not any(info.get(k) for k in bad)


def _polyalkenol_atoms(info: dict, oh_c: int) -> list[int]:
    atoms = {oh_c}
    for db in info.get("double_bonds") or []:
        atoms.add(db["c1"])
        atoms.add(db["c2"])
    return list(atoms)


def _polyalkenol_parent(
    info: dict, bad: tuple, cover_fn, parent_fn, pairs_fn,
) -> dict | None:
    if not _polyalkenol_ok(info, bad):
        return None
    oh_c = _aliphatic_entries(info, "hydroxyls")[0]["c_idx"]
    chain = cover_fn(info["mol"], _polyalkenol_atoms(info, oh_c))
    if not chain or oh_c not in chain:
        return None
    return parent_fn(
        chain, "alkenol", oh_c_idx=oh_c, double_bonds=pairs_fn(info), mol=info["mol"],
    )
