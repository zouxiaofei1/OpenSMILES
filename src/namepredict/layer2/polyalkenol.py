"""Open-chain polyunsaturated monoalcohols (polyalkenols, P-63.1.1 / P-31.1)."""
from __future__ import annotations

from namepredict.layer2.aliph_fg import _aliphatic_entries


def _open_poly_atoms(mol, oh_c: int, dbs: list) -> bool:
    """OH attach carbon and all C=C ends must be acyclic (aryl OK)."""
    atoms = {int(oh_c)}
    for db in dbs:
        atoms.add(int(db["c1"]))
        atoms.add(int(db["c2"]))
    return all(not mol.GetAtomWithIdx(a).IsInRing() for a in atoms)


def _polyalkenol_ok(info: dict, bad: tuple) -> bool:
    """≥2 open-chain C=C, one aliphatic OH; molecule rings OK if parent open."""
    if info.get("triple_bonds"):
        return False
    ohs = _aliphatic_entries(info, "hydroxyls")
    dbs = info.get("double_bonds") or []
    if len(ohs) != 1 or len(dbs) < 2:
        return False
    if not _open_poly_atoms(info["mol"], ohs[0]["c_idx"], dbs):
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
        chain, "alcohol", oh_c_idx=oh_c, double_bonds=pairs_fn(info), mol=info["mol"],
    )
