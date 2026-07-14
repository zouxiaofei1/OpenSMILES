"""Open-chain prefix-alkenoic acids (hydroxy/amino/oxo; P-65.1.2 / P-31.1)."""
from __future__ import annotations

from namepredict.layer2.aliph_fg import _aliphatic_entries


def _primary_amines(info: dict) -> list:
    """Aliphatic primary amines (degree 1 with c_idx)."""
    return [e for e in _aliphatic_entries(info, "amines") if e.get("degree") == 1]


def _prefix_counts(info: dict) -> tuple[int, int, int]:
    """(n_aliph_OH, n_primary_amine, n_ketone)."""
    return (
        len(_aliphatic_entries(info, "hydroxyls")),
        len(_primary_amines(info)),
        len(_aliphatic_entries(info, "ketones")),
    )


def _prefix_fg_ok(info: dict) -> bool:
    """≤1 each of OH/amine/ketone and at least one present."""
    oh, am, ox = _prefix_counts(info)
    return oh <= 1 and am <= 1 and ox <= 1 and (oh + am + ox) >= 1


def _open_parent_atoms(mol, cooh: int, dbs: list) -> bool:
    """Principal COOH carbon and all C=C ends must be acyclic (aryl OK)."""
    atoms = {int(cooh)}
    for db in dbs:
        atoms.add(int(db["c1"]))
        atoms.add(int(db["c2"]))
    return all(not mol.GetAtomWithIdx(a).IsInRing() for a in atoms)


def _prefix_ok(info: dict, bad: tuple) -> bool:
    """Mono COOH + ≥1 open-chain C=C + optional mono OH/amine/ketone prefixes."""
    if info.get("triple_bonds") or info.get("has_alkyne"):
        return False
    carbs = info.get("carboxyls") or []
    dbs = info.get("double_bonds") or []
    if len(carbs) != 1 or not dbs:
        return False
    if not _open_parent_atoms(info["mol"], carbs[0]["c_idx"], dbs):
        return False
    return _prefix_fg_ok(info) and not any(info.get(k) for k in bad)


def _add_fg_cs(atoms: set, info: dict) -> None:
    for key in ("hydroxyls", "ketones"):
        for e in _aliphatic_entries(info, key):
            atoms.add(e["c_idx"])
    for e in _primary_amines(info):
        atoms.add(e["c_idx"])


def _prefix_atoms(info: dict, cooh: int) -> list[int]:
    """COOH + all C=C ends + aliphatic OH/amine/ketone carbons."""
    atoms = {cooh}
    for db in info.get("double_bonds") or []:
        atoms.add(db["c1"])
        atoms.add(db["c2"])
    _add_fg_cs(atoms, info)
    return list(atoms)


def _db_kw(dbs: list) -> dict:
    if len(dbs) >= 2:
        return {"double_bonds": dbs}
    return {"double_bond": dbs[0]}


def _prefix_alkenoic_parent(
    info: dict, bad: tuple, cover_fn, parent_fn, pairs_fn,
) -> dict | None:
    if not _prefix_ok(info, bad):
        return None
    cooh = info["carboxyls"][0]["c_idx"]
    need = set(_prefix_atoms(info, cooh))
    chain = cover_fn(info["mol"], list(need))
    if not chain or not need.issubset(chain):
        return None
    return parent_fn(
        chain, "alkenoic_acid", cooh_c_idx=cooh, mol=info["mol"], **_db_kw(pairs_fn(info)),
    )
