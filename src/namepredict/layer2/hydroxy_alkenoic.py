"""Open-chain hydroxyalkenoic acids (P-65.1.2 / P-31.1)."""
from __future__ import annotations

from namepredict.layer2.aliph_fg import _aliphatic_entries


def _hy_alkenoic_ok(info: dict, bad: tuple) -> bool:
    """Mono COOH + mono aliphatic OH + ≥1 C=C; no ring/alkyne/higher FG."""
    if info.get("has_ring") or info.get("triple_bonds") or info.get("has_alkyne"):
        return False
    if len(info.get("carboxyls") or []) != 1:
        return False
    if len(_aliphatic_entries(info, "hydroxyls")) != 1:
        return False
    if not (info.get("double_bonds") or []):
        return False
    return not any(info.get(k) for k in bad)


def _hy_alkenoic_atoms(info: dict, cooh: int, oh_c: int) -> list[int]:
    atoms = {cooh, oh_c}
    for db in info.get("double_bonds") or []:
        atoms.add(db["c1"])
        atoms.add(db["c2"])
    return list(atoms)


def _hy_db_kw(dbs: list) -> dict:
    if len(dbs) >= 2:
        return {"double_bonds": dbs}
    return {"double_bond": dbs[0]}


def _hy_alkenoic_parent(
    info: dict, bad: tuple, cover_fn, parent_fn, pairs_fn,
) -> dict | None:
    if not _hy_alkenoic_ok(info, bad):
        return None
    cooh = info["carboxyls"][0]["c_idx"]
    oh_c = _aliphatic_entries(info, "hydroxyls")[0]["c_idx"]
    chain = cover_fn(info["mol"], _hy_alkenoic_atoms(info, cooh, oh_c))
    if not chain or cooh not in chain or oh_c not in chain:
        return None
    return parent_fn(
        chain, "alkenoic_acid", cooh_c_idx=cooh, mol=info["mol"], **_hy_db_kw(pairs_fn(info)),
    )
