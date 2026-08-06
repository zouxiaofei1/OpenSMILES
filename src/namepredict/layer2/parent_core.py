"""L2 parent assembly / chain-walk / gate helpers (single authority).

Producers import helpers from here. parent_selector keeps FG try + select_parent
and may thin-re-export for compatibility.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.fg_helpers import _aliphatic_entries, _c_idxs, _no_fgs
from namepredict.layer2.chain_walk import (
    _better,
    _carbon_neighbors,
    _chain_through,
    _chain_through_two,
    _longest_chain,
    _longest_from,
)
# Re-export chain / gate primitives used by producers.
__all__ = [
    "_parent_dict",
    "_no_fgs",
    "_longest_from",
    "_longest_chain",
    "_c_idxs",
    "_is_open_sat",
    "_hetero_open_chain",
    "_side_carbons",
    "_arm_ok",
    "_is_mono_fg",
    "_is_mono_alkene",
    "_is_mono_alkyne",
    "_covers",
    "_best_cover_pair",
    "_db_pairs",
    "_fg_chain",
    "_try_unsat_fg",
    "_try_polyunsat_fg",
    "_unsat_or_sat",
]


def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, **kw}


def _is_open_sat(info: dict) -> bool:
    return not (info.get("has_alkene") or info.get("has_alkyne"))


def _hetero_open_chain(mol: Mol, idx: int) -> bool:
    return not mol.GetAtomWithIdx(idx).IsInRing()


def _side_carbons(mol: Mol, start: int, forbid: set[int]) -> set[int]:
    seen, stack = set(), [start]
    while stack:
        cur = stack.pop()
        if cur not in seen and cur not in forbid:
            seen.add(cur)
            stack.extend(_carbon_neighbors(mol, cur))
    return seen


def _arm_ok(mol: Mol, arm: list[int], n_idx: int) -> bool:
    if not arm or len(arm) > 4 or mol.GetAtomWithIdx(arm[0]).GetIsAromatic():
        return False
    return len(_side_carbons(mol, arm[0], {n_idx})) == len(arm)


def _is_mono_fg(info: dict, flag: str, key: str) -> bool:
    xs = info.get(key) or []
    return bool(info.get(flag)) and len(xs) == 1


def _is_mono_alkene(info: dict) -> bool:
    return bool(info.get("has_alkene")) and len(info.get("double_bonds") or []) == 1


def _is_mono_alkyne(info: dict) -> bool:
    return len(info.get("triple_bonds") or []) == 1 and not (info.get("double_bonds") or [])


def _covers(chain: list[int], atoms: list[int]) -> bool:
    s = set(chain)
    return all(a in s for a in atoms)


def _best_cover_pair(mol: Mol, atoms: list[int]) -> list[int]:
    best: list[int] = []
    for i, a in enumerate(atoms):
        for b in atoms[i + 1 :]:
            chain = _chain_through_two(mol, a, b)
            if _covers(chain, atoms) and _better(mol, chain, best):
                best = chain
    return best


def _db_pairs(info: dict) -> list[tuple[int, int]]:
    return [(db["c1"], db["c2"]) for db in info.get("double_bonds") or []]


def _fg_chain(info: dict, ekey: str, kind: str, ckey: str, **extra) -> dict:
    c = info[ekey][0]["c_idx"]
    return _parent_dict(_chain_through(info, c), kind, **{ckey: c, **extra})


def _open_chain_unsat_atoms(mol, fg_c: int, db: dict) -> bool:
    """True iff principal FG attach carbon and both C=C ends are acyclic."""
    atoms = (int(fg_c), int(db["c1"]), int(db["c2"]))
    return all(not mol.GetAtomWithIdx(a).IsInRing() for a in atoms)


def _ok_unsat_fg(info: dict, flag: str, ekey: str, bad: tuple) -> bool:
    """Mono open-chain unsat FG: FG carbon + C=C ends not in ring (aryl OK)."""
    if info.get("has_alkyne") or not _is_mono_fg(info, flag, ekey) or not _is_mono_alkene(info):
        return False
    c, db, mol = info[ekey][0]["c_idx"], info["double_bonds"][0], info["mol"]
    return _open_chain_unsat_atoms(mol, c, db) and _no_fgs(info, bad)


def _ok_polyunsat_fg(info: dict, flag: str, ekey: str, bad: tuple) -> bool:
    """Multi open-chain unsat FG: FG + 2+ C=C ends all acyclic."""
    if info.get("has_alkyne") or not _is_mono_fg(info, flag, ekey):
        return False
    dbs = info.get("double_bonds") or []
    if len(dbs) < 2:
        return False
    c, mol = info[ekey][0]["c_idx"], info["mol"]
    return all(_open_chain_unsat_atoms(mol, c, db) for db in dbs) and _no_fgs(info, bad)


def _unsat_cover_atoms(info, c_idx: int, db: dict) -> list[int]:
    atoms = {c_idx, db["c1"], db["c2"]}
    for key in ("hydroxyls", "ketones", "amines"):
        for e in _aliphatic_entries(info, key):
            atoms.add(e["c_idx"])
    return list(atoms)


def _polyunsat_cover_atoms(info, c_idx: int) -> list[int]:
    atoms = {c_idx}
    for db in info.get("double_bonds") or []:
        atoms.update((db["c1"], db["c2"]))
    return list(atoms)


def _try_unsat_fg(info, flag, ekey, bad, kind, ckey, **extra) -> dict | None:
    if not _ok_unsat_fg(info, flag, ekey, bad):
        return None
    c_idx, db = info[ekey][0]["c_idx"], info["double_bonds"][0]
    chain = _best_cover_pair(info["mol"], _unsat_cover_atoms(info, c_idx, db))
    if not chain or c_idx not in chain:
        return None
    return _parent_dict(
        chain, kind, **{ckey: c_idx, "double_bond": (db["c1"], db["c2"]), "mol": info["mol"], **extra},
    )


def _try_polyunsat_fg(info, flag, ekey, bad, kind, ckey, **extra) -> dict | None:
    if not _ok_polyunsat_fg(info, flag, ekey, bad):
        return None
    c_idx = info[ekey][0]["c_idx"]
    chain = _best_cover_pair(info["mol"], _polyunsat_cover_atoms(info, c_idx))
    if not chain or c_idx not in chain:
        return None
    return _parent_dict(
        chain, kind,
        **{ckey: c_idx, "double_bonds": _db_pairs(info), "mol": info["mol"], **extra},
    )


def _unsat_or_sat(info, flag, ekey, bad, ukind, skind, ckey, **extra):
    from namepredict.layer2.alkynoic import ynsat_or_unsat_or_sat as _y
    return _y(info, flag, ekey, bad, ukind, skind, ckey, **extra)
