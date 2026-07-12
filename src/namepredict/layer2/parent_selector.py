from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _endocyclic_double,
    _is_simple_cycloalcohol,
    _is_simple_cycloalkane,
    _is_simple_cycloalkene,
    _is_simple_cycloamine,
    _is_simple_cycloketone,
)

def _carbon_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]

def _extend_best(mol: Mol, node: int, path: list[int], forbid: set[int], best: list[int]) -> list[int]:
    for nb in _carbon_neighbors(mol, node):
        if nb in path or nb in forbid:
            continue
        cand = _dfs_path(mol, nb, path + [nb], forbid)
        if len(cand) > len(best):
            best = cand
    return best

def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int]) -> list[int]:
    return _extend_best(mol, node, path, forbid, path)

def _longest_from(mol: Mol, start: int, forbidden: set[int] | None = None) -> list[int]:
    return _dfs_path(mol, start, [start], forbidden or set())

def _all_carbons(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]

def _side_count(mol: Mol, chain: list[int]) -> int:
    chain_set = set(chain)
    n = 0
    for c in chain:
        atom = mol.GetAtomWithIdx(c)
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in chain_set:
                n += 1
    return n

def _chain_key(mol: Mol, path: list[int]) -> tuple:
    return (len(path), _side_count(mol, path))

def _better(mol: Mol, cand: list[int], best: list[int]) -> bool:
    if not best:
        return True
    return _chain_key(mol, cand) > _chain_key(mol, best)

def _best_among(mol: Mol, seeds: list[int]) -> list[int]:
    best: list[int] = []
    for c in seeds:
        path = _longest_from(mol, c)
        if _better(mol, path, best):
            best = path
    return best

def _longest_chain(mol: Mol, seeds: list[int] | None = None) -> list[int]:
    return _best_among(mol, seeds or _all_carbons(mol))

def _arms_from(mol: Mol, center: int) -> list[list[int]]:
    forbid = {center}
    return [_longest_from(mol, nb, forbid) for nb in _carbon_neighbors(mol, center)]

def _join_through(center: int, arms: list[list[int]]) -> list[int]:
    arms = sorted(arms, key=len, reverse=True)
    if not arms:
        return [center]
    if len(arms) == 1:
        return list(reversed(arms[0])) + [center]
    return list(reversed(arms[0])) + [center] + arms[1]

def _chain_through(info: dict, c_idx: int) -> list[int]:
    mol: Mol = info["mol"]
    return _join_through(c_idx, _arms_from(mol, c_idx))

def _bfs_expand(mol: Mol, cur: int, prev: dict, q: list) -> None:
    for nb in _carbon_neighbors(mol, cur):
        if nb not in prev:
            prev[nb] = cur
            q.append(nb)

def _bfs_prev(mol: Mol, start: int, goal: int) -> dict | None:
    prev: dict = {start: None}
    q = [start]
    while q:
        cur = q.pop(0)
        if cur == goal:
            return prev
        _bfs_expand(mol, cur, prev, q)
    return None

def _rebuild_path(prev: dict, end: int) -> list[int]:
    path = [end]
    while prev[path[-1]] is not None:
        path.append(prev[path[-1]])
    return list(reversed(path))

def _path_between(mol: Mol, a: int, b: int) -> list[int]:
    if a == b:
        return [a]
    prev = _bfs_prev(mol, a, b)
    return _rebuild_path(prev, b) if prev else [a]

def _chain_through_two(mol: Mol, c1: int, c2: int) -> list[int]:
    path = _path_between(mol, c1, c2)
    left = _best_arm_away(mol, path[0], set(path[1:]))
    right = _best_arm_away(mol, path[-1], set(path[:-1]))
    return list(reversed(left)) + path + right

def _c_idxs(entries, n: int) -> list[int] | None:
    if not entries or len(entries) != n:
        return None
    xs = [e["c_idx"] for e in entries]
    return xs if len(set(xs)) == n else None

def _no_fgs(info: dict, keys: tuple) -> bool:
    return not any(info.get(k) for k in keys)

_DIOL_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_amine", "has_acyl_chloride",
    "has_anhydride",
)
_DIACID_BAD = (
    "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol",
    "has_anhydride",
)
_DIAMINE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_alcohol", "has_acyl_chloride",
    "has_anhydride",
)
_DIONE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_amine", "has_alcohol", "has_acyl_chloride",
    "has_anhydride",
)
_ANHYDRIDE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol",
    "has_acyl_chloride",
)

def _is_open_sat(info: dict) -> bool:
    return not (info.get("has_ring") or info.get("has_alkene") or info.get("has_alkyne"))

def _is_simple_n(info: dict, bad: tuple, ekey: str, n: int) -> bool:
    return _is_open_sat(info) and _no_fgs(info, bad) and _c_idxs(info.get(ekey) or [], n) is not None

def _cover_parent(info: dict, ekey: str, n: int, kind: str, key: str) -> dict:
    atoms = _c_idxs(info.get(ekey) or [], n) or []
    chain = _best_cover_pair(info["mol"], atoms) or _longest_chain(info["mol"])
    return _parent_dict(chain, kind, **{key: atoms})

def _is_simple_alkanediol(info: dict) -> bool:
    return _is_simple_n(info, _DIOL_BAD, "hydroxyls", 2)

def _diol_parent(info: dict) -> dict:
    return _cover_parent(info, "hydroxyls", 2, "diol", "oh_c_idxs")

def _is_simple_alkanetriol(info: dict) -> bool:
    return _is_simple_n(info, _DIOL_BAD, "hydroxyls", 3)

def _triol_parent(info: dict) -> dict:
    return _cover_parent(info, "hydroxyls", 3, "triol", "oh_c_idxs")

def _is_simple_alkanedioic(info: dict) -> bool:
    return _is_simple_n(info, _DIACID_BAD, "carboxyls", 2)

def _diacid_parent(info: dict) -> dict:
    return _cover_parent(info, "carboxyls", 2, "diacid", "cooh_c_idxs")

def _is_simple_alkanediamine(info: dict) -> bool:
    return _is_simple_n(info, _DIAMINE_BAD, "amines", 2)

def _diamine_parent(info: dict) -> dict:
    return _cover_parent(info, "amines", 2, "diamine", "amine_c_idxs")

def _is_simple_alkanedione(info: dict) -> bool:
    return _is_simple_n(info, _DIONE_BAD, "ketones", 2)

def _dione_parent(info: dict) -> dict:
    return _cover_parent(info, "ketones", 2, "dione", "ketone_c_idxs")

def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, **kw}

def _alcohol_parent(info: dict) -> dict:
    if _is_simple_cycloalcohol(info):
        return _cycloalcohol_parent(info)
    if _is_simple_alkanetriol(info):
        return _triol_parent(info)
    if _is_simple_alkanediol(info):
        return _diol_parent(info)
    return _unsat_or_sat(
        info, "has_alcohol", "hydroxyls", _ALKENOL_BAD, "alkenol", "alcohol", "oh_c_idx",
    )

def _thiol_parent(info: dict) -> dict:
    return _fg_chain(info, "thiols", "thiol", "sh_c_idx")

def _amine_parent(info: dict) -> dict:
    if _is_simple_cycloamine(info):
        return _cycloamine_parent(info)
    if _is_simple_alkanediamine(info):
        return _diamine_parent(info)
    return _fg_chain(info, "amines", "amine", "amine_c_idx")

_ALKENOIC_BAD = _DIACID_BAD + ("has_thiol",)
_UNSAT_FG_BASE = (
    "has_acid", "has_ester", "has_amide", "has_ketone", "has_amine",
    "has_acyl_chloride", "has_anhydride", "has_thiol",
)
_UNSAT_FG_CORE = _UNSAT_FG_BASE + ("has_alcohol",)
_ALKENAL_BAD = _UNSAT_FG_CORE + ("has_nitrile",)
_ALKENENITRILE_BAD = _UNSAT_FG_CORE + ("has_aldehyde",)
_ALKENOL_BAD = _UNSAT_FG_BASE + ("has_aldehyde", "has_nitrile")
_ALKENOATE_BAD = tuple(k for k in _UNSAT_FG_CORE + ("has_nitrile",) if k != "has_ester")

def _ok_unsat_fg(info: dict, flag: str, ekey: str, bad: tuple) -> bool:
    if info.get("has_ring") or info.get("has_alkyne"):
        return False
    return _is_mono_fg(info, flag, ekey) and _is_mono_alkene(info) and _no_fgs(info, bad)

def _try_unsat_fg(info, flag, ekey, bad, kind, ckey, **extra) -> dict | None:
    if not _ok_unsat_fg(info, flag, ekey, bad):
        return None
    c_idx, db = info[ekey][0]["c_idx"], info["double_bonds"][0]
    chain = _best_cover_pair(info["mol"], [c_idx, db["c1"], db["c2"]])
    if not chain or c_idx not in chain:
        return None
    return _parent_dict(chain, kind, **{ckey: c_idx, "double_bond": (db["c1"], db["c2"]), **extra})

def _fg_chain(info: dict, ekey: str, kind: str, ckey: str, **extra) -> dict:
    c = info[ekey][0]["c_idx"]
    return _parent_dict(_chain_through(info, c), kind, **{ckey: c, **extra})

def _unsat_or_sat(info, flag, ekey, bad, ukind, skind, ckey, **extra):
    u = _try_unsat_fg(info, flag, ekey, bad, ukind, ckey, **extra)
    return u or _fg_chain(info, ekey, skind, ckey, **extra)

def _acid_parent(info: dict) -> dict:
    if _is_simple_alkanedioic(info):
        return _diacid_parent(info)
    return _unsat_or_sat(
        info, "has_acid", "carboxyls", _ALKENOIC_BAD, "alkenoic_acid", "acid", "cooh_c_idx",
    )

def _ketone_parent(info: dict) -> dict:
    if _is_simple_cycloketone(info):
        return _cycloketone_parent(info)
    if _is_simple_alkanedione(info):
        return _dione_parent(info)
    return _fg_chain(info, "ketones", "ketone", "ketone_c_idx")

def _aldehyde_parent(info: dict) -> dict:
    return _unsat_or_sat(
        info, "has_aldehyde", "aldehydes", _ALKENAL_BAD, "alkenal", "aldehyde", "aldehyde_c_idx",
    )

def _amide_parent(info: dict) -> dict:
    return _fg_chain(info, "amides", "amide", "amide_c_idx")

def _is_mono_fg(info: dict, flag: str, key: str) -> bool:
    xs = info.get(key) or []
    return bool(info.get(flag)) and len(xs) == 1

def _acyl_chloride_parent(info: dict) -> dict:
    e = info["acyl_chlorides"][0]
    return _parent_dict(
        _chain_through(info, e["c_idx"]), "acyl_chloride",
        acyl_c_idx=e["c_idx"], cl_idx=e["cl_idx"],
    )

def _acyl_chain_len(info: dict, c_idx: int) -> int:
    return len(_chain_through(info, c_idx))

def _is_sym_anhydride(info: dict) -> bool:
    if not _is_open_sat(info) or not _no_fgs(info, _ANHYDRIDE_BAD):
        return False
    anhs = info.get("anhydrides") or []
    if len(anhs) != 1:
        return False
    e = anhs[0]
    return _acyl_chain_len(info, e["c1_idx"]) == _acyl_chain_len(info, e["c2_idx"])

def _anhydride_parent(info: dict) -> dict:
    e = info["anhydrides"][0]
    chain = _chain_through(info, e["c1_idx"])
    return _parent_dict(
        chain, "anhydride",
        acyl_c_idx=e["c1_idx"], o_idx=e["o_idx"], other_acyl_c_idx=e["c2_idx"],
    )

def _nitrile_parent(info: dict) -> dict:
    return _unsat_or_sat(
        info, "has_nitrile", "nitriles", _ALKENENITRILE_BAD,
        "alkenenitrile", "nitrile", "nitrile_c_idx",
    )

def _ester_meta(info: dict) -> dict:
    e = info["esters"][0]
    return dict(
        o_idx=e["o_idx"], alkoxy_c_idx=e["alkoxy_c_idx"],
        alkoxy_n=len(_longest_from(info["mol"], e["alkoxy_c_idx"])),
    )

def _ester_parent(info: dict) -> dict:
    return _unsat_or_sat(
        info, "has_ester", "esters", _ALKENOATE_BAD, "alkenoate", "ester",
        "ester_c_idx", **_ester_meta(info),
    )

def _best_arm_away(mol: Mol, from_c: int, forbid: set[int]) -> list[int]:
    best: list[int] = []
    for nb in _carbon_neighbors(mol, from_c):
        if nb in forbid:
            continue
        path = _longest_from(mol, nb, forbid | {from_c})
        if len(path) > len(best):
            best = path
    return best

def _chain_through_bond(mol: Mol, c1: int, c2: int) -> list[int]:
    left = _best_arm_away(mol, c1, {c2})
    right = _best_arm_away(mol, c2, {c1})
    return list(reversed(left)) + [c1, c2] + right

def _alkene_parent(info: dict) -> dict:
    if _is_simple_cycloalkene(info):
        return _cycloalkene_parent(info)
    db = info["double_bonds"][0]
    c1, c2 = db["c1"], db["c2"]
    chain = _chain_through_bond(info["mol"], c1, c2)
    return _parent_dict(chain, "alkene", double_bond=(c1, c2))

def _alkyne_parent(info: dict) -> dict:
    tb = info["triple_bonds"][0]
    c1, c2 = tb["c1"], tb["c2"]
    chain = _chain_through_bond(info["mol"], c1, c2)
    return _parent_dict(chain, "alkyne", triple_bond=(c1, c2))

def _is_mono_alkene(info: dict) -> bool:
    bonds = info.get("double_bonds") or []
    return bool(info.get("has_alkene")) and len(bonds) == 1

def _is_mono_alkyne(info: dict) -> bool:
    triples = info.get("triple_bonds") or []
    doubles = info.get("double_bonds") or []
    return len(triples) == 1 and len(doubles) == 0

def _no_main_fg(info: dict) -> bool:
    bad = (
        "has_acid", "has_ester", "has_amide", "has_nitrile",
        "has_aldehyde", "has_ketone", "has_amine", "has_alcohol",
        "has_acyl_chloride", "has_anhydride",
    )
    return not any(info.get(k) for k in bad)

def _is_polyene(info: dict) -> bool:
    bonds = info.get("double_bonds") or []
    if len(bonds) < 2 or info.get("has_ring"):
        return False
    if info.get("triple_bonds"):
        return False
    return _no_main_fg(info)

def _db_atoms(info: dict) -> list[int]:
    atoms: set[int] = set()
    for db in info.get("double_bonds") or []:
        atoms.add(db["c1"])
        atoms.add(db["c2"])
    return list(atoms)

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

def _polyene_chain(info: dict) -> list[int]:
    mol, atoms = info["mol"], _db_atoms(info)
    chain = _longest_chain(mol)
    if _covers(chain, atoms):
        return chain
    return _best_cover_pair(mol, atoms) or chain

def _db_pairs(info: dict) -> list[tuple[int, int]]:
    return [(db["c1"], db["c2"]) for db in info.get("double_bonds") or []]

def _polyene_parent(info: dict) -> dict:
    chain = _polyene_chain(info)
    return _parent_dict(chain, "polyene", double_bonds=_db_pairs(info))

def _ring_atoms(info: dict) -> list[int]:
    return list(info["rings"][0]["atom_ids"])

def _cycloalkane_parent(info: dict) -> dict:
    return _parent_dict(_ring_atoms(info), "cycloalkane")

def _cycloalkene_parent(info: dict) -> dict:
    chain = _ring_atoms(info)
    return _parent_dict(chain, "cycloalkene", double_bond=_endocyclic_double(info, set(chain)))

def _cyclo_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    return _parent_dict(_ring_atoms(info), kind, **{ckey: info[ekey][0]["c_idx"]})

def _cycloalcohol_parent(info: dict) -> dict:
    return _cyclo_fg_parent(info, "cycloalcohol", "hydroxyls", "oh_c_idx")

def _cycloamine_parent(info: dict) -> dict:
    return _cyclo_fg_parent(info, "cycloamine", "amines", "amine_c_idx")

def _cycloketone_parent(info: dict) -> dict:
    return _cyclo_fg_parent(info, "cycloketone", "ketones", "ketone_c_idx")

def _mono_acyl_or_ester(info: dict) -> dict | None:
    if _is_mono_fg(info, "has_acyl_chloride", "acyl_chlorides"):
        return _acyl_chloride_parent(info)
    if _is_mono_fg(info, "has_ester", "esters"):
        return _ester_parent(info)
    return None


def _mono_amide_or_nitrile(info: dict) -> dict | None:
    if _is_mono_fg(info, "has_amide", "amides"):
        return _amide_parent(info)
    if _is_mono_fg(info, "has_nitrile", "nitriles"):
        return _nitrile_parent(info)
    return None


def _mono_anhydride(info: dict) -> dict | None:
    if _is_sym_anhydride(info):
        return _anhydride_parent(info)
    return None


def _acid_ester_amide(info: dict) -> dict | None:
    if info.get("has_acid") and info.get("carboxyls"):
        return _acid_parent(info)
    top = _mono_anhydride(info) or _mono_acyl_or_ester(info)
    return top or _mono_amide_or_nitrile(info)

def _aldehyde_ketone(info: dict) -> dict | None:
    if info.get("has_aldehyde") and info.get("aldehydes"):
        return _aldehyde_parent(info)
    if info.get("has_ketone") and info.get("ketones"):
        return _ketone_parent(info)
    return None

def _carbonyl_parent(info: dict) -> dict | None:
    top = _acid_ester_amide(info)
    return top if top is not None else _aldehyde_ketone(info)

def _hetero_parent(info: dict) -> dict | None:
    if info.get("has_alcohol") and info.get("hydroxyls"):
        return _alcohol_parent(info)
    if info.get("has_thiol") and info.get("thiols"):
        return _thiol_parent(info)
    if info.get("has_amine") and info.get("amines"):
        return _amine_parent(info)
    return None

def _fg_parent(info: dict) -> dict | None:
    carb = _carbonyl_parent(info)
    return carb if carb is not None else _hetero_parent(info)

def _unsat_parent(info: dict) -> dict | None:
    if _is_mono_alkyne(info):
        return _alkyne_parent(info)
    if _is_polyene(info):
        return _polyene_parent(info)
    if _is_mono_alkene(info):
        return _alkene_parent(info)
    return None

def select_parent(info: dict) -> dict:
    fg = _fg_parent(info)
    if fg is not None:
        return fg
    if _is_simple_cycloalkane(info):
        return _cycloalkane_parent(info)
    unsat = _unsat_parent(info)
    if unsat is not None:
        return unsat
    return _parent_dict(_longest_chain(info["mol"]), "alkane")
