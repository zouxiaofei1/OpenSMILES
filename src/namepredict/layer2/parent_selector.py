from __future__ import annotations
from rdkit.Chem import Mol

from namepredict.layer2.alkenedioic import _alkenedioic_parent, _is_simple_alkenedioic
from namepredict.layer2.arene_carbonyl import (
    _try_acetophenone_parent, _try_arene_other_fg, _try_benzaldehyde_parent,
    _try_benzoic_parent,
)
from namepredict.layer2.cyclo_carboxylic import _try_cycloalkanecarboxylic_parent
from namepredict.layer2.heteroarene5 import _try_diazine_parent, _try_hetero5_parent
from namepredict.layer2.pyridine import (
    _try_pyridine_parent, _try_pyridinecarboxylic_parent,
    _try_pyridin_fg_parent,
)
from namepredict.layer2.ring_parent import (
    _benzenediol_parent, _endocyclic_double, _is_simple_aniline,
    _is_simple_benzene, _is_simple_benzenediol, _is_simple_cycloalcohol,
    _is_simple_cycloalkane, _is_simple_cycloalkene, _is_simple_cycloamine,
    _is_simple_cycloketone, _is_simple_phenol,
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
    xs = [e["c_idx"] for e in entries or [] if "c_idx" in e]
    return xs if len(xs) == n and len(set(xs)) == n else None
def _no_fgs(info: dict, keys: tuple) -> bool:
    return not any(info.get(k) for k in keys)
_CORE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_acyl_chloride", "has_anhydride",
)
_DIOL_BAD = _CORE_BAD + ("has_amine",)
_DIACID_BAD = (
    "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_anhydride",
)
_DIAMINE_BAD = _CORE_BAD + ("has_alcohol",)
_DIONE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_amine", "has_alcohol", "has_acyl_chloride", "has_anhydride",
)
_ANHYDRIDE_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol", "has_acyl_chloride",
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
def _ring_fg_try(info: dict, pairs, ekey: str, ckey: str) -> dict | None:
    for pred, kind in pairs:
        if pred(info):
            return _cyclo_fg_parent(info, kind, ekey, ckey)
    return None
def _ring_alcohol_parent(info: dict) -> dict | None:
    if _is_simple_benzenediol(info):
        return _benzenediol_parent(info)
    return _try_pyridin_fg_parent(info) or _ring_fg_try(
        info, ((_is_simple_phenol, "phenol"), (_is_simple_cycloalcohol, "cycloalcohol")),
        "hydroxyls", "oh_c_idx",
    )
def _alcohol_parent(info: dict) -> dict:
    ring = _ring_alcohol_parent(info)
    if ring is not None: return ring
    if _is_simple_alkanetriol(info): return _triol_parent(info)
    if _is_simple_alkanediol(info): return _diol_parent(info)
    return _unsat_or_sat(
        info, "has_alcohol", "hydroxyls", _ALKENOL_BAD, "alkenol", "alcohol", "oh_c_idx",
    )
def _thiol_parent(info: dict) -> dict:
    return _fg_chain(info, "thiols", "thiol", "sh_c_idx")
def _side_carbons(mol: Mol, start: int, forbid: set[int]) -> set[int]:
    seen, stack = set(), [start]
    while stack:
        cur = stack.pop()
        if cur not in seen and cur not in forbid:
            seen.add(cur)
            stack.extend(_carbon_neighbors(mol, cur))
    return seen
def _arm_ok(mol: Mol, arm: list[int], n_idx: int) -> bool:
    return bool(arm) and len(arm) <= 4 and len(_side_carbons(mol, arm[0], {n_idx})) == len(arm)
def _amine_of_deg(info: dict, deg: int) -> dict | None:
    ams = [a for a in info.get("amines") or [] if a.get("degree") == deg]
    return ams[0] if len(ams) == 1 and len(info.get("amines") or []) == 1 else None
def _n_arms(info: dict, deg: int) -> list[list[int]] | None:
    am = _amine_of_deg(info, deg)
    if am is None: return None
    mol, n_idx, cs = info["mol"], am["n_idx"], am["c_idxs"]
    arms = [_longest_from(mol, c, set()) for c in cs]
    return arms if all(_arm_ok(mol, a, n_idx) for a in arms) else None
def _sec_amine_parent(info: dict) -> dict | None:
    if not _is_open_sat(info) or not _no_fgs(info, _DIAMINE_BAD): return None
    arms = _n_arms(info, 2)
    if arms is None: return None
    parent, n_arm = (arms[0], arms[1]) if len(arms[0]) >= len(arms[1]) else (arms[1], arms[0])
    return _parent_dict(parent, "sec_amine", amine_c_idx=parent[0], n_alkyl_n=len(n_arm))
def _tert_amine_parent(info: dict) -> dict | None:
    if not _is_open_sat(info) or not _no_fgs(info, _DIAMINE_BAD): return None
    arms = _n_arms(info, 3)
    if arms is None: return None
    arms = sorted(arms, key=len, reverse=True)
    parent, ns = arms[0], [len(a) for a in arms[1:]]
    return _parent_dict(parent, "tert_amine", amine_c_idx=parent[0], n_alkyl_ns=ns)
def _primary_amine_parent(info: dict) -> dict:
    prim = next((a for a in info.get("amines") or [] if "c_idx" in a), None)
    if prim is None: return _parent_dict(_longest_chain(info["mol"]), "alkane")
    return _parent_dict(_chain_through(info, prim["c_idx"]), "amine", amine_c_idx=prim["c_idx"])
def _ring_amine_parent(info: dict) -> dict | None:
    return _try_pyridin_fg_parent(info) or _ring_fg_try(
        info, ((_is_simple_aniline, "aniline"), (_is_simple_cycloamine, "cycloamine")),
        "amines", "amine_c_idx",
    )
def _amine_parent(info: dict) -> dict:
    ring = _ring_amine_parent(info)
    if ring is not None: return ring
    if _is_simple_alkanediamine(info): return _diamine_parent(info)
    return _tert_amine_parent(info) or _sec_amine_parent(info) or _primary_amine_parent(info)
_ETHER_BAD = _CORE_BAD + ("has_amine", "has_alcohol", "has_thiol", "has_sulfide")
_SULFIDE_BAD = _CORE_BAD + ("has_amine", "has_alcohol", "has_thiol", "has_ether")
def _ether_arms(info: dict) -> tuple[list[int], list[int], dict] | None:
    ets = info.get("ethers") or []
    if len(ets) != 1:
        return None
    e, mol = ets[0], info["mol"]
    o_idx, cs = e["o_idx"], [e["c1"], e["c2"]]
    arms = [_longest_from(mol, c, set()) for c in cs]
    return (arms[0], arms[1], e) if all(_arm_ok(mol, a, o_idx) for a in arms) else None
def _ether_parent(info: dict) -> dict | None:
    if not _is_open_sat(info) or not _no_fgs(info, _ETHER_BAD):
        return None
    got = _ether_arms(info)
    if got is None:
        return None
    a1, a2, e = got
    parent, short = (a1, a2) if len(a1) >= len(a2) else (a2, a1)
    return _parent_dict(
        parent, "ether", ether_c_idx=parent[0], o_idx=e["o_idx"], alkoxy_n=len(short),
    )
def _sulfide_arms(info: dict) -> tuple[list[int], list[int], dict] | None:
    sfs = info.get("sulfides") or []
    if len(sfs) != 1:
        return None
    e, mol = sfs[0], info["mol"]
    s_idx, cs = e["s_idx"], [e["c1"], e["c2"]]
    arms = [_longest_from(mol, c, set()) for c in cs]
    return (arms[0], arms[1], e) if all(_arm_ok(mol, a, s_idx) for a in arms) else None
def _sulfide_parent(info: dict) -> dict | None:
    if not _is_open_sat(info) or not _no_fgs(info, _SULFIDE_BAD):
        return None
    got = _sulfide_arms(info)
    if got is None:
        return None
    a1, a2, e = got
    parent = a1 if len(a1) >= len(a2) else a2
    return _parent_dict(
        parent, "sulfide", s_idx=e["s_idx"], alkyl_ns=(len(a1), len(a2)),
    )
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
    return _parent_dict(chain, kind, **{ckey: c_idx, "double_bond": (db["c1"], db["c2"]), "mol": info["mol"], **extra})
def _fg_chain(info: dict, ekey: str, kind: str, ckey: str, **extra) -> dict:
    c = info[ekey][0]["c_idx"]
    return _parent_dict(_chain_through(info, c), kind, **{ckey: c, **extra})
def _unsat_or_sat(info, flag, ekey, bad, ukind, skind, ckey, **extra):
    u = _try_unsat_fg(info, flag, ekey, bad, ukind, ckey, **extra)
    return u or _fg_chain(info, ekey, skind, ckey, **extra)
def _acid_parent(info: dict) -> dict:
    for fn in (_try_pyridinecarboxylic_parent, _try_benzoic_parent, _try_cycloalkanecarboxylic_parent):
        b = fn(info)
        if b is not None: return b
    if _is_simple_alkanedioic(info): return _diacid_parent(info)
    if _is_simple_alkenedioic(info): return _alkenedioic_parent(info)
    return _unsat_or_sat(
        info, "has_acid", "carboxyls", _ALKENOIC_BAD, "alkenoic_acid", "acid", "cooh_c_idx",
    )
def _ketone_parent(info: dict) -> dict:
    a = _try_acetophenone_parent(info)
    if a is not None:
        return a
    if _is_simple_cycloketone(info):
        return _cyclo_fg_parent(info, "cycloketone", "ketones", "ketone_c_idx")
    if _is_simple_alkanedione(info):
        return _dione_parent(info)
    return _fg_chain(info, "ketones", "ketone", "ketone_c_idx")
def _aldehyde_parent(info: dict) -> dict:
    return _try_benzaldehyde_parent(info) or _unsat_or_sat(
        info, "has_aldehyde", "aldehydes", _ALKENAL_BAD, "alkenal", "aldehyde", "aldehyde_c_idx",
    )
def _amide_n_meta(info: dict) -> dict:
    ams = info.get("amides") or []
    if len(ams) != 1: return {}
    mol, am, cs = info["mol"], ams[0], ams[0].get("n_c_idxs") or []
    arms = [_longest_from(mol, c, set()) for c in cs]
    if not all(_arm_ok(mol, a, am["n_idx"]) for a in arms): return {}
    if len(arms) == 1: return {"n_alkyl_n": len(arms[0])}
    return {"n_alkyl_ns": [len(a) for a in arms]} if arms else {}
def _amide_parent(info: dict) -> dict:
    return _fg_chain(info, "amides", "amide", "amide_c_idx", **_amide_n_meta(info))
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
    return bool(info.get("has_alkene")) and len(info.get("double_bonds") or []) == 1
def _is_mono_alkyne(info: dict) -> bool:
    return len(info.get("triple_bonds") or []) == 1 and not (info.get("double_bonds") or [])
def _no_main_fg(info: dict) -> bool:
    return _no_fgs(info, _CORE_BAD + ("has_amine", "has_alcohol"))
def _is_polyene(info: dict) -> bool:
    bonds = info.get("double_bonds") or []
    if len(bonds) < 2 or info.get("has_ring") or info.get("triple_bonds"):
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
def _benzene_parent(info: dict) -> dict:
    return _parent_dict(_ring_atoms(info), "benzene")
def _cycloalkene_parent(info: dict) -> dict:
    chain = _ring_atoms(info)
    return _parent_dict(chain, "cycloalkene", double_bond=_endocyclic_double(info, set(chain)))
def _cyclo_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    return _parent_dict(_ring_atoms(info), kind, **{ckey: info[ekey][0]["c_idx"]})
def _pick_mono(info: dict, flag: str, key: str, fn):
    return fn(info) if _is_mono_fg(info, flag, key) else None
def _chain_carbonyl_fg(info: dict) -> dict | None:
    return (
        _pick_mono(info, "has_acyl_chloride", "acyl_chlorides", _acyl_chloride_parent)
        or _pick_mono(info, "has_ester", "esters", _ester_parent)
        or _pick_mono(info, "has_amide", "amides", _amide_parent)
        or _pick_mono(info, "has_nitrile", "nitriles", _nitrile_parent)
    )
def _mono_other_carbonyl(info: dict) -> dict | None:
    if _is_sym_anhydride(info):
        return _anhydride_parent(info)
    return _try_arene_other_fg(info) or _chain_carbonyl_fg(info)
def _acid_ester_amide(info: dict) -> dict | None:
    if info.get("has_acid") and info.get("carboxyls"):
        return _acid_parent(info)
    return _mono_other_carbonyl(info)
def _aldehyde_ketone(info: dict) -> dict | None:
    if info.get("has_aldehyde") and info.get("aldehydes"):
        return _aldehyde_parent(info)
    return _ketone_parent(info) if info.get("has_ketone") and info.get("ketones") else None
def _carbonyl_parent(info: dict) -> dict | None:
    top = _acid_ester_amide(info)
    return top if top is not None else _aldehyde_ketone(info)
def _hetero_parent(info: dict) -> dict | None:
    if info.get("has_alcohol") and info.get("hydroxyls"): return _alcohol_parent(info)
    if info.get("has_thiol") and info.get("thiols"): return _thiol_parent(info)
    if info.get("has_amine") and info.get("amines"): return _amine_parent(info)
    eth = _ether_parent(info)
    return eth if eth is not None else _sulfide_parent(info)
def _fg_parent(info: dict) -> dict | None:
    carb = _carbonyl_parent(info)
    return carb if carb is not None else _hetero_parent(info)
def _unsat_parent(info: dict) -> dict | None:
    if _is_mono_alkyne(info): return _alkyne_parent(info)
    if _is_polyene(info): return _polyene_parent(info)
    return _alkene_parent(info) if _is_mono_alkene(info) else None
def _ring_parent(info: dict) -> dict | None:
    for try_fn in (_try_pyridine_parent, _try_diazine_parent, _try_hetero5_parent):
        p = try_fn(info)
        if p is not None: return p
    if _is_simple_benzene(info): return _benzene_parent(info)
    return _cycloalkane_parent(info) if _is_simple_cycloalkane(info) else None
def select_parent(info: dict) -> dict:
    fg = _fg_parent(info)
    if fg is not None: return fg
    ring = _ring_parent(info)
    if ring is not None: return ring
    unsat = _unsat_parent(info)
    if unsat is not None: return unsat
    return _parent_dict(_longest_chain(info["mol"]), "alkane")
