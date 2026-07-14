from __future__ import annotations
from rdkit.Chem import Mol
from namepredict.layer2.alkenedioic import _alkenedioic_parent, _is_simple_alkenedioic
from namepredict.layer2.arene_carbonyl import (
    _try_acetophenone_parent, _try_arene_other_fg, _try_benzaldehyde_parent, _try_benzoic_parent,
)
from namepredict.layer2.cyclo_carboxylic import _try_cycloalkanecarboxylic_parent
from namepredict.layer2.hetero5_carboxylic import _try_hetero5carboxylic_parent as _try_h5cooh
from namepredict.layer2.sat_hetero_carboxylic import _try_sat_hetero_carboxylic_parent as _try_shcooh
from namepredict.layer2.benzofuran import _try_benzofuranamine_parent as _try_bfam
from namepredict.layer2.benzothiophene import _try_benzothiophenol_parent as _try_btol
from namepredict.layer2.benzothiazole import _try_benzothiazolamine_parent as _try_btzam
from namepredict.layer2.benzoxazole import _try_benzoxazolamine_parent as _try_boxam
from namepredict.layer2.benzimidazole import _try_benzimidazolamine_parent as _try_bimam
from namepredict.layer2.indazole import (
    _try_indazolecarbaldehyde_parent as _try_izald, _try_indazolecarbonitrile_parent as _try_izcn,
)
from namepredict.layer2.quinoline import (
    _try_quinolinecarboxylic_parent as _try_qcooh, _try_quinolinol_parent as _try_qol,
)
from namepredict.layer2.heteroarene5 import _try_pyrimidinamine_parent
from namepredict.layer2.pyridine import (
    _try_pyridinecarbonitrile_parent as _try_pycn, _try_pyridinecarboxylic_parent, _try_pyridin_fg_parent,
)
from namepredict.layer2.phenol_aniline import (
    _aniline_parent, _is_simple_aniline, _is_simple_phenol, _phenol_parent)
from namepredict.layer2.ring_parent import (
    _benzenediol_parent, _endocyclic_double, _is_simple_benzene, _is_simple_benzenediol,
    _is_simple_cycloalcohol, _is_simple_cycloalkane, _is_simple_cycloalkene,
    _is_simple_cycloamine, _is_simple_cycloketone, _unsub_phenyl_at)
from namepredict.layer2.scoring import _pick_best
from namepredict.layer2.aliph_fg import _aliph_c_idxs, _aliphatic_entries, _c_idxs
from namepredict.layer2.chain_walk import (
    _better, _best_arm_away, _carbon_neighbors, _chain_through, _chain_through_bond,
    _chain_through_two, _longest_chain, _longest_from,
)
from namepredict.layer2.polyalkenol import _polyalkenol_parent as _try_polyalkenol
from namepredict.layer2.prefix_alkenoic import _prefix_alkenoic_parent as _try_hy_alkenoic
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
def _is_open_sat(info: dict) -> bool: return not (info.get("has_alkene") or info.get("has_alkyne"))
def _hetero_open_chain(mol: Mol, idx: int) -> bool:
    """True iff hetero atom is not in a ring (open-chain ether/amine/sulfide)."""
    return not mol.GetAtomWithIdx(idx).IsInRing()
def _is_simple_n(info: dict, bad: tuple, ekey: str, n: int) -> bool:
    """Poly-FG only when exactly n *aliphatic* attachment carbons."""
    return (
        _is_open_sat(info) and _no_fgs(info, bad)
        and _aliph_c_idxs(info, ekey, n) is not None
    )
def _cover_parent(info: dict, ekey: str, n: int, kind: str, key: str) -> dict:
    atoms = _aliph_c_idxs(info, ekey, n) or _c_idxs(info.get(ekey) or [], n) or []
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
    if _is_simple_benzenediol(info): return _benzenediol_parent(info)
    top = _try_qol(info) or _try_btol(info) or _try_pyridin_fg_parent(info)
    if top is not None: return top
    if _is_simple_phenol(info): return _phenol_parent(info)
    return _ring_fg_try(info, ((_is_simple_cycloalcohol, "cycloalcohol"),),
                        "hydroxyls", "oh_c_idx")
def _chain_alcohol_parent(info: dict) -> dict | None:
    """Mono aliphatic OH → alcohol parent (ignore aromatic phenolic OH)."""
    aliph = _aliphatic_entries(info, "hydroxyls")
    if len(aliph) != 1:
        return None
    return _parent_dict(
        _chain_through(info, aliph[0]["c_idx"]), "alcohol", oh_c_idx=aliph[0]["c_idx"],
    )
def _polyol_or_chain_alcohol(info: dict) -> dict | None:
    if _is_simple_alkanetriol(info): return _triol_parent(info)
    if _is_simple_alkanediol(info): return _diol_parent(info)
    return _chain_alcohol_parent(info)


def _chain_or_unsat_alcohol(info: dict) -> dict | None:
    """Mono/poly alkenol first, else saturated polyol/alcohol."""
    poly = _try_polyalkenol(
        info, _ALKENOL_BAD, _best_cover_pair, _parent_dict, _db_pairs,
    )
    if poly is not None:
        return poly
    unsat = _try_unsat_fg(
        info, "has_alcohol", "hydroxyls", _ALKENOL_BAD, "alkenol", "oh_c_idx",
    )
    return unsat if unsat is not None else _polyol_or_chain_alcohol(info)


def _alcohol_parent(info: dict) -> dict | None:
    """Ring OH; poly/mono alkenol; else polyol/saturated alcohol."""
    ring = _ring_alcohol_parent(info)
    return ring if ring is not None else _chain_or_unsat_alcohol(info)
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
    if not arm or len(arm) > 4 or mol.GetAtomWithIdx(arm[0]).GetIsAromatic():
        return False
    return len(_side_carbons(mol, arm[0], {n_idx})) == len(arm)
def _arm_has_aryl(mol: Mol, arm: list[int]) -> bool:
    """True if any arm carbon is bonded to an aromatic carbon outside the arm."""
    arm_set = set(arm)
    for i in arm:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetIdx() not in arm_set and n.GetAtomicNum() == 6 and n.GetIsAromatic():
                return True
    return False
def _pick_amine_arms(mol: Mol, arms: list[list[int]]) -> tuple[list[int], list[list[int]]]:
    """Parent = longest arm; tie-break prefer aryl-bearing arm."""
    def key(a: list[int]) -> tuple:
        return (len(a), 1 if _arm_has_aryl(mol, a) else 0)
    ranked = sorted(arms, key=key, reverse=True)
    return ranked[0], ranked[1:]
def _amine_of_deg(info: dict, deg: int) -> dict | None:
    ams = [a for a in info.get("amines") or [] if a.get("degree") == deg]
    return ams[0] if len(ams) == 1 and len(info.get("amines") or []) == 1 else None
def _n_arms(info: dict, deg: int) -> list[list[int]] | None:
    am = _amine_of_deg(info, deg)
    if am is None: return None
    mol, n_idx, cs = info["mol"], am["n_idx"], am["c_idxs"]
    if not _hetero_open_chain(mol, n_idx):
        return None
    arms = [_longest_from(mol, c, set()) for c in cs]
    return arms if all(_arm_ok(mol, a, n_idx) for a in arms) else None
def _amine_sat_ok(info: dict) -> bool: return _no_fgs(info, _CORE_BAD) and not (info.get("has_alkene") or info.get("has_alkyne"))
def _sec_amine_parent(info: dict) -> dict | None:
    if not _amine_sat_ok(info): return None
    arms = _n_arms(info, 2)
    if arms is None: return None
    mol = info["mol"]
    parent, rest = _pick_amine_arms(mol, arms)
    return _parent_dict(parent, "sec_amine", amine_c_idx=parent[0], n_alkyl_n=len(rest[0]))
def _tert_amine_parent(info: dict) -> dict | None:
    if not _amine_sat_ok(info): return None
    arms = _n_arms(info, 3)
    if arms is None: return None
    parent, rest = _pick_amine_arms(info["mol"], arms)
    return _parent_dict(parent, "tert_amine", amine_c_idx=parent[0], n_alkyl_ns=[len(a) for a in rest])
def _primary_amine_parent(info: dict) -> dict:
    aliph = _aliphatic_entries(info, "amines")
    prim = next((a for a in aliph if "c_idx" in a), None)
    if prim is None:
        prim = next((a for a in info.get("amines") or [] if "c_idx" in a), None)
    if prim is None: return _parent_dict(_longest_chain(info["mol"]), "alkane")
    return _parent_dict(_chain_through(info, prim["c_idx"]), "amine", amine_c_idx=prim["c_idx"])
def _ring_amine_parent(info: dict) -> dict | None:
    top = (_try_btzam(info) or _try_boxam(info) or _try_bimam(info) or _try_bfam(info)
           or _try_pyrimidinamine_parent(info) or _try_pyridin_fg_parent(info))
    if top is not None: return top
    if _is_simple_aniline(info): return _aniline_parent(info)
    return _ring_fg_try(info, ((_is_simple_cycloamine, "cycloamine"),), "amines", "amine_c_idx")
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
    if not _hetero_open_chain(mol, o_idx):
        return None
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
    if not _hetero_open_chain(mol, s_idx):
        return None
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
# hydroxy/amino/oxo are prefixes on alkenoic acids (P-65.1.2), not competing
_ALKENOIC_OK = frozenset({"has_alcohol", "has_amine", "has_ketone"})
_ALKENOIC_BAD = tuple(k for k in _DIACID_BAD + ("has_thiol",) if k not in _ALKENOIC_OK)
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
def _unsat_cover_atoms(info, c_idx: int, db: dict) -> list[int]:
    atoms = {c_idx, db["c1"], db["c2"]}
    for key in ("hydroxyls", "ketones", "amines"):
        for e in _aliphatic_entries(info, key):
            atoms.add(e["c_idx"])
    return list(atoms)
def _try_unsat_fg(info, flag, ekey, bad, kind, ckey, **extra) -> dict | None:
    if not _ok_unsat_fg(info, flag, ekey, bad):
        return None
    c_idx, db = info[ekey][0]["c_idx"], info["double_bonds"][0]
    chain = _best_cover_pair(info["mol"], _unsat_cover_atoms(info, c_idx, db))
    if not chain or c_idx not in chain:
        return None
    return _parent_dict(chain, kind, **{ckey: c_idx, "double_bond": (db["c1"], db["c2"]), "mol": info["mol"], **extra})
def _fg_chain(info: dict, ekey: str, kind: str, ckey: str, **extra) -> dict:
    c = info[ekey][0]["c_idx"]
    return _parent_dict(_chain_through(info, c), kind, **{ckey: c, **extra})
def _unsat_or_sat(info, flag, ekey, bad, ukind, skind, ckey, **extra):
    u = _try_unsat_fg(info, flag, ekey, bad, ukind, ckey, **extra)
    return u or _fg_chain(info, ekey, skind, ckey, **extra)
def _with_anion(info: dict, parent: dict) -> dict:
    if any(c.get("anion") for c in info.get("carboxyls") or []):
        return {**parent, "anion": True}
    return parent
def _ring_acid_try(info: dict) -> dict | None:
    for fn in (_try_shcooh, _try_h5cooh, _try_qcooh, _try_pyridinecarboxylic_parent,
               _try_benzoic_parent, _try_cycloalkanecarboxylic_parent):
        if (b := fn(info)) is not None: return b
    return None
def _poly_acid_try(info: dict) -> dict | None:
    if _is_simple_alkanedioic(info): return _diacid_parent(info)
    if _is_simple_alkenedioic(info): return _alkenedioic_parent(info)
    return _try_hy_alkenoic(
        info, _ALKENOIC_BAD, _best_cover_pair, _parent_dict, _db_pairs,
    )
def _acid_parent_core(info: dict) -> dict:
    top = _ring_acid_try(info) or _poly_acid_try(info)
    if top is not None:
        return top
    return _unsat_or_sat(
        info, "has_acid", "carboxyls", _ALKENOIC_BAD, "alkenoic_acid", "acid", "cooh_c_idx",
    )
def _acid_parent(info: dict) -> dict:
    return _with_anion(info, _acid_parent_core(info))
def _ketone_parent(info: dict) -> dict:
    a = _try_acetophenone_parent(info)
    if a is not None: return a
    if _is_simple_cycloketone(info):
        return _cyclo_fg_parent(info, "cycloketone", "ketones", "ketone_c_idx")
    if _is_simple_alkanedione(info):
        return _dione_parent(info)
    return _fg_chain(info, "ketones", "ketone", "ketone_c_idx")
def _aldehyde_parent(info: dict) -> dict:
    return _try_izald(info) or _try_benzaldehyde_parent(info) or _unsat_or_sat(
        info, "has_aldehyde", "aldehydes", _ALKENAL_BAD, "alkenal", "aldehyde", "aldehyde_c_idx",
    )
def _amide_n_alkyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    arms = [_longest_from(mol, c, set()) for c in cs]
    if not all(_arm_ok(mol, a, am["n_idx"]) for a in arms): return {}
    if len(arms) == 1: return {"n_alkyl_n": len(arms[0])}
    return {"n_alkyl_ns": [len(a) for a in arms]} if arms else {}
def _amide_n_phenyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    if len(cs) != 1: return {}
    return {"n_phenyl": True} if _unsub_phenyl_at(mol, cs[0], am["n_idx"]) else {}
def _amide_n_benzyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    """N–CH2–Ph (Ph may carry simple leaves) → n_benzyl flag + ch2 idx."""
    from namepredict.layer2.aryl_sub import _ch2_ph_at
    if len(cs) != 1: return {}
    ch2 = cs[0]
    ph = _ch2_ph_at(mol, ch2, am["n_idx"])
    return {"n_benzyl": True, "n_benzyl_ch2": ch2} if ph is not None else {}
def _amide_n_meta(info: dict) -> dict:
    ams = info.get("amides") or []
    if len(ams) != 1: return {}
    mol, am, cs = info["mol"], ams[0], ams[0].get("n_c_idxs") or []
    return (
        _amide_n_benzyl(mol, am, cs)
        or _amide_n_phenyl(mol, am, cs)
        or _amide_n_alkyl(mol, am, cs)
    )
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
    return _try_izcn(info) or _try_pycn(info) or _unsat_or_sat(
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
    from namepredict.layer2.ring_parent import _pick_cycloalkane_ring as _pcr
    return _parent_dict(list(_pcr(info) or _ring_atoms(info)), "cycloalkane")
def _benzene_parent(info: dict) -> dict:
    from namepredict.layer2.benzene_pick import _pick_benzene_ring as _pbr
    return _parent_dict(_pbr(info) or _ring_atoms(info), "benzene")
def _cycloalkene_parent(info: dict) -> dict:
    c = _ring_atoms(info); return _parent_dict(c, "cycloalkene", double_bond=_endocyclic_double(info, set(c)))
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
    from namepredict.layer2.candidates import _fg_parent as _fg_best
    return _fg_best(info)
def _unsat_parent(info: dict) -> dict | None:
    if _is_mono_alkyne(info): return _alkyne_parent(info)
    if _is_polyene(info): return _polyene_parent(info)
    return _alkene_parent(info) if _is_mono_alkene(info) else None
def _ring_parent(info: dict) -> dict | None:
    from namepredict.layer2.candidates import _ring_parent as _ring_best
    return _ring_best(info)
def select_parent(info: dict) -> dict:
    from namepredict.layer2.candidates import _alkane_fallback, _collect_candidates
    best = _pick_best(info, _collect_candidates(info))
    parent = best if best is not None else _alkane_fallback(info)
    return parent if parent.get("mol") is not None else {**parent, "mol": info.get("mol")}
