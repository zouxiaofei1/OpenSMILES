from __future__ import annotations
from rdkit.Chem import Mol
from namepredict.layer2.alkenedioic import _alkenedioic_parent, _is_simple_alkenedioic
from namepredict.layer2.arene_carbonyl import (
    _try_acetophenone_parent, _try_arene_other_fg, _try_benzaldehyde_parent,
    _try_benzene_polycarboxylic, _try_benzoic_parent,
)
from namepredict.layer2.arene_fg_parent import (
    _try_arene_fg_aldehyde, _try_arene_fg_amine, _try_arene_fg_nitrile,
    _try_naphthalenediol_parent, _try_naphthalenol_parent,
    _try_quinolinediol_parent,
)
from namepredict.layer2.cyclo_carboxylic import _try_cycloalkanecarboxylic_parent
from namepredict.layer2.hetero5_carboxylic import _try_hetero5carboxylic_parent as _try_h5cooh
from namepredict.layer2.sat_hetero_carboxylic import _try_sat_hetero_carboxylic_parent as _try_shcooh
from namepredict.layer2.scaffold.benzofuran import _try_benzofuranamine_parent as _try_bfam
from namepredict.layer2.scaffold.benzothiophene import _try_benzothiophenol_parent as _try_btol
from namepredict.layer2.scaffold.benzothiazole import _try_benzothiazolamine_parent as _try_btzam
from namepredict.layer2.scaffold.benzoxazole import _try_benzoxazolamine_parent as _try_boxam
from namepredict.layer2.scaffold.benzimidazole import _try_benzimidazolamine_parent as _try_bimam
from namepredict.layer2.scaffold.indazole import (
    _try_indazolecarbaldehyde_parent as _try_izald, _try_indazolecarbonitrile_parent as _try_izcn,
)
from namepredict.layer2.scaffold.indole import _try_indolecarboxylic_parent as _try_indcooh
from namepredict.layer2.scaffold.naphthalene import _try_naphthalenecarboxylic_parent as _try_naphcooh
from namepredict.layer2.scaffold.quinoline import _try_quinolinecarboxylic_parent as _try_qcooh, _try_quinolinol_parent as _try_qol
from namepredict.layer2.scaffold.heteroarene5 import _try_pyrimidinamine_parent
from namepredict.layer2.scaffold.pyridine import (
    _try_pyridinecarbonitrile_parent as _try_pycn, _try_pyridinecarboxylic_parent, _try_pyridin_fg_parent)
from namepredict.layer2.phenol_aniline import (
    _aniline_parent, _is_simple_aniline, _is_simple_phenol, _phenol_parent)
from namepredict.layer2.cyclo_poly_fg import try_cycloalkanediol, try_cycloalkanedione
from namepredict.layer2.cyclo_ene_fg import (
    _is_simple_cycloalkenol, _is_simple_cycloalkenone,
)
from namepredict.layer2.scaffold.ring_parent import (
    _benzenediol_parent, _endocyclic_double, _is_simple_benzene, _is_simple_benzenediol,
    _is_simple_cycloalcohol, _is_simple_cycloalkane, _is_simple_cycloalkene,
    _is_simple_cycloamine, _is_simple_cycloketone)
from namepredict.layer2.aliph_fg import _aliph_c_idxs, _aliphatic_entries, _c_idxs
from namepredict.layer2.chain_walk import (
    _carbon_neighbors, _chain_through, _chain_through_bond,
    _longest_chain, _longest_from,
)
from namepredict.layer2.polyalkenol import _polyalkenol_parent as _try_polyalkenol
from namepredict.layer2.prefix_alkenoic import _prefix_alkenoic_parent as _try_hy_alkenoic
from namepredict.layer2.polycarboxylic import try_polycarboxylic_parent
from namepredict.layer2.carboxymethyl_diacid import carboxymethyl_diacid_parent
from namepredict.constants import Br, C, Cl
from namepredict.layer2.parent_selector_common import _ANHYDRIDE_BAD, _CORE_BAD, _DIACID_BAD, _DIAMINE_BAD, _DIOL_BAD, _DIONE_BAD, _no_fgs
from namepredict.layer2.parent_core import (
    _arm_ok, _best_cover_pair, _covers, _db_pairs, _fg_chain, _hetero_open_chain,
    _is_mono_alkene, _is_mono_alkyne, _is_mono_fg, _is_open_sat, _parent_dict,
    _try_unsat_fg, _unsat_or_sat,
)
def _is_simple_n(info: dict, bad: tuple, ekey: str, n: int) -> bool:
    return (
        _is_open_sat(info) and _no_fgs(info, bad)
        and _aliph_c_idxs(info, ekey, n) is not None
    )
def _cover_parent(info: dict, ekey: str, n: int, kind: str, key: str) -> dict:
    atoms = _aliph_c_idxs(info, ekey, n) or _c_idxs(info.get(ekey) or [], n) or []
    chain = _best_cover_pair(info["mol"], atoms) or _longest_chain(info["mol"])
    return _parent_dict(chain, kind, **{key: atoms})
def _polyol_kind(n: int) -> str:
    return {1: "alcohol", 2: "diol", 3: "triol"}.get(n, "polyol")
def _polyol_parent(info: dict, n: int) -> dict | None:
    if not _is_simple_n(info, _DIOL_BAD, "hydroxyls", n):
        return None
    p = _cover_parent(info, "hydroxyls", n, _polyol_kind(n), "oh_c_idxs")
    p["n_oh"] = n
    return p
def _is_simple_alkanedioic(info: dict) -> bool:
    return _is_simple_n(info, _DIACID_BAD, "carboxyls", 2)
def _diacid_parent(info: dict) -> dict:
    return _cover_parent(info, "carboxyls", 2, "diacid", "cooh_c_idxs")
def _polyamine_kind(n: int) -> str:
    return {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"}.get(n, "polyamine")
def _polyamine_parent(info: dict, n: int) -> dict | None:
    if not _is_simple_n(info, _DIAMINE_BAD, "amines", n):
        return None
    p = _cover_parent(info, "amines", n, _polyamine_kind(n), "amine_c_idxs")
    p["n_amine"] = n
    return p
def _is_simple_alkanedione(info: dict) -> bool:
    return _is_simple_n(info, _DIONE_BAD, "ketones", 2)
def _dione_parent(info: dict) -> dict:
    return _cover_parent(info, "ketones", 2, "dione", "ketone_c_idxs")
def _has_n_oh(info: dict, n: int) -> bool:
    ohs = _aliphatic_entries(info, "hydroxyls"); return len(ohs) == n
def _alkenediol_good(info: dict, bad: tuple) -> bool:
    dbs = info.get("double_bonds") or []
    return (_has_n_oh(info, 2) and len(dbs) == 1
            and not info.get("has_alkyne") and _no_fgs(info, bad))
def _polyenediol_good(info: dict, bad: tuple) -> bool:
    dbs = info.get("double_bonds") or []
    return (_has_n_oh(info, 2) and len(dbs) >= 2
            and not info.get("has_alkyne") and _no_fgs(info, bad))
def _alkynediol_good(info: dict, bad: tuple) -> bool:
    tbs = info.get("triple_bonds") or []
    return (_has_n_oh(info, 2) and len(tbs) == 1
            and not info.get("has_alkene") and _no_fgs(info, bad))
def _alkenetriol_good(info: dict, bad: tuple) -> bool:
    dbs = info.get("double_bonds") or []
    return (_has_n_oh(info, 3) and len(dbs) == 1
            and not info.get("has_alkyne") and _no_fgs(info, bad))
def _uniq_atoms(*lists) -> list[int]:
    """Deduplicated atom id list preserving first-occurrence order."""
    seen: set[int] = set()
    out: list[int] = []
    for x in lists:
        for a in x:
            if a not in seen:
                seen.add(a); out.append(a)
    return out

def _ring_fg_try(info: dict, pairs, ekey: str, ckey: str) -> dict | None:
    for pred, kind in pairs:
        if pred(info):
            return _cyclo_fg_parent(info, kind, ekey, ckey)
    return None
def _try_arene_amine(info: dict) -> dict | None:
    """Generalized arene FG amine parents (NH2 on fused/monocyclic aromatic rings)."""
    return _try_arene_fg_amine(info)


def _ring_alcohol_parent(info: dict) -> dict | None:
    if _is_simple_benzenediol(info): return _benzenediol_parent(info)
    return (_try_qol(info) or _try_btol(info) or _try_pyridin_fg_parent(info)
            or _try_naphthalenol_parent(info) or _try_naphthalenediol_parent(info)
            or _try_quinolinediol_parent(info)
            or (_phenol_parent(info) if _is_simple_phenol(info) else None)
            or try_cycloalkanediol(info)
            or (_cyclo_ene_fg_parent(info, "cycloalcohol", "hydroxyls", "oh_c_idx")
                if _is_simple_cycloalkenol(info) else
                _ring_fg_try(info, ((_is_simple_cycloalcohol, "cycloalcohol"),), "hydroxyls", "oh_c_idx")))
def _oh_c_in_ring(info: dict, c_idx: int) -> bool:
    return info["mol"].GetAtomWithIdx(int(c_idx)).IsInRing()
def _chain_alcohol_parent(info: dict) -> dict | None:
    aliph = _aliphatic_entries(info, "hydroxyls")
    if len(aliph) != 1 or _oh_c_in_ring(info, aliph[0]["c_idx"]):
        return None
    c = aliph[0]["c_idx"]
    return _parent_dict(_chain_through(info, c), "alcohol", oh_c_idx=c, oh_c_idxs=[c], n_oh=1)
def _polyol_or_chain_alcohol(info: dict) -> dict | None:
    for n in (3, 2):
        p = _polyol_parent(info, n)
        if p is not None:
            return p
    return _chain_alcohol_parent(info)
def _chain_or_unsat_alcohol(info: dict) -> dict | None:
    from namepredict.layer2.alkynoic import chain_or_unsat_alcohol as _c
    b = _ALKENOL_BAD
    up = _try_unsat_polyol(info, b)
    return up if up is not None else _c(
        info, b,
        lambda: _try_polyalkenol(info, b, _best_cover_pair, _parent_dict, _db_pairs),
        lambda: _try_unsat_fg(info, "has_alcohol", "hydroxyls", b, "alcohol", "oh_c_idx"),
        lambda: _polyol_or_chain_alcohol(info),
    )
def _alcohol_parent(info: dict) -> dict | None:
    ring = _ring_alcohol_parent(info)
    return ring if ring is not None else _chain_or_unsat_alcohol(info)
def _thiol_parent(info: dict) -> dict:
    return _fg_chain(info, "thiols", "thiol", "sh_c_idx")
def _arm_has_aryl(mol: Mol, arm: list[int]) -> bool:
    """True if any arm carbon is bonded to an aromatic carbon outside the arm."""
    arm_set = set(arm)
    for i in arm:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetIdx() not in arm_set and n.GetAtomicNum() == C and n.GetIsAromatic():
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
def _amine_arm_fields(info: dict, degree: int, parent: list[int], rest: list[list[int]]) -> dict:
    from namepredict.layer2.amine_expression import build_amine_arm_facts
    amine = _amine_of_deg(info, degree)
    facts = build_amine_arm_facts(amine["n_idx"], degree, parent, rest)
    return {"amine_c_idx": parent[0], "neutral_amine_arm_facts": facts}


def _sec_amine_parent(info: dict) -> dict | None:
    if not _amine_sat_ok(info): return None
    arms = _n_arms(info, 2)
    if arms is None: return None
    parent, rest = _pick_amine_arms(info["mol"], arms)
    fields = _amine_arm_fields(info, 2, parent, rest)
    return _parent_dict(parent, "sec_amine", n_alkyl_n=len(rest[0]), **fields)
def _tert_amine_parent(info: dict) -> dict | None:
    if not _amine_sat_ok(info): return None
    arms = _n_arms(info, 3)
    if arms is None: return None
    parent, rest = _pick_amine_arms(info["mol"], arms)
    fields = _amine_arm_fields(info, 3, parent, rest)
    return _parent_dict(parent, "tert_amine", n_alkyl_ns=[len(a) for a in rest], **fields)
def _primary_amine_parent(info: dict) -> dict:
    aliph = _aliphatic_entries(info, "amines")
    prim = next((a for a in aliph if "c_idx" in a), None)
    if prim is None:
        prim = next((a for a in info.get("amines") or [] if "c_idx" in a), None)
    if prim is None: return _parent_dict(_longest_chain(info["mol"]), "alkane")
    # Gate: reject open-chain parent when amine carbon is in a ring
    # (mirrors _oh_c_in_ring gate in _chain_alcohol_parent, and _open_mono_fg_ok
    # gate for carboxyls/ketones/aldehydes/esters — prevents ring-NH₂ collapsing
    # to methanamine via a single-atom _chain_through result).
    if info["mol"].GetAtomWithIdx(prim["c_idx"]).IsInRing():
        return _parent_dict(_longest_chain(info["mol"]), "alkane")
    return _parent_dict(_chain_through(info, prim["c_idx"]), "amine", amine_c_idx=prim["c_idx"])
def _ring_amine_parent(info: dict) -> dict | None:
    top = (_try_btzam(info) or _try_boxam(info) or _try_bimam(info) or _try_bfam(info)
           or _try_pyrimidinamine_parent(info) or _try_pyridin_fg_parent(info))
    if top is not None: return top
    if _is_simple_aniline(info): return _aniline_parent(info)
    top = _try_arene_amine(info)
    if top is not None: return top
    return _ring_fg_try(info, ((_is_simple_cycloamine, "cycloamine"),), "amines", "amine_c_idx")
def _amine_parent(info: dict) -> dict:
    ring = _ring_amine_parent(info)
    if ring is not None: return ring
    for n in (4, 3, 2):
        p = _polyamine_parent(info, n)
        if p is not None:
            return p
    return _tert_amine_parent(info) or _sec_amine_parent(info) or _primary_amine_parent(info)
_ETHER_BAD = _CORE_BAD + ("has_amine", "has_alcohol", "has_thiol", "has_sulfide")
_SULFIDE_BAD = _CORE_BAD + ("has_amine", "has_alcohol", "has_thiol", "has_ether")
def _ether_parent(info: dict) -> dict | None:
    """Dialkyl ether parent (linear + isopropyl/HFIP specials)."""
    from namepredict.layer2.ether_parent import _ether_parent as _try
    return _try(info)
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
_ALKENOIC_BAD = _DIACID_BAD + ("has_thiol",)
_UNSAT_FG_BASE = (
    "has_acid", "has_ester", "has_amide", "has_ketone", "has_amine",
    "has_acyl_chloride", "has_anhydride", "has_thiol",
)
_UNSAT_FG_CORE = _UNSAT_FG_BASE + ("has_alcohol",)
_ALKENAL_BAD = _UNSAT_FG_CORE + ("has_nitrile",)
_ALKENENITRILE_BAD = _UNSAT_FG_CORE + ("has_aldehyde",)
_ALKENOL_BAD = _UNSAT_FG_BASE + ("has_aldehyde", "has_nitrile")
_ALKENOATE_BAD = tuple(k for k in _UNSAT_FG_BASE + ("has_nitrile",) if k != "has_ester")
_ALKENONE_BAD = tuple(k for k in _ALKENAL_BAD + ("has_aldehyde",) if k != "has_ketone")
def _with_anion(info: dict, parent: dict) -> dict:
    acids = info.get("carboxyls") or []
    return {**parent, "anion": True} if acids and all(c.get("anion") for c in acids) else parent
def _ring_acid_try(info: dict) -> dict | None:
    for fn in (_try_shcooh, _try_h5cooh, _try_indcooh, _try_naphcooh, _try_qcooh,
               _try_pyridinecarboxylic_parent, _try_benzene_polycarboxylic,
               _try_benzoic_parent,
               _try_cycloalkanecarboxylic_parent):
        if (b := fn(info)) is not None: return b
    return None
def _poly_acid_try(info: dict) -> dict | None:
    cm = carboxymethyl_diacid_parent(info)
    if cm is not None:
        return cm
    top = try_polycarboxylic_parent(info)
    if top is not None: return top
    if _is_simple_alkanedioic(info): return _diacid_parent(info)
    if _is_simple_alkenedioic(info): return _alkenedioic_parent(info)
    return _try_hy_alkenoic(info, _ALKENOIC_BAD, _best_cover_pair, _parent_dict, _db_pairs)
def _open_mono_fg_ok(info: dict, ekey: str) -> bool:
    """Open mono-FG only if FG carbon has open-chain C arm (or bare CX2)."""
    entries = info.get(ekey) or []
    if len(entries) != 1:
        return True
    mol, c = info["mol"], entries[0]["c_idx"]
    if _carbon_neighbors(mol, c):
        return True
    return not any(n.GetAtomicNum() == C for n in mol.GetAtomWithIdx(c).GetNeighbors())
def _acid_parent_core(info: dict) -> dict | None:
    top = _ring_acid_try(info) or _poly_acid_try(info)
    if top is not None:
        return top
    if not _open_mono_fg_ok(info, "carboxyls"):
        return None
    return _unsat_or_sat(
        info, "has_acid", "carboxyls", _ALKENOIC_BAD, "acid", "acid", "cooh_c_idx",
    )
def _acid_parent(info: dict) -> dict | None:
    core = _acid_parent_core(info)
    return None if core is None else _with_anion(info, core)
def _ring_or_poly_ketone(info: dict) -> dict | None:
    cyc = try_cycloalkanedione(info)
    if cyc is not None: return cyc
    if _is_simple_cycloalkenone(info):
        return _cyclo_ene_fg_parent(info, "cycloketone", "ketones", "ketone_c_idx")
    if _is_simple_cycloketone(info):
        return _cyclo_fg_parent(info, "cycloketone", "ketones", "ketone_c_idx")
    return _dione_parent(info) if _is_simple_alkanedione(info) else None
def _ketone_parent(info: dict) -> dict | None:
    a = _try_acetophenone_parent(info)
    if a is not None: return a
    top = _ring_or_poly_ketone(info)
    if top is not None: return top
    if not _open_mono_fg_ok(info, "ketones"): return None
    return _unsat_or_sat(
        info, "has_ketone", "ketones", _ALKENONE_BAD, "ketone", "ketone", "ketone_c_idx",
    )
def _aldehyde_parent(info: dict) -> dict | None:
    top = _try_izald(info) or _try_arene_fg_aldehyde(info) or _try_benzaldehyde_parent(info)
    if top is not None:
        return top
    if not _open_mono_fg_ok(info, "aldehydes"):
        return None
    return _unsat_or_sat(
        info, "has_aldehyde", "aldehydes", _ALKENAL_BAD, "aldehyde", "aldehyde", "aldehyde_c_idx",
    )
def _amide_parent(info: dict) -> dict:
    from namepredict.layer2.alkenamide import _amide_parent as _ap
    return _ap(info)
def _acyl_chloride_parent(info: dict) -> dict:
    e = info["acyl_chlorides"][0]
    hz = int(e.get("hal_z") or Cl)
    kind = "acyl_bromide" if hz == Br else "acyl_chloride"
    return _parent_dict(
        _chain_through(info, e["c_idx"]), kind, acyl_c_idx=e["c_idx"],
        cl_idx=e["cl_idx"], hal_idx=e.get("hal_idx", e["cl_idx"]), hal_z=hz,
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
    return (_try_izcn(info) or _try_pycn(info) or _try_arene_fg_nitrile(info)
            or _unsat_or_sat(
                info, "has_nitrile", "nitriles", _ALKENENITRILE_BAD,
                "nitrile", "nitrile", "nitrile_c_idx",
            ))
def _ester_meta(info: dict) -> dict:
    from namepredict.layer2.alkoxy_side import classify_alkoxy
    e = info["esters"][0]
    side = classify_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"])
    return dict(o_idx=e["o_idx"], alkoxy_c_idx=e["alkoxy_c_idx"], **side)
def _ester_parent(info: dict) -> dict | None:
    if not _open_mono_fg_ok(info, "esters"):
        return None
    return _unsat_or_sat(
        info, "has_ester", "esters", _ALKENOATE_BAD, "ester", "ester",
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
def _polyene_chain(info: dict) -> list[int]:
    mol, atoms = info["mol"], _db_atoms(info)
    chain = _longest_chain(mol)
    if _covers(chain, atoms):
        return chain
    return _best_cover_pair(mol, atoms) or chain
def _polyene_parent(info: dict) -> dict:
    chain = _polyene_chain(info)
    return _parent_dict(chain, "polyene", double_bonds=_db_pairs(info))
def _ring_atoms(info: dict) -> list[int]:
    return list(info["rings"][0]["atom_ids"])
def _cycloalkane_parent(info: dict) -> dict:
    from namepredict.layer2.scaffold.ring_parent import _pick_cycloalkane_ring as _pcr
    return _parent_dict(list(_pcr(info) or _ring_atoms(info)), "cycloalkane")
def _benzene_parent(info: dict) -> dict:
    from namepredict.layer2.scaffold.benzene_pick import _pick_benzene_ring as _pbr
    return _parent_dict(_pbr(info) or _ring_atoms(info), "benzene")
def _cycloalkene_parent(info: dict) -> dict:
    c = _ring_atoms(info); return _parent_dict(c, "cycloalkene", double_bond=_endocyclic_double(info, set(c)))
def _cyclo_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    return _parent_dict(_ring_atoms(info), kind, **{ckey: info[ekey][0]["c_idx"]})
def _cyclo_ene_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    c = _ring_atoms(info)
    db = _endocyclic_double(info, set(c))
    return _parent_dict(c, kind, double_bond=db, **{ckey: info[ekey][0]["c_idx"]})
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
def _unsat_parent(info: dict) -> dict | None:
    if _is_mono_alkyne(info): return _alkyne_parent(info)
    if _is_polyene(info): return _polyene_parent(info)
    return _alkene_parent(info) if _is_mono_alkene(info) else None
def _ring_parent(info: dict) -> dict | None:
    from namepredict.layer2.candidates import _ring_parent as _ring_best
    return _ring_best(info)
def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)


def _finalize_ranked(info: dict, cands: list[dict]) -> list[dict]:
    from namepredict.layer2.kind_registry import pack_parent_stem
    from namepredict.layer2.parent_candidate import with_principal_group_contract
    from namepredict.layer2.parent_ownership import finalize_parent_ownership
    mol = info.get("mol")
    return [
        finalize_parent_ownership(
            pack_parent_stem(with_principal_group_contract(c), mol), mol,
        )
        for c in _rank_candidates(info, cands)
    ]


def iter_parent_candidates(info: dict) -> list[dict]:
    """Ranked parent candidates, each finalized with immutable owned_atoms."""
    from namepredict.layer2.candidates import _alkane_fallback, _collect_candidates
    cands = _collect_candidates(info) or [_alkane_fallback(info)]
    return _finalize_ranked(info, cands)


def select_parent(info: dict) -> dict:
    """Legacy: first finalized ranked candidate (or alkane fallback)."""
    return iter_parent_candidates(info)[0]
