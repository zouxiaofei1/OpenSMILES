from __future__ import annotations
from namepredict.layer2.arene_carbonyl import _try_acetophenone_parent
from namepredict.layer2.cyclo_poly_fg import try_cycloalkanedione
from namepredict.layer2.cyclo_ene_fg import _is_simple_cycloalkenone
from namepredict.layer2.scaffold.ring_parent import _endocyclic_double, _is_simple_cycloketone
from namepredict.layer2.aliph_fg import _aliph_c_idxs, _c_idxs
from namepredict.layer2.chain_walk import _carbon_neighbors, _longest_chain
from namepredict.constants import C
from namepredict.layer2.parent_selector_common import _DIONE_BAD, _no_fgs
from namepredict.layer2.parent_core import _best_cover_pair, _is_open_sat, _parent_dict, _unsat_or_sat
def _is_simple_n(info: dict, bad: tuple, ekey: str, n: int) -> bool:
    return (
        _is_open_sat(info) and _no_fgs(info, bad)
        and _aliph_c_idxs(info, ekey, n) is not None
    )
def _cover_parent(info: dict, ekey: str, n: int, kind: str, key: str) -> dict:
    atoms = _aliph_c_idxs(info, ekey, n) or _c_idxs(info.get(ekey) or [], n) or []
    chain = _best_cover_pair(info["mol"], atoms) or _longest_chain(info["mol"])
    return _parent_dict(chain, kind, **{key: atoms})
def _is_simple_alkanedione(info: dict) -> bool:
    return _is_simple_n(info, _DIONE_BAD, "ketones", 2)
def _dione_parent(info: dict) -> dict:
    return _cover_parent(info, "ketones", 2, "dione", "ketone_c_idxs")

_UNSAT_FG_BASE = (
    "has_acid", "has_ester", "has_amide", "has_ketone", "has_amine",
    "has_acyl_chloride", "has_anhydride", "has_thiol",
)
_UNSAT_FG_CORE = _UNSAT_FG_BASE + ("has_alcohol",)
_ALKENAL_BAD = _UNSAT_FG_CORE + ("has_nitrile",)
_ALKENONE_BAD = tuple(k for k in _ALKENAL_BAD + ("has_aldehyde",) if k != "has_ketone")
def _open_mono_fg_ok(info: dict, ekey: str) -> bool:
    """Open mono-FG only if FG carbon has open-chain C arm (or bare CX2)."""
    entries = info.get(ekey) or []
    if len(entries) != 1:
        return True
    mol, c = info["mol"], entries[0]["c_idx"]
    if _carbon_neighbors(mol, c):
        return True
    return not any(n.GetAtomicNum() == C for n in mol.GetAtomWithIdx(c).GetNeighbors())
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
def _ring_atoms(info: dict) -> list[int]:
    return list(info["rings"][0]["atom_ids"])
def _cyclo_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    return _parent_dict(_ring_atoms(info), kind, **{ckey: info[ekey][0]["c_idx"]})
def _cyclo_ene_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    c = _ring_atoms(info)
    db = _endocyclic_double(info, set(c))
    return _parent_dict(c, kind, double_bond=db, **{ckey: info[ekey][0]["c_idx"]})
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
