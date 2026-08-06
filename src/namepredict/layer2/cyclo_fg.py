"""Monocarbocycle mono/di ene + FG parents (P-31.1 / P-63.1.2 / P-64.2.1).

Cycloalkenol / cycloalkenone (mono ene + mono FG) and cycloalkanediol /
cycloalkanedione (two ring-carbon OH or two ring ketone carbonyls on an
unfused C3–C10 sat carbocycle).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.parent_core import _parent_dict
from namepredict.layer2.scaffold.ring_parent import (
    _cyclo_ene_fg_ok,
    _cyclo_fg_parent_ok,
    _dbl_o_idx,
    _di_oh_on_ring,
    _mono_ketone_on_ring,
    _mono_oh_on_ring,
)


def _ring_set0(info: dict) -> set[int]:
    return set((info.get("rings") or [{}])[0].get("atom_ids") or [])


def _is_simple_cycloalkenol(info: dict) -> bool:
    ring = _ring_set0(info)
    oh = _mono_oh_on_ring(info, ring)
    return bool(oh) and _cyclo_ene_fg_ok(info, {oh["o_idx"]})


def _is_simple_cycloalkenone(info: dict) -> bool:
    ring = _ring_set0(info)
    ket = _mono_ketone_on_ring(info, ring)
    if ket is None:
        return False
    o_idx = _dbl_o_idx(info["mol"], ket["c_idx"])
    return o_idx is not None and _cyclo_ene_fg_ok(info, {o_idx})


def _di_ketone_on_ring(info: dict, ring_set: set[int]) -> list[dict] | None:
    ketones = info.get("ketones") or []
    if len(ketones) != 2 or any(k["c_idx"] not in ring_set for k in ketones):
        return None
    return ketones


def _ketone_o_set(mol: Mol, kets: list[dict]) -> set[int] | None:
    out: set[int] = set()
    for k in kets:
        o = _dbl_o_idx(mol, k["c_idx"])
        if o is None:
            return None
        out.add(o)
    return out


def _is_simple_cycloalkanediol(info: dict) -> bool:
    ring = _ring_set0(info)
    ohs = _di_oh_on_ring(info, ring)
    return bool(ohs) and _cyclo_fg_parent_ok(info, {h["o_idx"] for h in ohs})


def _is_simple_cycloalkanedione(info: dict) -> bool:
    ring = _ring_set0(info)
    kets = _di_ketone_on_ring(info, ring)
    if not kets:
        return False
    allowed = _ketone_o_set(info["mol"], kets)
    return bool(allowed) and _cyclo_fg_parent_ok(info, allowed)


def _cyclo_poly_fg_parent(info: dict, kind: str, ekey: str, ckey: str) -> dict:
    cs = [e["c_idx"] for e in info.get(ekey) or []]
    return _parent_dict(list(info["rings"][0]["atom_ids"]), kind, **{ckey: cs})


def try_cycloalkanediol(info: dict) -> dict | None:
    if not _is_simple_cycloalkanediol(info):
        return None
    return _cyclo_poly_fg_parent(info, "cycloalkanediol", "hydroxyls", "oh_c_idxs")


def try_cycloalkanedione(info: dict) -> dict | None:
    if not _is_simple_cycloalkanedione(info):
        return None
    return _cyclo_poly_fg_parent(info, "cycloalkanedione", "ketones", "ketone_c_idxs")
