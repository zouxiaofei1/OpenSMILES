"""Monocarbocycle mono-ene + mono FG (cycloalkenol / cycloalkenone; P-31.1)."""
from __future__ import annotations

from namepredict.layer2.ring_parent import (
    _cyclo_ene_fg_ok,
    _dbl_o_idx,
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
