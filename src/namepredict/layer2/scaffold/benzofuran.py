"""Retained benzofuran / benzofuranamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]furan. O = 1; mono-methyl / mono-halo;
mono primary amine → benzofuranamine.

Engine via fused56; MonoSpec + ScaffoldSpec from scaffold.builders.fused56.
"""
from __future__ import annotations

from namepredict.layer2.scaffold.fused56 import (
    _fg_block,
    _mono_parent_dict,
    _mono_parts,
    _ring_set_parts,
    _subs_ok_cap,
    _try_mono_fused56,
)
from namepredict.layer2.scaffold.ring_parent import _mono_amine_on_ring, _outside_ok
from namepredict.layer2.scaffold.builders.fused56 import BF_AMINE, BF_MONO

_AMINE_BLOCK = (
    "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
    "has_ester", "has_amide", "has_nitrile", "has_thiol",
    "has_nitro", "has_acyl_chloride", "has_anhydride",
)


def _try_benzofuran_parent(info: dict) -> dict | None:
    return _try_mono_fused56(info, BF_MONO)


def _bf_primary_amine(info: dict, ring: set[int]) -> dict | None:
    am = _mono_amine_on_ring(info, ring)
    return am if am is not None and am.get("degree") == 1 else None


def _is_simple_benzofuranamine(info: dict) -> bool:
    if _mono_parts(info, BF_MONO) is None or _fg_block(info, _AMINE_BLOCK):
        return False
    parts = _mono_parts(info, BF_MONO)
    assert parts is not None
    mol, ring = info["mol"], _ring_set_parts(parts)
    am = _bf_primary_amine(info, ring)
    if am is None or not _outside_ok(mol, ring, {am["n_idx"]}):
        return False
    return _subs_ok_cap(mol, ring, BF_MONO.sub_cap)


def _benzofuranamine_parent(info: dict) -> dict:
    return _mono_parent_dict(
        info, BF_AMINE, amine_c_idx=info["amines"][0]["c_idx"],
    )


def _try_benzofuranamine_parent(info: dict) -> dict | None:
    return (
        _benzofuranamine_parent(info)
        if _is_simple_benzofuranamine(info) else None
    )
