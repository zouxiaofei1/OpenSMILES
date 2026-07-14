"""Retained 1-benzothiophene / benzothiophenol parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[b]thiophene. S = 1; mono-methyl / mono-halo;
mono ring OH → 1-benzothiophen-n-ol / 苯并[b]噻吩-n-醇.

Core gate/subs/parent dict via fused56.Fused56MonoSpec engine.
"""
from __future__ import annotations

from namepredict.layer2.fused56 import (
    Fused56MonoSpec,
    _fg_block,
    _mono_parent_dict,
    _mono_parts,
    _ring_set_parts,
    _subs_ok_cap,
    _try_mono_fused56,
)
from namepredict.layer2.ring_parent import _mono_oh_on_ring, _outside_ok

_BT = Fused56MonoSpec(
    kind="benzothiophene", hetero_z=16, hetero_key="s_idx", sub_cap=1,
)

_OL_BLOCK = (
    "has_acid", "has_aldehyde", "has_ketone", "has_ester",
    "has_amide", "has_nitrile", "has_amine", "has_thiol",
    "has_nitro", "has_acyl_chloride", "has_anhydride",
)


def _try_benzothiophene_parent(info: dict) -> dict | None:
    return _try_mono_fused56(info, _BT)


def _is_simple_benzothiophenol(info: dict) -> bool:
    if _mono_parts(info, _BT) is None or _fg_block(info, _OL_BLOCK):
        return False
    parts = _mono_parts(info, _BT)
    assert parts is not None
    mol, ring = info["mol"], _ring_set_parts(parts)
    oh = _mono_oh_on_ring(info, ring)
    if oh is None or not _outside_ok(mol, ring, {oh["o_idx"]}):
        return False
    return _subs_ok_cap(mol, ring, _BT.sub_cap)


def _benzothiophenol_parent(info: dict) -> dict:
    return _mono_parent_dict(
        info, Fused56MonoSpec(
            kind="benzothiophenol", hetero_z=16, hetero_key="s_idx", sub_cap=1,
        ),
        oh_c_idx=info["hydroxyls"][0]["c_idx"],
    )


def _try_benzothiophenol_parent(info: dict) -> dict | None:
    return (
        _benzothiophenol_parent(info)
        if _is_simple_benzothiophenol(info) else None
    )
