"""Five-membered heteroarene carboxylic acids (IUPAC P-65.1.1 / P-22.2.1).

Monocyclic furan / thiophene / 1H-pyrrole / 1H-imidazole / 1H-pyrazole with
exactly one ring-attached COOH. Optional ≤2 simple prefixes: halo / methyl /
primary amino. Chinese retained suffix 甲酸.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_conflict,
    _carboxyl_ring_c,
    _cooh_oxygen_idxs,
)
from namepredict.layer2.heteroarene5 import (
    _hetero5_idx,
    _hetero5_kind,
    _imidazole_n_pair,
    _is_hetero5_core,
    _mono_methyl_only,
    _pyrazole_n_pair,
)
from namepredict.layer2.ring_parent import (
    _hetero_or_ring_halo,
    _ring_halo_n,
    _ring_primary_amines,
    _ring_side_starts,
)

# base hetero kind → carboxylic parent kind
_CARBOXY_KIND = {
    "furan": "furancarboxylic",
    "thiophene": "thiophenecarboxylic",
    "pyrrole": "pyrrolecarboxylic",
    "imidazole": "imidazolecarboxylic",
    "pyrazole": "pyrazolecarboxylic",
}


def _h5cooh_prefix_n(info: dict, mol: Mol, ring: set[int], fg_c: int) -> int | None:
    """Count halo + methyl + primary amino prefixes; None if invalid side chains."""
    starts = _ring_side_starts(mol, ring, {fg_c})
    if not _mono_methyl_only(mol, ring | {fg_c}, starts):
        return None
    ams = _ring_primary_amines(info, ring)
    if len(info.get("amines") or []) != len(ams):
        return None
    return _ring_halo_n(mol, ring) + len(starts) + len(ams)


def _h5cooh_subs_ok(info: dict, mol: Mol, ring: set[int], fg_c: int) -> bool:
    """Allow ≤2 simple ring prefixes; COOH oxygens are principal FG atoms."""
    n = _h5cooh_prefix_n(info, mol, ring, fg_c)
    if n is None or n > 2:
        return False
    allowed = _cooh_oxygen_idxs(mol, fg_c) | {a["n_idx"] for a in _ring_primary_amines(info, ring)}
    return _hetero_or_ring_halo(mol, ring, allowed)


def _h5cooh_fg_ok(info: dict) -> bool:
    """Block non-acid principal FGs (amine is prefix, not principal)."""
    return not _arene_fg_conflict(info, "has_aldehyde", "has_ketone", "has_alcohol")


def _h5cooh_ctx(info: dict) -> tuple[Mol, set[int], int] | None:
    """Return (mol, ring_set, fg_c) when one ring COOH is present and FG ok."""
    rings = info.get("rings") or []
    if not _h5cooh_fg_ok(info) or len(rings) != 1 or not info.get("carboxyls"):
        return None
    mol, ring = info["mol"], set(rings[0]["atom_ids"])
    if _carboxyl_ring_c(info, ring) is None:
        return None
    return mol, ring, info["carboxyls"][0]["c_idx"]


def _is_simple_h5cooh(info: dict, core_ok: bool) -> bool:
    if not core_ok:
        return False
    ctx = _h5cooh_ctx(info)
    if ctx is None:
        return False
    mol, ring, fg_c = ctx
    return _h5cooh_subs_ok(info, mol, ring, fg_c)


def _h5cooh_base(info: dict, kind: str, **extra) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": ring, "n_carbons": 5, "kind": kind,
        "cooh_c_idx": cooh_c,
        "ring_attach_idx": _carboxyl_ring_c(info, set(ring)),
        **extra,
    }


def _is_simple_mono_h5cooh(info: dict) -> bool:
    return _is_simple_h5cooh(info, _is_hetero5_core(info) and _hetero5_kind(info) is not None)


def _mono_h5cooh_parent(info: dict) -> dict:
    base = _hetero5_kind(info)
    return _h5cooh_base(
        info, _CARBOXY_KIND[base], hetero_idx=_hetero5_idx(info), base_kind=base,
    )


def _try_mono_h5cooh_parent(info: dict) -> dict | None:
    return _mono_h5cooh_parent(info) if _is_simple_mono_h5cooh(info) else None


def _is_simple_imidazolecarboxylic(info: dict) -> bool:
    return _is_simple_h5cooh(info, _imidazole_n_pair(info) is not None)


def _imidazolecarboxylic_parent(info: dict) -> dict:
    nh, n = _imidazole_n_pair(info)
    return _h5cooh_base(
        info, "imidazolecarboxylic",
        nh_idx=nh, n_idx=n, n_idxs=[nh, n], base_kind="imidazole",
    )


def _try_imidazolecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_imidazolecarboxylic(info):
        return None
    return _imidazolecarboxylic_parent(info)


def _is_simple_pyrazolecarboxylic(info: dict) -> bool:
    return _is_simple_h5cooh(info, _pyrazole_n_pair(info) is not None)


def _pyrazolecarboxylic_parent(info: dict) -> dict:
    nh, n = _pyrazole_n_pair(info)
    return _h5cooh_base(
        info, "pyrazolecarboxylic",
        nh_idx=nh, n_idx=n, n_idxs=[nh, n], base_kind="pyrazole",
    )


def _try_pyrazolecarboxylic_parent(info: dict) -> dict | None:
    if not _is_simple_pyrazolecarboxylic(info):
        return None
    return _pyrazolecarboxylic_parent(info)


def _try_hetero5carboxylic_parent(info: dict) -> dict | None:
    """Prefer diazole carboxylic, else mono-hetero5 carboxylic."""
    for fn in (
        _try_imidazolecarboxylic_parent,
        _try_pyrazolecarboxylic_parent,
        _try_mono_h5cooh_parent,
    ):
        if (p := fn(info)) is not None:
            return p
    return None
