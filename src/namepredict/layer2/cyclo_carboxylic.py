"""Monocyclic cycloalkane + one exocyclic carbonyl FG (P-65 / P-66).

Covers cycloalkanecarboxylic acid and horizontal peers: carbaldehyde,
carbonitrile, carboxamide, carboxylate (alkyl ester), carbonyl halide.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _benzoate_alkoxy,
    _cooh_oxygen_idxs,
    _ester_exclude,
    _fg_ring_c,
    _nitrile_n_idx,
)
from namepredict.layer2.scaffold.ring_parent import (
    _dbl_o_idx,
    _hetero_or_ring_halo,
    _is_cycloalkane_core,
    _outside_carbons,
    _ring_halo_n,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _disjoint_cover, _side_sets

# flag, ekey, kind, parent_c_key
_FG_SPECS: tuple[tuple[str, str, str, str], ...] = (
    ("has_acid", "carboxyls", "cycloalkanecarboxylic", "cooh_c_idx"),
    ("has_aldehyde", "aldehydes", "cycloalkanecarbaldehyde", "aldehyde_c_idx"),
    ("has_nitrile", "nitriles", "cycloalkanecarbonitrile", "nitrile_c_idx"),
    ("has_amide", "amides", "cycloalkanecarboxamide", "amide_c_idx"),
    ("has_ester", "esters", "cycloalkanecarboxylate", "ester_c_idx"),
    ("has_acyl_chloride", "acyl_chlorides", "cycloalkanecarbonyl_halide", "acyl_c_idx"),
)


def _subs_ok(mol: Mol, ring: set[int], exclude: set[int]) -> bool:
    """Allow recognized side-chain topologies (linear, branched, CF3, etc.) via probes."""
    starts = _ring_side_starts(mol, ring, exclude)
    outside = set(_outside_carbons(mol, ring, exclude))
    if not starts:
        return not outside
    sets = _side_sets(mol, ring, starts)
    return sets is not None and _disjoint_cover(sets, outside)


def _fg_entry(info: dict, ekey: str) -> dict | None:
    entries = info.get(ekey) or []
    return entries[0] if len(entries) == 1 else None


def _with_dbl_o(mol: Mol, c_idx: int, base: set[int]) -> set[int]:
    o = _dbl_o_idx(mol, c_idx)
    return base | ({o} if o is not None else set())


def _simple_pack(mol: Mol, e: dict, extra: set[int]) -> tuple[set[int], set[int]]:
    return {e["c_idx"]}, _with_dbl_o(mol, e["c_idx"], extra)


def _nitrile_pack(mol: Mol, e: dict) -> tuple[set[int], set[int]] | None:
    n = _nitrile_n_idx(mol, e["c_idx"])
    return ({e["c_idx"]}, {n}) if n is not None else None


def _amide_pack(mol: Mol, e: dict) -> tuple[set[int], set[int]] | None:
    if e.get("n_c_idxs"):
        return None  # primary amide only this round
    return _simple_pack(mol, e, {e["n_idx"]})


def _acyl_pack(mol: Mol, e: dict) -> tuple[set[int], set[int]]:
    h = e.get("hal_idx", e.get("cl_idx"))
    return _simple_pack(mol, e, {h} if h is not None else set())


def _acid_pack(mol: Mol, e: dict) -> tuple[set[int], set[int]]:
    return {e["c_idx"]}, _cooh_oxygen_idxs(mol, e["c_idx"])


def _aldehyde_pack(mol: Mol, e: dict) -> tuple[set[int], set[int]]:
    return _simple_pack(mol, e, set())


_PACKERS = {
    "carboxyls": _acid_pack, "aldehydes": _aldehyde_pack, "nitriles": _nitrile_pack,
    "amides": _amide_pack, "esters": lambda m, e: _ester_exclude(m, e),
    "acyl_chlorides": _acyl_pack,
}


def _fg_pack(info: dict, ekey: str) -> tuple[set[int], set[int]] | None:
    """(exclude carbons, allowed hetero idxs) for one exocyclic FG."""
    e, fn = _fg_entry(info, ekey), _PACKERS.get(ekey)
    return None if e is None or fn is None else fn(info["mol"], e)


def _pick_mono_sat_ring(info: dict) -> set[int] | None:
    from namepredict.layer2.scaffold.cyclo_pick import _cyclo_parent_candidates
    cands = _cyclo_parent_candidates(info)
    return cands[0] if len(cands) == 1 else None


def _ester_side_ok(info: dict) -> bool:
    e = info["esters"][0]
    return _benzoate_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"]) is not None


def _is_simple_cyclo_exo_fg(info: dict, flag: str, ekey: str) -> bool:
    if not info.get(flag) or not _is_cycloalkane_core(info):
        return False
    ring, pack = _pick_mono_sat_ring(info), _fg_pack(info, ekey)
    if ring is None or pack is None or _fg_ring_c(info, ring, ekey) is None:
        return False
    if ekey == "esters" and not _ester_side_ok(info):
        return False
    excl, allowed = pack
    mol = info["mol"]
    return _hetero_or_ring_halo(mol, ring, allowed) and _subs_ok(mol, ring, excl)


def _ester_fields(info: dict, e: dict) -> dict:
    side = _benzoate_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"]) or {}
    return {
        "o_idx": e["o_idx"], "alkoxy_c_idx": e["alkoxy_c_idx"],
        "alkoxy_n": side.get("alkoxy_n"), "alkoxy_en": side.get("alkoxy_en") or "",
        "alkoxy_zh": side.get("alkoxy_zh") or "",
    }


def _acyl_fields(e: dict) -> dict:
    hz = int(e.get("hal_z") or 17)
    kind = "cycloalkanecarbonyl_bromide" if hz == 35 else "cycloalkanecarbonyl_chloride"
    return {
        "cl_idx": e.get("cl_idx"), "hal_idx": e.get("hal_idx", e.get("cl_idx")),
        "hal_z": hz, "kind": kind,
    }


def _extra_fields(info: dict, ekey: str, e: dict) -> dict:
    if ekey == "esters":
        return _ester_fields(info, e)
    return _acyl_fields(e) if ekey == "acyl_chlorides" else {}


def _cyclo_exo_parent(info: dict, ekey: str, kind: str, ckey: str) -> dict:
    ring = list(_pick_mono_sat_ring(info) or ())
    e = info[ekey][0]
    base = {
        "chain": ring, "n_carbons": len(ring), "kind": kind,
        ckey: e["c_idx"], "ring_attach_idx": _fg_ring_c(info, set(ring), ekey),
    }
    return {**base, **_extra_fields(info, ekey, e)}


def _try_cyclo_exo_fg(info: dict, flag: str, ekey: str, kind: str, ckey: str) -> dict | None:
    if not _is_simple_cyclo_exo_fg(info, flag, ekey):
        return None
    return _cyclo_exo_parent(info, ekey, kind, ckey)


def _try_specs(info: dict, specs) -> dict | None:
    for flag, ekey, kind, ckey in specs:
        if (p := _try_cyclo_exo_fg(info, flag, ekey, kind, ckey)) is not None:
            return p
    return None


def _is_simple_cycloalkanecarboxylic(info: dict) -> bool:
    return _is_simple_cyclo_exo_fg(info, "has_acid", "carboxyls")


def _cycloalkanecarboxylic_parent(info: dict) -> dict:
    return _cyclo_exo_parent(info, "carboxyls", "cycloalkanecarboxylic", "cooh_c_idx")


def _try_cycloalkanecarboxylic_parent(info: dict) -> dict | None:
    return _try_cyclo_exo_fg(
        info, "has_acid", "carboxyls", "cycloalkanecarboxylic", "cooh_c_idx",
    )


def _try_cycloalkane_exocyclic_fg(info: dict) -> dict | None:
    """Non-acid exocyclic FG (acid already via _ring_acid_try)."""
    return _try_specs(info, _FG_SPECS[1:])
