"""Retained benzamide parent (Ph–C(=O)–N; IUPAC P-66.1.1 / P-65.1.1.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import (
    _arene_fg_ctx,
    _arene_subs_ok,
    _fg_ring_c,
    _pick_fg_ring,
)
from namepredict.layer2.ring_parent import _dbl_o_idx


def _amide_n_carbons(mol: Mol, am: dict) -> set[int]:
    """N-alkyl chain carbons or N-phenyl ring carbons (exclude from ring sides)."""
    from namepredict.layer2.aryl_sub import _phenyl_at
    from namepredict.layer2.chain_walk import _longest_from
    n_idx, out = am["n_idx"], set()
    for c in am.get("n_c_idxs") or []:
        ph = _phenyl_at(mol, c, n_idx)
        out |= ph if ph is not None else set(_longest_from(mol, c, set()) or [c])
    return out


def _amide_pack_benz(mol: Mol, am: dict) -> tuple[set[int], set[int]]:
    """(exclude carbons, allowed heteros) for Ph–C(=O)–N."""
    o_dbl = _dbl_o_idx(mol, am["c_idx"])
    allowed = {am["n_idx"]} | ({o_dbl} if o_dbl is not None else set())
    return {am["c_idx"]} | _amide_n_carbons(mol, am), allowed


def _benzamide_n_ok(info: dict) -> bool:
    """Primary or L3-simple N (C1–C4 alkyl / N,N / N-phenyl); reject complex N."""
    from namepredict.layer2.alkenamide import _amide_n_meta
    ams = info.get("amides") or []
    if len(ams) != 1:
        return False
    cs = ams[0].get("n_c_idxs") or []
    if not cs:
        return True
    meta = _amide_n_meta(info)
    return bool(meta) and not meta.get("n_benzyl")


def _is_simple_benzamide(info: dict) -> bool:
    ctx = _arene_fg_ctx(info, "has_amide", "amides")
    if ctx is None or not _benzamide_n_ok(info):
        return False
    mol, ring = ctx
    excl, allowed = _amide_pack_benz(mol, info["amides"][0])
    return _arene_subs_ok(info, mol, ring, excl, allowed)


def _benzamide_parent(info: dict) -> dict:
    from namepredict.layer2.alkenamide import _amide_n_meta
    ring = _pick_fg_ring(info, "amides") or set()
    am = info["amides"][0]
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "benzamide",
        "amide_c_idx": am["c_idx"], "ring_attach_idx": _fg_ring_c(info, ring, "amides"),
        **_amide_n_meta(info),
    }


def _try_benzamide_parent(info: dict) -> dict | None:
    return _benzamide_parent(info) if _is_simple_benzamide(info) else None
