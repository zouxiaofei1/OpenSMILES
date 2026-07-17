"""Saturated monohetero lactones / lactams (IUPAC P-65.6.3.5.1 / P-66.1.5.1).

PIN: heterocyclic pseudoketone — oxolan-2-one / pyrrolidin-2-one. Ring core is
a retained sat_hetero mono-O or mono-N; carbonyl C and hetero are both in-ring
and adjacent (cyclic ester O–C(=O) or amide N–C(=O)). Optional ≤1 simple ring
prefix: mono-methyl / mono-halo. N-unsub lactams only (first-cut).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _dbl_o_idx,
    _hetero_or_ring_halo,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)
from namepredict.layer2.sat_hetero import (
    _extra_amine,
    _hetero_asym,
    _hetero_pair,
    _kind_of,
    _n_unsub,
    _ring_atoms_if_mono,
    _ring_n_set,
    _sides_cover,
)

# base sat_hetero kind → lactone/lactam parent kind
_ONE_KIND = {
    "oxolane": "oxolanone",
    "oxane": "oxanone",
    "pyrrolidine": "pyrrolidinone",
    "piperidine": "piperidinone",
}


def _one_base_kind(info: dict) -> str | None:
    kind = _kind_of(info)
    return kind if kind in _ONE_KIND else None


def _one_fg_ok(info: dict) -> bool:
    """Exactly one cyclic ester or one cyclic amide; no competing principal FG."""
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return False
    if info.get("has_nitrile") or info.get("has_thiol") or info.get("has_nitro"):
        return False
    if info.get("has_alcohol"):
        return False
    n_e, n_a = len(info.get("esters") or []), len(info.get("amides") or [])
    return (n_e == 1 and n_a == 0) or (n_e == 0 and n_a == 1)


def _lactone_entry(info: dict, ring: set[int]) -> tuple[dict, int, int] | None:
    e = (info.get("esters") or [None])[0]
    if e is None:
        return None
    c, h = e["c_idx"], e["o_idx"]
    if c in ring and h in ring and e.get("alkoxy_c_idx") in ring:
        return e, c, h
    return None


def _lactam_entry(info: dict, ring: set[int]) -> tuple[dict, int, int] | None:
    am = (info.get("amides") or [None])[0]
    if am is None or am["c_idx"] not in ring or am["n_idx"] not in ring:
        return None
    if any(x not in ring for x in (am.get("n_c_idxs") or [])):
        return None
    return am, am["c_idx"], am["n_idx"]


def _entry_and_hetero(info: dict, ring: set[int]) -> tuple[dict, int, int] | None:
    """Return (entry, carbonyl_c, hetero_idx) when cyclic ester/amide on ring."""
    return _lactone_entry(info, ring) or _lactam_entry(info, ring)


def _one_ring_ctx(info: dict) -> tuple[Mol, list[int], set[int], str] | None:
    kind = _one_base_kind(info) if _one_fg_ok(info) else None
    atom_ids = _ring_atoms_if_mono(info) if kind else None
    if atom_ids is None:
        return None
    return info["mol"], atom_ids, set(atom_ids), kind


def _one_ctx(info: dict) -> tuple[Mol, list[int], set[int], int, int, str] | None:
    """(mol, atom_ids, ring, co_c, hetero, base_kind) for simple sat_hetero-one."""
    base = _one_ring_ctx(info)
    if base is None:
        return None
    mol, atom_ids, ring, kind = base
    got = _entry_and_hetero(info, ring)
    if got is None:
        return None
    _, co_c, hetero = got
    return mol, atom_ids, ring, co_c, hetero, kind


def _one_amine_ok(info: dict, mol: Mol, atom_ids: list[int]) -> bool:
    ring_ns = _ring_n_set(mol, atom_ids)
    return not _extra_amine(info, ring_ns) and _n_unsub(mol, ring_ns)


def _one_prefix_n(mol: Mol, ring: set[int]) -> int | None:
    """Count mono-halo + claimable ring-C alkyl sides (C1–C4 / branched)."""
    starts = _ring_side_starts(mol, ring)
    if not _sides_cover(mol, ring, starts):
        return None
    return _ring_halo_n(mol, ring) + len(starts)


def _one_subs_ok(info: dict, mol: Mol, ring: set[int], co_c: int) -> bool:
    """Allow ≤1 ring-C alkyl or mono-halo; carbonyl O outside is allowed."""
    o_idx = _dbl_o_idx(mol, co_c)
    n = _one_prefix_n(mol, ring)
    if o_idx is None or n is None or n > 1:
        return False
    allow = {o_idx}
    return _hetero_or_ring_halo(mol, ring, allow) and _outside_ok(mol, ring, allow)


def _is_simple_sat_hetero_one(info: dict) -> bool:
    ctx = _one_ctx(info)
    if ctx is None:
        return False
    mol, atom_ids, ring, co_c, _, _ = ctx
    if not _one_amine_ok(info, mol, atom_ids):
        return False
    return _one_subs_ok(info, mol, ring, co_c)


def _sat_hetero_one_parent(info: dict) -> dict:
    mol, atom_ids, _, co_c, hetero, base = _one_ctx(info)
    hs = _hetero_pair(info)
    return {
        "chain": atom_ids, "n_carbons": len(atom_ids), "kind": _ONE_KIND[base],
        "base_kind": base, "hetero_idxs": hs, "hetero_asym": _hetero_asym(info, hs),
        "hetero_idx": hetero, "ketone_c_idx": co_c, "ring_attach_idx": co_c,
        "one_o_idx": _dbl_o_idx(mol, co_c),
    }


def _try_sat_hetero_one_parent(info: dict) -> dict | None:
    if not _is_simple_sat_hetero_one(info):
        return None
    return _sat_hetero_one_parent(info)
