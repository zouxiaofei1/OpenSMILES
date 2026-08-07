"""1,4-Benzoquinone parent: cyclohexa-2,5-diene-1,4-dione (IUPAC P-64.2).

Single C6 carbocycle + exactly two para ring ketones + two endocyclic C=C.
PIN is systematic (retained 1,4-benzoquinone is general-only).
"""
from __future__ import annotations

from namepredict.layer2.scaffold.ring_parent import (
    _is_carbocycle_ring,
)


def _ketone_idxs(info: dict) -> list[int]:
    return [e["c_idx"] for e in (info.get("ketones") or []) if "c_idx" in e]


def _bq_ring(info: dict) -> list[int] | None:
    atom_ids = _is_carbocycle_ring(info)
    return list(atom_ids) if atom_ids is not None and len(atom_ids) == 6 else None


def _ring_dist(ring: list[int], a: int, b: int) -> int:
    i, j = ring.index(a), ring.index(b)
    d = abs(i - j)
    return min(d, len(ring) - d)


def _ketones_para(ring: list[int], ket: list[int]) -> bool:
    if len(ket) != 2 or any(k not in ring for k in ket):
        return False
    return _ring_dist(ring, ket[0], ket[1]) == 3


_BQ_BLOCK = (
    "has_acid", "has_aldehyde", "has_ester", "has_amide",
    "has_nitrile", "has_amine", "has_thiol", "has_nitro",
    "has_acyl_chloride", "has_anhydride",
)


# ── 1,2-benzoquinone (ortho) — same substitution gates, ring distance only ──

def _ketones_ortho(ring: list[int], ket: list[int]) -> bool:
    if len(ket) != 2 or any(k not in ring for k in ket):
        return False
    return _ring_dist(ring, ket[0], ket[1]) == 1


