"""1,2-Benzoquinone parent: cyclohexa-3,5-diene-1,2-dione (IUPAC P-64.2).

Single C6 carbocycle + exactly two ortho ring ketones + two endocyclic C=C.
PIN is systematic (retained o-benzoquinone is general-only). Reuses 1,4-BQ
substitution gates; ortho vs para is ring distance only.
"""
from __future__ import annotations

from namepredict.layer2.scaffold.benzoquinone import (
    _bq_doubles_ok,
    _bq_fg_block,
    _bq_ring,
    _bq_subs_ok,
    _ketone_idxs,
    _ring_dist,
)
from namepredict.layer2.scaffold.ring_parent import _endocyclic_doubles


def _ketones_ortho(ring: list[int], ket: list[int]) -> bool:
    if len(ket) != 2 or any(k not in ring for k in ket):
        return False
    return _ring_dist(ring, ket[0], ket[1]) == 1


def _is_simple_ortho_benzoquinone(info: dict) -> bool:
    ring, ket = _bq_ring(info), _ketone_idxs(info)
    if ring is None or _bq_fg_block(info) or not _ketones_ortho(ring, ket):
        return False
    if not _bq_doubles_ok(info, ring):
        return False
    return _bq_subs_ok(info, info["mol"], ring, ket)


def _ortho_benzoquinone_parent(info: dict) -> dict:
    ring = _bq_ring(info) or []
    ket = _ketone_idxs(info)
    bonds = _endocyclic_doubles(info, set(ring))
    return {
        "chain": ring, "n_carbons": 6, "kind": "ortho_benzoquinone",
        "scaffold_id": "ortho_benzoquinone", "ketone_c_idxs": ket,
        "double_bonds": bonds,
    }


def _try_ortho_benzoquinone_parent(info: dict) -> dict | None:
    if not _is_simple_ortho_benzoquinone(info):
        return None
    return _ortho_benzoquinone_parent(info)
