"""Bootstrap unsat parent producers into kind_registry (single registration site).

candidates._unsat_candidates iterates kind_registry.unsat_try_fns(); add new
open-chain unsaturation hydrocarbon parents here instead of editing a hand
tuple in candidates. Keep try order stable: alkyne → polyene → alkene.
"""
from __future__ import annotations

from namepredict.layer2.kind_registry import register_unsat_try
from namepredict.layer2.parent_selector import (
    _alkene_parent,
    _alkyne_parent,
    _is_mono_alkene,
    _is_mono_alkyne,
    _is_polyene,
    _polyene_parent,
)


def _try_alkyne(info: dict) -> dict | None:
    return _alkyne_parent(info) if _is_mono_alkyne(info) else None


def _try_polyene(info: dict) -> dict | None:
    return _polyene_parent(info) if _is_polyene(info) else None


def _try_alkene(info: dict) -> dict | None:
    return _alkene_parent(info) if _is_mono_alkene(info) else None


_UNSAT_PRODUCERS = (
    _try_alkyne,
    _try_polyene,
    _try_alkene,
)


def _bootstrap() -> None:
    for fn in _UNSAT_PRODUCERS:
        register_unsat_try(fn)


_bootstrap()
