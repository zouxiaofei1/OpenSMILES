"""Scaffold specs registry (L2 data only; no naming assembly)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NumberingPolicy:
    mode: str
    standard_path: tuple = ()
    anchors: tuple[str, ...] = ()
    substitutable: frozenset[str] = frozenset()


@dataclass(frozen=True)
class ScaffoldSpec:
    id: str
    naming_class: str
    stem_en: str | None
    stem_zh: str | None
    n_rings: int
    ring: str
    retained: bool
    fg_rank: int
    numbering: NumberingPolicy
    sub_rules: object | None = None
    principal_slots: object | None = None


def _carbo(sid: str, mode: str) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="carbocycle", stem_en=None, stem_zh=None,
        n_rings=1, ring="carbo", retained=False, fg_rank=0, numbering=pol,
    )


CARBOCYCLE_SPECS: tuple[ScaffoldSpec, ...] = (
    _carbo("cycloalkane", "carbocycle_free"),
    _carbo("cycloalkene", "carbocycle_free"),
    _carbo("cyclopolyene", "poly_unsat"),
)

_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in CARBOCYCLE_SPECS}


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    return CARBOCYCLE_SPECS
