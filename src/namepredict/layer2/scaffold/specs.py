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


# Shared fused 5+6 path labels (IUPAC P-22.2.1 / P-25): hetero=1 … 7a.
FUSED56_LABELS: tuple[str, ...] = (
    "1", "2", "3", "3a", "4", "5", "6", "7", "7a",
)


def _carbo(sid: str, mode: str) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="carbocycle", stem_en=None, stem_zh=None,
        n_rings=1, ring="carbo", retained=False, fg_rank=0, numbering=pol,
    )


def _fused56(
    sid: str, stem_en: str, stem_zh: str, fg_rank: int = 0,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fused56_fixed", standard_path=FUSED56_LABELS)
    return ScaffoldSpec(
        id=sid, naming_class="fused56", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring="hetero", retained=True, fg_rank=fg_rank, numbering=pol,
    )


CARBOCYCLE_SPECS: tuple[ScaffoldSpec, ...] = (
    _carbo("cycloalkane", "carbocycle_free"),
    _carbo("cycloalkene", "carbocycle_free"),
    _carbo("cyclopolyene", "poly_unsat"),
)

FUSED56_SPECS: tuple[ScaffoldSpec, ...] = (
    _fused56("benzofuran", "benzofuran", "苯并呋喃"),
    _fused56("benzofuranamine", "benzofuranamine", "苯并呋喃胺", fg_rank=3),
    _fused56("benzothiophene", "1-benzothiophene", "苯并[b]噻吩"),
    _fused56("benzothiophenol", "1-benzothiophenol", "苯并[b]噻吩酚", fg_rank=5),
    _fused56("benzothiazole", "1,3-benzothiazole", "1,3-苯并噻唑"),
    _fused56("benzothiazolamine", "benzothiazolamine", "苯并噻唑胺", fg_rank=3),
    _fused56("benzoxazole", "1,3-benzoxazole", "1,3-苯并噁唑"),
    _fused56("benzoxazolamine", "benzoxazolamine", "苯并噁唑胺", fg_rank=3),
)

_ALL_SPECS: tuple[ScaffoldSpec, ...] = CARBOCYCLE_SPECS + FUSED56_SPECS
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    return _ALL_SPECS


def fused56_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids that use fused56 1…7a labels."""
    return frozenset(s.id for s in FUSED56_SPECS)
