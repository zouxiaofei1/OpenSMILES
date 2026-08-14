"""Scaffold specs registry (L2 data only; no naming assembly)."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2.identity import ScaffoldIdentity, identity_of


@dataclass(frozen=True)
class NumberingPolicy:
    mode: str
    standard_path: tuple = ()
    materialize_plan: bool = True
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

    @property
    def identity(self) -> ScaffoldIdentity:
        return identity_of(self)


# Shared fused 5+6 path labels (IUPAC P-22.2.1 / P-25): hetero=1 … 7a.
FUSED56_LABELS: tuple[str, ...] = (
    # "1", "2", "3", "3a", "4", "5", "6", "7", "7a",
)
# Naphthalene / quinoline family path labels (P-25): 1…4a…8a.
NAPH_LABELS: tuple[str, ...] = (
    # "1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a",
)
ANTHRA_LABELS: tuple[str, ...] = (
    # "1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a",
)


def _carbo(sid: str, mode: str) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="carbocycle", stem_en=None, stem_zh=None,
        n_rings=1, ring="carbo", retained=False, fg_rank=0, numbering=pol,
    )


def _fused56(
    sid: str, stem_en: str | None, stem_zh: str | None, fg_rank: int = 0,
    *, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fused56_fixed", standard_path=FUSED56_LABELS)
    return ScaffoldSpec(
        id=sid, naming_class="fused56", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring="hetero", retained=retained, fg_rank=fg_rank,
        numbering=pol,
    )


def _naph(
    sid: str, stem_en: str | None, stem_zh: str | None, *,
    ring: str = "hetero", fg_rank: int = 0, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="naph_family", standard_path=NAPH_LABELS)
    return ScaffoldSpec(
        id=sid, naming_class="naph_family", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=2, ring=ring, retained=retained, fg_rank=fg_rank, numbering=pol,
    )


def _monohetero(
    sid: str, stem_en: str, stem_zh: str, *, fg_rank: int = 0,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode="fixed_hetero")
    return ScaffoldSpec(
        id=sid, naming_class="monohetero", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="hetero", retained=True, fg_rank=fg_rank,
        numbering=pol,
    )


def _mono_carbo(
    sid: str, stem_en: str, stem_zh: str, *,
    mode: str = "fixed_roles", fg_rank: int = 0, retained: bool = True,
) -> ScaffoldSpec:
    pol = NumberingPolicy(mode=mode)
    return ScaffoldSpec(
        id=sid, naming_class="mono_carbo", stem_en=stem_en, stem_zh=stem_zh,
        n_rings=1, ring="carbo", retained=retained, fg_rank=fg_rank, numbering=pol,
    )




_ALL_SPECS: tuple[ScaffoldSpec, ...] = (
    _mono_carbo("benzene", "benzene", "苯"),
    _monohetero("pyridine", "pyridine", "吡啶"),
    _naph("naphthalene", "naphthalene", "萘", ring="carbo"),
    _fused56("indole", "1H-indole", "吲哚"),
)
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}
_IDENTITIES: dict[str, ScaffoldIdentity] = {s.id: s.identity for s in _ALL_SPECS}


def get_identity(spec_id: str) -> ScaffoldIdentity | None:
    return _IDENTITIES.get(spec_id)


def all_identities() -> tuple[ScaffoldIdentity, ...]:
    return tuple(_IDENTITIES.values())


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    return _BY_ID.get(spec_id)



def all_specs() -> tuple[ScaffoldSpec, ...]:
    return _ALL_SPECS


def numbering_scaffold_facts(spec_id: str | None, atom_count: int) -> dict | None:
    """Materialize pure parent facts; ``relative_stereo`` reserves ring-face constraints."""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.numbering.materialize_plan:
        return None
    labels = spec.numbering.standard_path
    if spec.id in {"cycloalkane", "cycloalkane_polycarboxylic"}:
        labels = tuple(str(i) for i in range(1, atom_count + 1)) if 3 <= atom_count <= 10 else ()
    return None if not labels or len(labels) != atom_count else {
        "scaffold_id": spec.id, "labels": labels, "relative_stereo": None,
    }


def fused56_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids that use fused56 1…7a labels."""
    return kind_ids_for("fused56")


def naph_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids that use naph 1…8a labels."""
    return kind_ids_for("naph_family")


def monohetero_kind_ids() -> frozenset[str]:
    """Kind / scaffold ids for mono-hetero retained scaffolds."""
    return kind_ids_for("monohetero")


def kind_ids_for(naming_class: str) -> frozenset[str]:
    """Ids for a naming_class (fused56 / naph_family / monohetero / …)."""
    return frozenset(s.id for s in _ALL_SPECS if s.naming_class == naming_class)
