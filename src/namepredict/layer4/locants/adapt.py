"""Adapt retained fused-ring chain + kind → NumberingPlan (zero behavior).

L4 must not import L2. Fused56 labels live here for Plan build; L2
ScaffoldSpec.standard_path must stay equal (contract-tested in unit tests).
"""
from __future__ import annotations

from namepredict.layer4.locants.plan import NumberingPlan, locant, make_plan

# Shared tables: single source for numbering.py (no drift).
NAPH_LOCANTS = (1, 2, 3, 4, None, 5, 6, 7, 8, None)
INDOLE_LOCANTS = (1, 2, 3, None, 4, 5, 6, 7, None)
NAPH_LABELS = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")
# Fused 5+6 path labels (must match scaffold.specs.FUSED56_LABELS).
INDOLE_LABELS = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")

_Q_KINDS = frozenset({
    "quinoline", "isoquinoline", "quinolinol", "quinolinecarboxylic",
})
_NAPH_KINDS = frozenset({"naphthalene", "naphthalenecarboxylic"})
_INDOLE_ORIENT_KINDS = (
    "indole", "indolecarboxylic", "indazole", "indazolecarbonitrile",
    "indazolecarbaldehyde", "benzofuran", "benzofuranamine", "benzothiophene",
    "benzothiophenol", "benzothiazole", "benzothiazolamine", "benzoxazole",
    "benzoxazolamine", "benzimidazole", "benzimidazolamine", "quinoline",
    "isoquinoline", "quinolinol", "quinolinecarboxylic",
)
_FUSED56_KINDS = frozenset(_INDOLE_ORIENT_KINDS) - _Q_KINDS


def _labels_for(kind: str | None) -> tuple[str, ...] | None:
    """Resolve labels from L4 tables only (no L2 import)."""
    if kind in _NAPH_KINDS or kind in _Q_KINDS:
        return NAPH_LABELS
    if kind in _FUSED56_KINDS:
        return INDOLE_LABELS
    return None


def plan_from_chain(chain: list[int], kind: str | None) -> NumberingPlan | None:
    """Build NumberingPlan from oriented fused chain; None if not retained fused."""
    labs = _labels_for(kind)
    if labs is None or len(chain) != len(labs):
        return None
    return make_plan(kind or "", tuple(chain), labs)


def effective_sub_locant(plan: NumberingPlan, atom: int) -> int | None:
    """Plain int matching legacy _sub_locant (bridgehead → chain index + 1)."""
    lab = locant(plan, atom)
    if lab is None:
        return None
    if lab.isdigit():
        return int(lab)
    try:
        return plan.atom_order.index(atom) + 1
    except ValueError:
        return None
