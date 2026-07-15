"""NumberingPlan: pure atom→locant labels (L4). No choose_numbering engine."""
from __future__ import annotations

import re
from dataclasses import dataclass

# Locant total-order convention (bridgeheads):
#   plain "n"  → n * 10
#   "na"       → n * 10 + 1  (between n and n+1)
#   "nb"       → n * 10 + 2
# so "3" < "3a" < "4" as 30 < 31 < 40.
_LABEL_RE = re.compile(r"^(\d+)([a-z]?)$", re.I)
_SCALE = 10


@dataclass(frozen=True)
class NumberingPlan:
    scaffold_id: str
    atom_order: tuple[int, ...]
    labels: tuple[str, ...]
    atom_to_label: dict[int, str]
    label_to_atom: dict[str, int]
    sub_atoms: frozenset[int]
    constraints_applied: tuple[str, ...]


def _maps(
    atom_order: tuple[int, ...], labels: tuple[str, ...],
) -> tuple[dict[int, str], dict[str, int]]:
    a2l = dict(zip(atom_order, labels, strict=True))
    l2a = dict(zip(labels, atom_order, strict=True))
    return a2l, l2a


def make_plan(
    scaffold_id: str,
    atom_order: list[int] | tuple[int, ...],
    labels: list[str] | tuple[str, ...],
    sub_atoms: frozenset[int] = frozenset(),
    constraints_applied: tuple[str, ...] = (),
) -> NumberingPlan:
    """Build a NumberingPlan from aligned atom_order and labels."""
    order, labs = tuple(atom_order), tuple(labels)
    a2l, l2a = _maps(order, labs)
    return NumberingPlan(
        scaffold_id, order, labs, a2l, l2a,
        frozenset(sub_atoms), tuple(constraints_applied),
    )


def locant(plan: NumberingPlan, atom: int) -> str | None:
    """Return locant label for atom, or None if not in plan."""
    return plan.atom_to_label.get(atom)


def _label_to_int(label: str) -> int:
    m = _LABEL_RE.fullmatch(label.strip())
    if not m:
        raise ValueError(f"unsupported locant label: {label!r}")
    base = int(m.group(1))
    letter = m.group(2).lower()
    extra = (ord(letter) - ord("a") + 1) if letter else 0
    return base * _SCALE + extra


def locant_int(plan: NumberingPlan, atom: int) -> int | None:
    """Comparable int for atom's locant; '3a' ranks between '3' and '4'."""
    lab = locant(plan, atom)
    if lab is None:
        return None
    return _label_to_int(lab)
