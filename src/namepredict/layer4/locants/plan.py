"""NumberingPlan：纯 atom→locant 标签（L4）。"""
from __future__ import annotations

from dataclasses import dataclass


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
    """由对齐的 atom_order 与 labels 构建 NumberingPlan。"""
    order, labs = tuple(atom_order), tuple(labels)
    a2l, l2a = _maps(order, labs)
    return NumberingPlan(
        scaffold_id, order, labs, a2l, l2a,
        frozenset(sub_atoms), tuple(constraints_applied),
    )


def locant(plan: NumberingPlan, atom: int) -> str | None:
    """返回 atom 的 locant 标签；不在 plan 中则返回 None。"""
    return plan.atom_to_label.get(atom)
