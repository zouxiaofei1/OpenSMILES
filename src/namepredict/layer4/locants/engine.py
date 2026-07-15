"""choose_numbering: pick best NumberingPlan under mode constraints (L4)."""
from __future__ import annotations

from namepredict.layer4.locants.constraints import (
    constraint_key,
    constraints_applied,
)
from namepredict.layer4.locants.generate import candidates_for, labels_for
from namepredict.layer4.locants.plan import NumberingPlan, make_plan


def _key(
    cand: list[int],
    mode: str,
    double_bonds,
    sub_attach,
) -> tuple:
    return constraint_key(
        cand, mode, double_bonds=double_bonds, sub_attach=sub_attach,
    )


def _best_candidate(
    cands: list[list[int]],
    mode: str,
    double_bonds,
    sub_attach,
) -> list[int]:
    best = cands[0]
    best_k = _key(best, mode, double_bonds, sub_attach)
    for cand in cands[1:]:
        k = _key(cand, mode, double_bonds, sub_attach)
        if k < best_k:
            best, best_k = cand, k
    return best


def _to_plan(
    order: list[int],
    scaffold_id: str,
    mode: str,
    sub_attach: list[int] | None,
) -> NumberingPlan:
    subs = frozenset(sub_attach or ())
    return make_plan(
        scaffold_id, order, labels_for(len(order)),
        sub_atoms=subs,
        constraints_applied=constraints_applied(mode),
    )


def choose_numbering(
    chain: list[int],
    mode: str,
    *,
    double_bonds: list[tuple[int, int]] | None = None,
    sub_attach: list[int] | None = None,
    scaffold_id: str = "",
) -> NumberingPlan:
    """Select ring numbering plan for carbocycle_free / poly_unsat."""
    cands = candidates_for(mode, list(chain))
    if not cands:
        return _to_plan([], scaffold_id, mode, sub_attach)
    best = _best_candidate(cands, mode, double_bonds, sub_attach)
    return _to_plan(best, scaffold_id, mode, sub_attach)
