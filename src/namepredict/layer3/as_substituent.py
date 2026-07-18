"""name_as_substituent: cut submol → free-name pipeline → P-29 -yl form.

No host (benzamide gate / n_block extract) wiring — pure cut→pipeline→yl.
"""
from __future__ import annotations

from namepredict.layer2.submol_build import build_cut_submol
from namepredict.layer3.yl_form import yl_form
from namepredict.namer import _name_mol


def _locant_from_result(result, attach_new: int) -> int | None:
    chain = (result.meta or {}).get("parent_chain") or []
    if attach_new not in chain:
        return None
    return chain.index(attach_new) + 1


def _yl_from_sub(sub, *, depth: int) -> tuple[str, str, bool] | None:
    result = _name_mol(sub.mol, depth=depth)
    if not result.success or not result.en:
        return None
    loc = _locant_from_result(result, sub.attach_new)
    return None if loc is None else yl_form(result.en, result.zh, loc)


def name_as_substituent(
    mol, attach_old: int, atoms, *, depth: int = 0, max_depth: int = 4,
) -> tuple[str, str, bool] | None:
    """Cut atoms at attach_old, free-name the submol, emit -yl dual names."""
    if depth >= max_depth:
        return None
    sub = build_cut_submol(mol, frozenset(atoms), attach_old)
    return None if sub is None else _yl_from_sub(sub, depth=depth + 1)
