"""Orient sat_hetero_repl chains via multi_hetero locant engine (L4)."""
from __future__ import annotations

from namepredict.layer4.locants.engine import choose_numbering


def _attach_idxs(substituents: list) -> list[int]:
    return [s["attach_idx"] for s in substituents if s.get("attach_idx") is not None]


def orient_sat_hetero_repl(
    chain: list[int], parent: dict, substituents: list,
) -> list[int]:
    """multi_hetero: hetero set + a-order + sub set."""
    hs = list(parent.get("hetero_idxs") or [])
    plan = choose_numbering(
        chain, "multi_hetero", hetero_atoms=hs, hetero_z=parent.get("hetero_z"),
        sub_attach=_attach_idxs(substituents), scaffold_id="sat_hetero_repl",
    )
    return list(plan.atom_order) if plan.atom_order else chain
