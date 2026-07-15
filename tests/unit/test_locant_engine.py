# IUPAC: P-14.4 / P-31.1
# Layer: L4
"""Ring locant engine: carbocycle_free + poly_unsat (no namer wire)."""
from __future__ import annotations

from namepredict.layer4.locants import (
    NumberingPlan,
    choose_numbering,
    locant,
    locant_int,
    ring_candidates,
)
from namepredict.layer4.locants.constraints import constraint_key
from namepredict.layer4.locants.generate import labels_for


# --- helpers (synthetic ring, no mol) ---

def _subs(attach: int) -> list[dict]:
    return [{"attach_idx": attach, "en": "methyl"}]


def _plan_free(chain, sub_attach=None, scaffold_id="cycloalkane"):
    return choose_numbering(
        chain, "carbocycle_free",
        sub_attach=sub_attach, scaffold_id=scaffold_id,
    )


def _plan_poly(chain, double_bonds, sub_attach=None, scaffold_id="cyclopolyene"):
    return choose_numbering(
        chain, "poly_unsat",
        double_bonds=double_bonds, sub_attach=sub_attach,
        scaffold_id=scaffold_id,
    )


# --- generate ---

def test_ring_candidates_n6_count():
    """n=6 free mode: 6 rotations × 2 directions = 12."""
    chain = [0, 1, 2, 3, 4, 5]
    cands = ring_candidates(chain)
    assert len(cands) == 12
    assert [0, 1, 2, 3, 4, 5] in cands
    assert [0, 5, 4, 3, 2, 1] in cands


def test_labels_for_sequential():
    assert labels_for(6) == ("1", "2", "3", "4", "5", "6")
    assert labels_for(5) == ("1", "2", "3", "4", "5")


# --- carbocycle_free: methylcyclohexane ---

def test_methylcyclohexane_sub_at_locant_1():
    """Single ring methyl: lowest set places attach at locant 1."""
    # chain order arbitrary; methyl on atom 3
    chain = [10, 11, 12, 13, 14, 15]
    plan = _plan_free(chain, sub_attach=[13])
    assert isinstance(plan, NumberingPlan)
    assert locant(plan, 13) == "1"
    assert locant_int(plan, 13) == 10
    assert plan.atom_order[0] == 13
    assert plan.scaffold_id == "cycloalkane"
    assert "substituent" in plan.constraints_applied


def test_dimethyl_lowest_set_1_3_not_1_4():
    """Two subs: prefer 1,3 over 1,4 when both possible by rotation."""
    # atoms 0..5; subs at 0 and 2 → under free numbering → 1,3
    chain = [0, 1, 2, 3, 4, 5]
    plan = _plan_free(chain, sub_attach=[0, 2])
    locs = sorted(int(locant(plan, a)) for a in (0, 2))
    assert locs == [1, 3]


def test_free_no_sub_identity_order():
    """No substituents: first candidate / original direction wins stably."""
    chain = [0, 1, 2, 3, 4, 5]
    plan = _plan_free(chain, sub_attach=None)
    assert plan.atom_order == (0, 1, 2, 3, 4, 5)
    assert plan.labels == ("1", "2", "3", "4", "5", "6")


# --- poly_unsat: 1,3 vs 1,4 diene ---

def test_cyclohexa_1_3_diene_ene_locants():
    """Double bonds at 0=1 and 2=3 → ene locants (1,3), not (1,5)."""
    chain = [0, 1, 2, 3, 4, 5]
    bonds = [(0, 1), (2, 3)]
    plan = _plan_poly(chain, double_bonds=bonds)
    assert plan.atom_order[0] in (0, 2)  # either double start at 1
    # lower endpoint locants of each C=C, sorted
    ene = _ene_mins(plan, bonds)
    assert ene == (1, 3)
    assert "unsaturation" in plan.constraints_applied


def test_cyclohexa_1_4_diene_ene_locants():
    """Double bonds at 0=1 and 3=4 → ene locants (1,4)."""
    chain = [0, 1, 2, 3, 4, 5]
    bonds = [(0, 1), (3, 4)]
    plan = _plan_poly(chain, double_bonds=bonds)
    ene = _ene_mins(plan, bonds)
    assert ene == (1, 4)


def test_poly_unsat_distinguishes_13_vs_14():
    """Same ring size, different bond placement → different ene key."""
    chain = [0, 1, 2, 3, 4, 5]
    p13 = _plan_poly(chain, double_bonds=[(0, 1), (2, 3)])
    p14 = _plan_poly(chain, double_bonds=[(0, 1), (3, 4)])
    assert _ene_mins(p13, [(0, 1), (2, 3)]) == (1, 3)
    assert _ene_mins(p14, [(0, 1), (3, 4)]) == (1, 4)
    # constraint keys must differ on primary ene component
    k13 = constraint_key(
        p13.atom_order, "poly_unsat",
        double_bonds=[(0, 1), (2, 3)], sub_attach=None,
    )
    k14 = constraint_key(
        p14.atom_order, "poly_unsat",
        double_bonds=[(0, 1), (3, 4)], sub_attach=None,
    )
    assert k13[0] == (1, 3)
    assert k14[0] == (1, 4)
    assert k13[0] < k14[0]


def test_poly_unsat_sub_secondary():
    """Ene set primary; when ene tied, substituent set breaks ties."""
    # 1,4-diene is symmetric; methyl should go to lowest free locant
    chain = [0, 1, 2, 3, 4, 5]
    bonds = [(0, 1), (3, 4)]
    plan = _plan_poly(chain, double_bonds=bonds, sub_attach=[2])
    ene = _ene_mins(plan, bonds)
    assert ene == (1, 4)
    # sub on atom between the doubles can be 2, 3, 5, or 6 depending on orient
    sub_loc = int(locant(plan, 2))
    assert sub_loc in (2, 3, 5, 6)
    # but must be the lowest achievable among candidates with ene (1,4)
    assert sub_loc == min(
        _sub_loc_for_cand(c, 2)
        for c in ring_candidates(chain)
        if _ene_mins_order(c, bonds) == (1, 4)
    )


# --- constraint_key pure ---

def test_constraint_key_free_subs():
    order = (13, 14, 15, 10, 11, 12)
    key = constraint_key(order, "carbocycle_free", sub_attach=[13])
    assert key[0] == (1,)  # sub at locant 1


def test_constraint_key_poly_ene_primary():
    order = (0, 1, 2, 3, 4, 5)
    key = constraint_key(
        order, "poly_unsat",
        double_bonds=[(0, 1), (2, 3)], sub_attach=None,
    )
    assert key[0] == (1, 3)


# --- helpers used by tests ---

def _ene_mins(plan: NumberingPlan, bonds) -> tuple[int, ...]:
    return _ene_mins_order(plan.atom_order, bonds)


def _bond_loc_order(order, a: int, b: int) -> int:
    idx = {x: i + 1 for i, x in enumerate(order)}
    pa, pb = idx[a], idx[b]
    lo, hi = min(pa, pb), max(pa, pb)
    n = len(order)
    if hi - lo == 1:
        return lo
    if lo == 1 and hi == n:
        return n
    return 999


def _ene_mins_order(order, bonds) -> tuple[int, ...]:
    return tuple(sorted(_bond_loc_order(order, a, b) for a, b in bonds))


def _sub_loc_for_cand(cand, atom: int) -> int:
    return list(cand).index(atom) + 1
