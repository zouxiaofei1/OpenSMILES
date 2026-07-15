# IUPAC: P-14 / P-25
# Layer: L4
"""Adapt retained fused-ring chains → NumberingPlan; parity with legacy tables."""
from __future__ import annotations

import pytest

from namepredict.layer4.locants.adapt import (
    INDOLE_LABELS,
    INDOLE_LOCANTS,
    NAPH_LABELS,
    NAPH_LOCANTS,
    effective_sub_locant,
    plan_from_chain,
)
from namepredict.layer4.locants.plan import locant
from namepredict.layer4.numbering import (
    _FUSED56_KINDS,
    _NAPH_KINDS,
    _Q_KINDS,
    _indole_sub_locant,
    _naph_sub_locant,
    _sub_locant,
)

_INDOLE_KINDS = sorted(_FUSED56_KINDS)
_NAPH_ALL = sorted(_NAPH_KINDS | _Q_KINDS)


def _legacy_from_table(table: tuple, chain: list[int], attach: int) -> int:
    """Mirror legacy: table None → chain.index+1; missing → 0."""
    if attach not in chain or len(chain) != len(table):
        return chain.index(attach) + 1 if attach in chain else 0
    loc = table[chain.index(attach)]
    return loc if loc is not None else chain.index(attach) + 1


def test_indole_plan_labels_and_length():
    chain = list(range(9))
    plan = plan_from_chain(chain, "indole")
    assert plan is not None
    assert plan.atom_order == tuple(chain)
    assert plan.labels == INDOLE_LABELS
    assert locant(plan, chain[3]) == "3a"
    assert locant(plan, chain[8]) == "7a"


def test_naph_plan_labels_and_length():
    chain = list(range(10))
    plan = plan_from_chain(chain, "naphthalene")
    assert plan is not None
    assert plan.labels == NAPH_LABELS
    assert locant(plan, chain[4]) == "4a"
    assert locant(plan, chain[9]) == "8a"


@pytest.mark.parametrize("idx", range(9))
def test_indole_effective_parity_vs_legacy_table(idx: int):
    chain = list(range(9))
    atom = chain[idx]
    plan = plan_from_chain(chain, "indole")
    assert plan is not None
    expected = _legacy_from_table(INDOLE_LOCANTS, chain, atom)
    assert effective_sub_locant(plan, atom) == expected
    assert expected == _indole_sub_locant(chain, atom)


@pytest.mark.parametrize("idx", range(10))
def test_naph_effective_parity_vs_legacy_table(idx: int):
    chain = list(range(10))
    atom = chain[idx]
    plan = plan_from_chain(chain, "naphthalene")
    assert plan is not None
    expected = _legacy_from_table(NAPH_LOCANTS, chain, atom)
    assert effective_sub_locant(plan, atom) == expected
    assert expected == _naph_sub_locant(chain, atom)


@pytest.mark.parametrize("kind", _INDOLE_KINDS)
@pytest.mark.parametrize("idx", range(9))
def test_sub_locant_prefers_plan_indole_family(kind: str, idx: int):
    chain = list(range(100, 109))
    atom = chain[idx]
    plan = plan_from_chain(chain, kind)
    assert plan is not None
    assert _sub_locant(chain, atom, kind) == effective_sub_locant(plan, atom)
    assert _sub_locant(chain, atom, kind) == _indole_sub_locant(chain, atom)


@pytest.mark.parametrize("kind", _NAPH_ALL)
@pytest.mark.parametrize("idx", range(10))
def test_sub_locant_prefers_plan_naph_family(kind: str, idx: int):
    chain = list(range(200, 210))
    atom = chain[idx]
    plan = plan_from_chain(chain, kind)
    assert plan is not None
    assert _sub_locant(chain, atom, kind) == effective_sub_locant(plan, atom)
    assert _sub_locant(chain, atom, kind) == _naph_sub_locant(chain, atom)


def test_plan_none_for_non_fused_kind():
    assert plan_from_chain(list(range(6)), "benzene") is None
    assert plan_from_chain(list(range(9)), "alkane") is None


def test_plan_none_for_wrong_chain_length():
    assert plan_from_chain(list(range(8)), "indole") is None
    assert plan_from_chain(list(range(9)), "naphthalene") is None


def test_effective_unknown_atom_none():
    plan = plan_from_chain(list(range(9)), "indole")
    assert plan is not None
    assert effective_sub_locant(plan, 999) is None


def test_label_to_atom_unique_indole():
    plan = plan_from_chain(list(range(9)), "benzothiazole")
    assert plan is not None
    assert len(plan.label_to_atom) == 9
    assert plan.label_to_atom["3a"] == 3
    assert plan.label_to_atom["7a"] == 8
