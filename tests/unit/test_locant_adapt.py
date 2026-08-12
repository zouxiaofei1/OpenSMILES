# IUPAC: P-14 / P-25
# Layer: L4
"""Retained fused locants: plan_from_chain / effective_sub_locant single track."""
from __future__ import annotations

import pytest

from namepredict.layer2.specs import get_spec, numbering_scaffold_facts
from namepredict.layer4.locants.adapt import effective_sub_locant, plan_from_chain


def _plan(chain: list[int], kind: str):
    return plan_from_chain(chain, kind, numbering_scaffold_facts(kind, len(chain)))

from namepredict.layer4.locants.plan import locant
from namepredict.layer4.numbering import _amine_locant, _oh_locant, _sub_locant

_INDOLE_KINDS = ["indole", "indazole", "benzofuran", "benzofuranamine", "benzothiophene", "benzothiophenol", "benzothiazole", "benzothiazolamine", "benzoxazole", "benzoxazolamine", "benzimidazole", "benzimidazolamine"]
_NAPH_ALL = ["quinoline", "isoquinoline", "naphthalene"]
_ANTHRA_ALL = ["anthracene", "anthraquinone"]
INDOLE_LABELS = get_spec("indole").numbering.standard_path
NAPH_LABELS = get_spec("naphthalene").numbering.standard_path
ANTHRA_LABELS = get_spec("anthracene").numbering.standard_path
INDOLE_LOCANTS = (1, 2, 3, None, 4, 5, 6, 7, None)
NAPH_LOCANTS = (1, 2, 3, 4, None, 5, 6, 7, 8, None)
ANTHRA_LOCANTS = (1, 2, 3, 4, None, 10, None, 5, 6, 7, 8, None, 9, None)


def _legacy_from_table(table: tuple, chain: list[int], attach: int) -> int:
    """Mirror legacy table: None → chain.index+1; missing → 0."""
    if attach not in chain or len(chain) != len(table):
        return chain.index(attach) + 1 if attach in chain else 0
    loc = table[chain.index(attach)]
    return loc if loc is not None else chain.index(attach) + 1


def test_indole_plan_labels_and_length():
    chain = list(range(9))
    plan = _plan(chain, "indole")
    assert plan is not None
    assert plan.atom_order == tuple(chain)
    assert plan.labels == INDOLE_LABELS
    assert locant(plan, chain[3]) == "3a"
    assert locant(plan, chain[8]) == "7a"


def test_naph_plan_labels_and_length():
    chain = list(range(10))
    plan = _plan(chain, "naphthalene")
    assert plan is not None
    assert plan.labels == NAPH_LABELS
    assert locant(plan, chain[4]) == "4a"
    assert locant(plan, chain[9]) == "8a"


def test_anthra_plan_labels_and_length():
    chain = list(range(14))
    plan = _plan(chain, "anthraquinone")
    assert plan is not None
    assert plan.labels == ANTHRA_LABELS
    assert locant(plan, chain[4]) == "4a"
    assert locant(plan, chain[5]) == "10"


@pytest.mark.parametrize("idx", range(9))
def test_indole_effective_parity_vs_legacy_table(idx: int):
    chain = list(range(9))
    atom = chain[idx]
    plan = _plan(chain, "indole")
    assert plan is not None
    expected = _legacy_from_table(INDOLE_LOCANTS, chain, atom)
    assert effective_sub_locant(plan, atom) == expected
    assert _sub_locant(chain, atom, "indole", numbering_scaffold_facts("indole", len(chain))) == expected


@pytest.mark.parametrize("idx", range(10))
def test_naph_effective_parity_vs_legacy_table(idx: int):
    chain = list(range(10))
    atom = chain[idx]
    plan = _plan(chain, "naphthalene")
    assert plan is not None
    expected = _legacy_from_table(NAPH_LOCANTS, chain, atom)
    assert effective_sub_locant(plan, atom) == expected
    assert _sub_locant(chain, atom, "naphthalene", numbering_scaffold_facts("naphthalene", len(chain))) == expected


@pytest.mark.parametrize("idx", range(14))
def test_anthra_effective_parity_vs_legacy_table(idx: int):
    chain = list(range(14))
    atom = chain[idx]
    plan = _plan(chain, "anthraquinone")
    assert plan is not None
    expected = _legacy_from_table(ANTHRA_LOCANTS, chain, atom)
    assert effective_sub_locant(plan, atom) == expected
    assert _sub_locant(chain, atom, "anthraquinone", numbering_scaffold_facts("anthraquinone", len(chain))) == expected


@pytest.mark.parametrize("kind", _INDOLE_KINDS)
@pytest.mark.parametrize("idx", range(9))
def test_sub_locant_prefers_plan_indole_family(kind: str, idx: int):
    chain = list(range(100, 109))
    atom = chain[idx]
    plan = _plan(chain, kind)
    assert plan is not None
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == effective_sub_locant(plan, atom)
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == _legacy_from_table(INDOLE_LOCANTS, chain, atom)


@pytest.mark.parametrize("kind", _NAPH_ALL)
@pytest.mark.parametrize("idx", range(10))
def test_sub_locant_prefers_plan_naph_family(kind: str, idx: int):
    chain = list(range(200, 210))
    atom = chain[idx]
    plan = _plan(chain, kind)
    assert plan is not None
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == effective_sub_locant(plan, atom)
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == _legacy_from_table(NAPH_LOCANTS, chain, atom)


@pytest.mark.parametrize("kind", _ANTHRA_ALL)
@pytest.mark.parametrize("idx", range(14))
def test_sub_locant_prefers_plan_anthra_family(kind: str, idx: int):
    chain = list(range(300, 314))
    atom = chain[idx]
    plan = _plan(chain, kind)
    assert plan is not None
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == effective_sub_locant(plan, atom)
    assert _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain))) == _legacy_from_table(ANTHRA_LOCANTS, chain, atom)


# FG locants on retained fused parents must use the same plan path as _sub_locant.
_FG_CASES = [
    # (kind, chain_len, atom_key, locant_fn, chain_idx, expected_loc)
    ("benzothiophenol", 9, "oh_c_idx", _oh_locant, 1, 2),
    ("benzothiophenol", 9, "oh_c_idx", _oh_locant, 4, 4),
    ("benzothiophenol", 9, "oh_c_idx", _oh_locant, 6, 6),
    ("benzofuranamine", 9, "amine_c_idx", _amine_locant, 1, 2),
    ("benzothiazolamine", 9, "amine_c_idx", _amine_locant, 1, 2),
    ("benzoxazolamine", 9, "amine_c_idx", _amine_locant, 1, 2),
    ("benzimidazolamine", 9, "amine_c_idx", _amine_locant, 1, 2),
]


@pytest.mark.parametrize("kind,n,key,fn,idx,expected", _FG_CASES)
def test_fg_locant_matches_plan(kind, n, key, fn, idx, expected):
    chain = list(range(400, 400 + n))
    atom = chain[idx]
    oriented = {"kind": kind, "chain": chain, key: atom, "numbering_scaffold": numbering_scaffold_facts(kind, len(chain))}
    plan = _plan(chain, kind)
    assert plan is not None
    assert fn(oriented) == expected
    assert fn(oriented) == effective_sub_locant(plan, atom)
    assert fn(oriented) == _sub_locant(chain, atom, kind, numbering_scaffold_facts(kind, len(chain)))


def test_sub_locant_non_fused_uses_index_plus_one():
    """Negative: non-retained chain still uses chain.index + 1 (no plan)."""
    chain = [10, 20, 30, 40]
    assert plan_from_chain(chain, "alkane") is None
    assert plan_from_chain(chain, "alcohol") is None
    assert _sub_locant(chain, 30, "alkane") == 3
    assert _sub_locant(chain, 10, "alcohol") == 1
    oriented = {"kind": "alcohol", "chain": chain, "oh_c_idx": 20}
    assert _oh_locant(oriented) == 2


def test_plan_none_for_non_fused_kind():
    assert plan_from_chain(list(range(6)), "benzene") is None
    assert plan_from_chain(list(range(9)), "alkane") is None


def test_plan_none_for_wrong_chain_length():
    assert plan_from_chain(list(range(8)), "indole") is None
    assert _plan(list(range(9)), "naphthalene") is None
    assert plan_from_chain(list(range(10)), "anthraquinone") is None


def test_effective_unknown_atom_none():
    plan = _plan(list(range(9)), "indole")
    assert plan is not None
    assert effective_sub_locant(plan, 999) is None


def test_label_to_atom_unique_indole():
    plan = _plan(list(range(9)), "benzothiazole")
    assert plan is not None
    assert len(plan.label_to_atom) == 9
    assert plan.label_to_atom["3a"] == 3
    assert plan.label_to_atom["7a"] == 8
