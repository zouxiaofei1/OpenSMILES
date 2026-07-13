# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
"""Saturated monocyclic hetero carboxylic acids (IUPAC P-65.1.1 / P-22.2.2).

Parent = sat monoheterocycle (sat_hetero core) with exactly one ring-C COOH.
Hetero = locant 1 (pair for dihetero); COOH locant from ring attach.
Chinese retained suffix 甲酸. Optional ≤1 simple ring prefix: hydroxy / methyl / halo.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted piperidinecarboxylic (all COOH positions)
    ("O=C(O)C1CCCCN1", "piperidine-2-carboxylic acid", "哌啶-2-甲酸"),
    ("O=C(O)C1CCCNC1", "piperidine-3-carboxylic acid", "哌啶-3-甲酸"),
    ("O=C(O)C1CCNCC1", "piperidine-4-carboxylic acid", "哌啶-4-甲酸"),
    # positive: pyrrolidine / piperazine / morpholine / oxolane
    ("O=C(O)C1CCCN1", "pyrrolidine-2-carboxylic acid", "吡咯烷-2-甲酸"),
    ("O=C(O)C1CCNC1", "pyrrolidine-3-carboxylic acid", "吡咯烷-3-甲酸"),
    ("O=C(O)C1CNCCN1", "piperazine-2-carboxylic acid", "哌嗪-2-甲酸"),
    ("O=C(O)C1COCCN1", "morpholine-3-carboxylic acid", "吗啉-3-甲酸"),
    ("O=C(O)C1CNCCO1", "morpholine-2-carboxylic acid", "吗啉-2-甲酸"),
    ("O=C(O)C1CCCO1", "oxolane-2-carboxylic acid", "氧杂环戊烷-2-甲酸"),
    # positive: optional mono hydroxy prefix on ring
    ("O=C(O)C1CCC(O)CN1", "5-hydroxypiperidine-2-carboxylic acid", "5-羟基哌啶-2-甲酸"),
    # negative: must not steal existing correct parents
    ("O=C(O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),
    ("C1CCNCC1", "piperidine", "哌啶"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("O=C(O)CCCCC", "hexanoic acid", "己酸"),
    ("O=C(O)c1ccncc1", "pyridine-4-carboxylic acid", "吡啶-4-甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_piperidine_2_carboxylic_not_hexanoic() -> None:
    """Ring COOH on piperidine must not collapse into open-chain hexanoic acid."""
    r = SMILESNNamer().name("O=C(O)C1CCCCN1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "piperidine-2-carboxylic acid"
    assert "hexanoic" not in en


def test_pyrrolidine_2_not_butanoic() -> None:
    """Pyrrolidine-2-carboxylic acid must not collapse to butanoic acid."""
    r = SMILESNNamer().name("O=C(O)C1CCCN1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "pyrrolidine-2-carboxylic acid"
    assert "butanoic" not in en
