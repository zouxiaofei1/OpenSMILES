# IUPAC: P-65.6.3.5.1 / P-66.1.5.1
# Layer: L2,L4,L5
"""Saturated monocyclic lactones / lactams as heterocyclic pseudoketones.

PIN method (1): oxolan-2-one / pyrrolidin-2-one etc. (Hantzsch–Widman / retained
stem + -one). Hetero = 1; carbonyl locant usually 2. First-cut: unsubstituted
or simple mono-methyl/halo on ring C; N-unsub lactams only; 5/6-membered mono-O
or mono-N cores.

Must not steal open-chain ester/amide, plain sat_hetero, or cycloketone parents.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted lactones
    ("O=C1CCCO1", "oxolan-2-one", "氧杂环戊烷-2-酮"),
    ("O=C1CCCCO1", "oxan-2-one", "氧杂环己烷-2-酮"),
    # positive: unsubstituted lactams (N–H)
    ("O=C1CCCN1", "pyrrolidin-2-one", "吡咯烷-2-酮"),
    ("O=C1CCCCN1", "piperidin-2-one", "哌啶-2-酮"),
    # positive: simple mono-alkyl / mono-halo on ring C (hetero=1 → lowest set)
    ("CC1CCC(=O)O1", "5-methyloxolan-2-one", "5-甲基氧杂环戊烷-2-酮"),
    ("CC1CCC(=O)N1", "5-methylpyrrolidin-2-one", "5-甲基吡咯烷-2-酮"),
    ("CCC1CCC(=O)O1", "5-ethyloxolan-2-one", "5-乙基氧杂环戊烷-2-酮"),
    ("CCCC1COC(=O)C1", "4-propyloxolan-2-one", "4-丙基氧杂环戊烷-2-酮"),
    ("ClC1CCC(=O)O1", "5-chlorooxolan-2-one", "5-氯氧杂环戊烷-2-酮"),
    # negative: open-chain ester / amide keep open-chain names
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    # negative: plain sat_hetero / cycloketone must not break
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
    ("C1CCNC1", "pyrrolidine", "吡咯烷"),
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_one(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_lactone_not_open_ester_or_methane() -> None:
    """γ-Butyrolactone must not collapse to methane / open ester."""
    r = SMILESNNamer().name("O=C1CCCO1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "oxolan-2-one"
    assert "methane" not in en
    assert "formate" not in en
    assert "acetate" not in en


def test_lactam_not_formamide() -> None:
    """2-Pyrrolidone must not collapse to N-methylformamide."""
    r = SMILESNNamer().name("O=C1CCCN1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "pyrrolidin-2-one"
    assert "formamide" not in en
    assert "acetamide" not in en
