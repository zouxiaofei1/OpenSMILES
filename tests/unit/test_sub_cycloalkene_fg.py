# IUPAC: P-31.1 / P-22.1.1 / P-63.1.1 / P-64.2.1 / P-14.3.4
# Layer: L2,L4,L5
"""Alkyl/halo monocyclic monoalkenes + cycloalkenol / cycloalkenone.

1) Parent cycloalkene with claimable n-alkyl / ring halo:
   double bond fixed at 1; lowest sub set → 1-methylcyclohexene.
2) Ring mono-OH + endocyclic C=C: FG@1, lowest ene → cyclohex-2-en-1-ol.
3) Ring mono-ketone + endocyclic C=C: FG@1, lowest ene → cyclohex-2-en-1-one.
Unsubstituted cyclohexene / sat mono cyclo FG / open alkenol·one stay intact.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: alkyl / halo cycloalkene
    ("CC1=CCCCC1", "1-methylcyclohexene", "1-甲基环己烯"),
    ("CCC1=CCCCC1", "1-ethylcyclohexene", "1-乙基环己烯"),
    ("CC1=CCCC1", "1-methylcyclopentene", "1-甲基环戊烯"),
    ("CC1CCC=CC1", "4-methylcyclohexene", "4-甲基环己烯"),
    ("ClC1=CCCCC1", "1-chlorocyclohexene", "1-氯环己烯"),
    # positive: cycloalkenol
    ("OC1C=CCCC1", "cyclohex-2-en-1-ol", "环己-2-烯-1-醇"),
    ("OC1CC=CCC1", "cyclohex-3-en-1-ol", "环己-3-烯-1-醇"),
    ("OC1C=CCC1", "cyclopent-2-en-1-ol", "环戊-2-烯-1-醇"),
    # positive: cycloalkenone
    ("O=C1C=CCCC1", "cyclohex-2-en-1-one", "环己-2-烯-1-酮"),
    ("O=C1CC=CCC1", "cyclohex-3-en-1-one", "环己-3-烯-1-酮"),
    # negative: unsubstituted / sat mono / open-chain must not regress
    ("C1=CCCCC1", "cyclohexene", "环己烯"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("O=C1CCCCC1", "cyclohexanone", "环己酮"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("CC1CCCCC1O", "2-methylcyclohexan-1-ol", "2-甲基环己-1-醇"),
    ("CC=CCO", "but-2-en-1-ol", "丁-2-烯-1-醇"),
    ("CC(=O)C=C", "but-3-en-2-one", "丁-3-烯-2-酮"),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sub_cycloalkene_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methylcyclohexene_not_propene() -> None:
    """Alkyl cycloalkene must not collapse to open-chain alkene."""
    r = SMILESNNamer().name("CC1=CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-methylcyclohexene"
    assert en != "propene"
    assert "cyclohexene" in en


def test_cyclohexenol_not_ethene() -> None:
    """Cycloalkenol must not be named as ethene / open alkenol wrongly."""
    r = SMILESNNamer().name("OC1C=CCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohex-2-en-1-ol"
    assert en != "ethene"
    assert "en-1-ol" in en
