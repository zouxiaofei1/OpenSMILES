# IUPAC: P-63.2.1
# Layer: L1,L2,L5
"""Open-chain simple dialkyl sulfides (alkyl alkyl sulfide / 硫醚).

Exactly one sulfide S (two C neighbors, no H). Both arms unsubstituted linear
alkyl C1–C4. Symmetric: diethyl sulfide / 二乙硫醚. Asymmetric alphabetical:
ethyl methyl sulfide / 乙基甲基硫醚.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: asymmetric alphabetical alkyl alkyl sulfide
    ("CCSC", "ethyl methyl sulfide", "乙基甲基硫醚"),
    ("CCCSC", "methyl propyl sulfide", "甲基丙基硫醚"),
    # positive: symmetric retained di-alkyl sulfide
    ("CSC", "dimethyl sulfide", "二甲硫醚"),
    ("CCSCC", "diethyl sulfide", "二乙硫醚"),
    ("CCCSCCC", "dipropyl sulfide", "二丙硫醚"),
    # negative: thiol / ether / alcohol must not break
    ("CCS", "ethanethiol", "乙硫醇"),
    ("CCOCC", "diethyl ether", "二乙基醚"),
    ("CCO", "ethanol", "乙醇"),
]

# sulfoxide / thioester must not be named as dialkyl sulfide
NEG_NOT_SULFIDE = [
    ("CS(C)=O", "dimethyl sulfide"),
    ("CC(=O)SC", "methyl sulfide"),
]


@pytest.mark.parametrize("smiles,forbidden_en", NEG_NOT_SULFIDE)
def test_not_dialkyl_sulfide(smiles: str, forbidden_en: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert en != normalize_en(forbidden_en)
    assert "sulfide" not in en
    assert "硫醚" not in normalize_zh(r.zh)



@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_dialkyl_sulfide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_ethyl_methyl_sulfide_not_ethane() -> None:
    """CCSC must not collapse to ethane parent."""
    r = SMILESNNamer().name("CCSC")
    assert r.success
    assert normalize_en(r.en) == "ethyl methyl sulfide"
    assert normalize_en(r.en) != "ethane"
