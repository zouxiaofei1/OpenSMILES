# IUPAC: P-63.3
# Layer: L1,L2,L5
"""Open-chain simple dialkyl sulfoxides (alkyl alkyl sulfoxide / 亚砜).

Exactly one sulfoxide S (one =O, two C neighbors, not sulfone). Both arms
unsubstituted linear alkyl C1–C4. Symmetric: dimethyl sulfoxide / 二甲基亚砜
(gold). Asymmetric alphabetical EN: ethyl methyl sulfoxide / 乙基甲基亚砜.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: symmetric dialkyl sulfoxide (gold ZH uses 基)
    ("CS(C)=O", "dimethyl sulfoxide", "二甲基亚砜"),
    ("CCS(=O)CC", "diethyl sulfoxide", "二乙基亚砜"),
    ("CCCS(=O)CCC", "dipropyl sulfoxide", "二丙基亚砜"),
    # positive: asymmetric alphabetical alkyl alkyl sulfoxide
    ("CS(=O)CC", "ethyl methyl sulfoxide", "乙基甲基亚砜"),
    ("CS(=O)CCC", "methyl propyl sulfoxide", "甲基丙基亚砜"),
    # negative: sulfide must not regress
    ("CSC", "dimethyl sulfide", "二甲硫醚"),
    ("CCSCC", "diethyl sulfide", "二乙硫醚"),
]

# sulfone / sulfide must not be named as dialkyl sulfoxide
NEG_NOT_SULFOXIDE = [
    ("CS(=O)(=O)C", "sulfoxide"),
    ("CSC", "sulfoxide"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_dialkyl_sulfoxide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFOXIDE)
def test_not_dialkyl_sulfoxide(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "亚砜" not in zh


def test_dmso_not_methane() -> None:
    """CS(C)=O must not collapse to methane parent."""
    r = SMILESNNamer().name("CS(C)=O")
    assert r.success
    assert normalize_en(r.en) == "dimethyl sulfoxide"
    assert normalize_en(r.en) != "methane"
