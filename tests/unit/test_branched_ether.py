# IUPAC: P-63.2.2
# Layer: L2, L5
"""Branched dialkyl ethers: isopropyl / hexafluoroisopropyl functional class.

Linear C1–C4 arms keep alkoxyalkane / retained sym ether. Branched specials
use functional-class names (… methyl ether / …基甲醚).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # hexafluoroisopropyl methyl ether (benchmark tiers-15987)
    (
        "COC(C(F)(F)F)C(F)(F)F",
        "hexafluoroisopropyl methyl ether",
        "六氟异丙基甲醚",
    ),
    (
        "FC(F)(F)C(OC)C(F)(F)F",
        "hexafluoroisopropyl methyl ether",
        "六氟异丙基甲醚",
    ),
    # unsubstituted isopropyl methyl ether
    ("COC(C)C", "isopropyl methyl ether", "异丙基甲醚"),
    # linear regressions
    ("COC", "dimethyl ether", "二甲基醚"),
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCOCC", "diethyl ether", "二乙基醚"),
    # negatives
    ("CCO", "ethanol", "乙醇"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_branched_ether(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_hfip_not_propane_collapse() -> None:
    r = SMILESNNamer().name("COC(C(F)(F)F)C(F)(F)F")
    assert r.success
    en = normalize_en(r.en)
    assert "ether" in en
    assert "propane" not in en
    assert "hexafluoroisopropyl" in en
