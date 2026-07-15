# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Simple mono-sulfonyl chloride functional parents (P-65.3).

Scope (first cut):
- open-chain alkanesulfonyl chloride R–SO2–Cl, R = n-alkyl C1–C4
  with terminal F/CF3 on the chain (L3 halo prefixes)
- unfused arenesulfonyl chloride Ar–SO2–Cl, Ph with ≤2 Me/halo

Out of scope: sulfonyl fluoride/bromide, sulfonic anhydride, complex
branching beyond terminal CF3/F, poly-sulfonyl chlorides.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # unsubstituted alkanesulfonyl chloride
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride", "甲磺酰氯"),
    ("CCS(=O)(=O)Cl", "ethanesulfonyl chloride", "乙磺酰氯"),
    # terminal CF3 on propane chain (benchmark gold)
    (
        "FC(CCS(=O)(=O)Cl)(F)F",
        "3,3,3-trifluoropropanesulfonyl chloride",
        "3,3,3-三氟丙磺酰氯",
    ),
    # unsubstituted / simple arenesulfonyl chloride
    ("c1ccc(S(=O)(=O)Cl)cc1", "benzenesulfonyl chloride", "苯磺酰氯"),
    (
        "Cc1ccc(S(=O)(=O)Cl)cc1",
        "4-methylbenzenesulfonyl chloride",
        "4-甲基苯磺酰氯",
    ),
    # negative near-miss: must keep existing correct names
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    ("CC(=O)Cl", "acetyl chloride", "乙酰氯"),
    ("CS(C)=O", "dimethyl sulfoxide", "二甲基亚砜"),
    (
        "CC1=CC=C(C=C1)S(=O)(=O)OCCCC",
        "butyl 4-methylbenzenesulfonate",
        "对甲苯磺酸正丁酯",
    ),
]


# Must not be named as sulfonyl chloride
NEG_NOT_SULFONYL_CHLORIDE = [
    ("c1ccc(S(=O)(=O)N)cc1", "sulfonyl chloride"),
    ("CC(=O)Cl", "sulfonyl chloride"),
    ("CS(C)=O", "sulfonyl chloride"),
    ("CC1=CC=C(C=C1)S(=O)(=O)OCCCC", "sulfonyl chloride"),
    ("CS(=O)(=O)N", "sulfonyl chloride"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_sulfonyl_chloride(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONYL_CHLORIDE)
def test_not_sulfonyl_chloride(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰氯" not in zh


def test_methanesulfonyl_chloride_not_methane() -> None:
    r = SMILESNNamer().name("CS(=O)(=O)Cl")
    assert r.success
    assert normalize_en(r.en) == "methanesulfonyl chloride"
    assert normalize_en(r.en) != "methane"
