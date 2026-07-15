# IUPAC: P-65.3.2
# Layer: L1,L2,L5
"""Simple mono sulfonate ester functional parents (P-65.3.2).

Scope (first cut):
- alkyl arenesulfonate Ar–SO2–OR: unfused Ph with ≤2 Me/halo; R = n-alkyl C1–C4
- optional: unsubstituted methyl benzenesulfonate; methyl methanesulfonate

Out of scope: sulfonate salts ([O-]), sulfonic acids, sulfonamides, sulfonyl
chlorides, lactone-side tosylates, sulfate esters.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # benchmark gold: butyl tosylate
    (
        "CC1=CC=C(C=C1)S(=O)(=O)OCCCC",
        "butyl 4-methylbenzenesulfonate",
        "对甲苯磺酸正丁酯",
    ),
    # unsubstituted arenesulfonate
    (
        "c1ccc(S(=O)(=O)OC)cc1",
        "methyl benzenesulfonate",
        "苯磺酸甲酯",
    ),
    # ethyl tosylate
    (
        "CC1=CC=C(C=C1)S(=O)(=O)OCC",
        "ethyl 4-methylbenzenesulfonate",
        "对甲苯磺酸乙酯",
    ),
    # alkanesulfonate ester
    (
        "CS(=O)(=O)OC",
        "methyl methanesulfonate",
        "甲磺酸甲酯",
    ),
    # negative near-miss: must keep existing correct names
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    ("CC(=O)OCCCC", "butyl acetate", "乙酸丁酯"),
    ("CS(C)=O", "dimethyl sulfoxide", "二甲基亚砜"),
]


# Must not be named as sulfonate ester
NEG_NOT_SULFONATE = [
    ("c1ccc(S(=O)(=O)N)cc1", "sulfonate"),
    ("CC(=O)OCCCC", "sulfonate"),
    ("CS(C)=O", "sulfonate"),
    ("CS(=O)(=O)N", "sulfonate"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_sulfonate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONATE)
def test_not_sulfonate(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酸酯" not in zh
    assert not zh.endswith("磺酸")


def test_butyl_tosylate_not_butane() -> None:
    r = SMILESNNamer().name("CC1=CC=C(C=C1)S(=O)(=O)OCCCC")
    assert r.success
    assert normalize_en(r.en) == "butyl 4-methylbenzenesulfonate"
    assert normalize_en(r.en) != "butane"
