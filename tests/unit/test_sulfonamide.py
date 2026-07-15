# IUPAC: P-65.3
# Layer: L1,L2,L5
"""Simple mono-sulfonamide functional parents (P-65.3).

Scope (first cut):
- unsubstituted alkanesulfonamide R-SO2-NH2, R = n-alkyl C1–C4
- unsubstituted / simple arenesulfonamide Ar-SO2-NH2 (Ph ≤2 halo/Me)
- N-monoaryl alkanesulfonamide (S: n-alkyl C1–C4; N: simple Ph ≤2 halo/Me)
- N-monocycloalkyl arenesulfonamide (S: simple Ph; N: cycloalkyl C3–C6)

Out of scope: di-sulfonamide, N,N-disub, piperazine-sulfonyl, heteroaryl,
sulfonate ester, sulfonyl chloride.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # unsubstituted alkanesulfonamide
    ("CS(=O)(=O)N", "methanesulfonamide", "甲磺酰胺"),
    ("CCS(=O)(=O)N", "ethanesulfonamide", "乙磺酰胺"),
    # unsubstituted arenesulfonamide
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    # N-monoaryl alkanesulfonamide (benchmark gold)
    (
        "BrC1=C(C=C(C=C1)NS(=O)(=O)CC)C",
        "N-(4-bromo-3-methylphenyl)ethanesulfonamide",
        "N-(4-溴-3-甲基苯基)乙磺酰胺",
    ),
    # N-monocycloalkyl arenesulfonamide (benchmark gold)
    (
        "ClC1=C(C=CC(=C1)Cl)S(=O)(=O)NC1CCCC1",
        "2,4-dichloro-N-cyclopentylbenzenesulfonamide",
        "2,4-二氯-N-环戊基苯磺酰胺",
    ),
    # negative near-miss: must not become sulfonamide / keep correct names
    ("CS(C)=O", "dimethyl sulfoxide", "二甲基亚砜"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
]


# Must not be named as sulfonamide
NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),
    ("CC(=O)N", "sulfonamide"),
    ("c1ccccc1N", "sulfonamide"),
    ("CS(=O)(=O)C", "sulfonamide"),  # sulfone
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_simple_sulfonamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰胺" not in zh


def test_methanesulfonamide_not_methane() -> None:
    r = SMILESNNamer().name("CS(=O)(=O)N")
    assert r.success
    assert normalize_en(r.en) == "methanesulfonamide"
    assert normalize_en(r.en) != "methane"
