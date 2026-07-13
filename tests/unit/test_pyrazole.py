# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained parent 1H-pyrazole (IUPAC P-22.2.1):
unsubstituted + at most one ring monomethyl / monohalo (F/Cl/Br/I).
NH = 1; other N = 2 (1,2-diazole); mono-sub uses standard locants.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("c1cn[nH]c1", "1H-pyrazole", "吡唑"),
    # positive: monomethyl (NH=1, N=2)
    ("Cc1cn[nH]c1", "4-methyl-1H-pyrazole", "4-甲基吡唑"),
    ("Cc1n[nH]cc1", "3-methyl-1H-pyrazole", "3-甲基吡唑"),
    ("Cc1ccn[nH]1", "5-methyl-1H-pyrazole", "5-甲基吡唑"),
    # positive: monohalo
    ("Clc1cn[nH]c1", "4-chloro-1H-pyrazole", "4-氯吡唑"),
    # negative: must not misclassify 1,3-diazole / diazine / mono-hetero5
    ("c1cnc[nH]1", "1H-imidazole", "咪唑"),
    ("c1cncnc1", "pyrimidine", "嘧啶"),
    ("c1ccoc1", "furan", "呋喃"),
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyrazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_pyrazole_not_imidazole() -> None:
    r = SMILESNNamer().name("c1cn[nH]c1")
    assert r.success
    en = normalize_en(r.en)
    assert "pyrazole" in en
    assert "imidazole" not in en


def test_methylpyrazole_locants() -> None:
    """With NH=1 and N=2, 4-methyl must not collapse to wrong parent."""
    r = SMILESNNamer().name("Cc1cn[nH]c1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "4-methyl-1h-pyrazole"
    assert "imidazole" not in en
