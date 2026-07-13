# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained parent 1H-imidazole (IUPAC P-22.2.1):
unsubstituted + at most one ring monomethyl / monohalo (F/Cl/Br/I).
NH = 1; mono-sub takes lowest locant set (2- or 4-, not 5- when rotatable).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("c1cnc[nH]1", "1H-imidazole", "咪唑"),
    # positive: monomethyl / monohalo (NH=1 → 2- or 4-)
    ("Cc1ncc[nH]1", "2-methyl-1H-imidazole", "2-甲基-1H-咪唑"),
    ("Clc1ncc[nH]1", "2-chloro-1H-imidazole", "2-氯-1H-咪唑"),
    ("Cc1c[nH]cn1", "4-methyl-1H-imidazole", "4-甲基-1H-咪唑"),
    ("Clc1c[nH]cn1", "4-chloro-1H-imidazole", "4-氯-1H-咪唑"),
    # negative: must not regress mono-hetero5 / arene / open chain
    ("c1cc[nH]c1", "1H-pyrrole", "吡咯"),
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CCCCC", "pentane", "戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_imidazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_imidazole_not_alkane() -> None:
    r = SMILESNNamer().name("c1cnc[nH]1")
    assert r.success
    assert "imidazole" in normalize_en(r.en)
    assert normalize_en(r.en) != "ethane"


def test_methylimidazole_locant_not_5() -> None:
    """With NH fixed as 1, monomethyl must not be reported as 5-."""
    r = SMILESNNamer().name("Cc1c[nH]cn1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "4-methyl-1h-imidazole"
    assert "5-methyl" not in en
