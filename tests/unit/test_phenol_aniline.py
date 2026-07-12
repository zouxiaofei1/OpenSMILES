# IUPAC: P-63.1.4 / P-62.2.1.1.1
# Layer: L2,L4,L5
"""Retained parent names phenol / aniline (aromatic mono-OH / mono-NH2).

Parent = benzene with principal characteristic OH (phenol) or NH2 (aniline).
Allow 0–2 extra ring halo and/or methyl substituents; FG fixed as locant 1.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted retained parents
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    # positive: monohalo / monomethyl
    ("Oc1ccc(Cl)cc1", "4-chlorophenol", "4-氯苯酚"),
    ("Cc1ccccc1O", "2-methylphenol", "2-甲基苯酚"),
    ("Nc1ccc(C)cc1", "4-methylaniline", "4-甲基苯胺"),
    # positive: dimethyl phenol
    ("Cc1cc(C)cc(O)c1", "3,5-dimethylphenol", "3,5-二甲基苯酚"),
    # negative: bare / monohalo benzene, chain alcohol, saturated cyclo
    ("c1ccccc1", "benzene", "苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("CCO", "ethanol", "乙醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("NC1CCCCC1", "cyclohexanamine", "环己胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_phenol_aniline(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_phenol_not_chain_alcohol() -> None:
    """Aromatic OH must not expand ring into chain alcohol parent."""
    r = SMILESNNamer().name("Oc1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "phenol"
    assert "hexanol" not in en
    assert "nonan" not in en
