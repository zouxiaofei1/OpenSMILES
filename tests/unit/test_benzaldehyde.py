# IUPAC: P-66.6.1
# Layer: L2,L4,L5
"""Retained parent name benzaldehyde (benzene + one formyl on ring carbon).

Parent = benzene with principal characteristic CHO (benzaldehyde).
Allow 0–2 extra ring halo / methyl / phenolic hydroxy; CHO attach = locant 1.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted retained parent
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    # positive: phenolic hydroxy as prefix
    ("O=Cc1ccc(O)cc1", "4-hydroxybenzaldehyde", None),
    # positive: hydroxy + methyl
    ("Cc1cccc(C=O)c1O", "2-hydroxy-3-methylbenzaldehyde", None),
    # negative: benzoic acid, acetaldehyde, phenol, benzene, ethanol
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("c1ccccc1", "benzene", "苯"),
    ("CCO", "ethanol", "乙醇"),
    # negative: acetophenone retained parent (not benzaldehyde)
    ("CC(=O)c1ccccc1", "acetophenone", "苯乙酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzaldehyde(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzaldehyde_not_chain_al() -> None:
    """Aromatic formyl must not expand ring into chain aldehyde parent."""
    r = SMILESNNamer().name("O=Cc1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzaldehyde"
    assert "heptanal" not in en
    assert "hexanal" not in en
