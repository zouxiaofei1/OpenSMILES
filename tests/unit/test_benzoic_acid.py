# IUPAC: P-65.1.1.1
# Layer: L2,L3,L4,L5
"""Retained parent name benzoic acid (benzene + one carboxyl on ring carbon).

Parent = benzene with principal characteristic COOH (benzoic acid).
Allow 0–2 extra ring halo / methyl / phenolic hydroxy; COOH attach = locant 1.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted retained parent
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    # positive: monohalo / monomethyl
    ("OC(=O)c1ccc(Cl)cc1", "4-chlorobenzoic acid", "4-氯苯甲酸"),
    ("Cc1cccc(C(=O)O)c1", "3-methylbenzoic acid", None),
    # positive: phenolic hydroxy as prefix
    ("OC(=O)c1ccc(O)cc1", "4-hydroxybenzoic acid", "4-羟基苯甲酸"),
    ("O=C(O)c1cccc(O)c1", "3-hydroxybenzoic acid", None),
    # negative: phenol, acetic acid, benzene, diacid, ethanol
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("c1ccccc1", "benzene", "苯"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzoic_not_chain_acid() -> None:
    """Aromatic carboxyl must not expand ring into chain acid parent."""
    r = SMILESNNamer().name("OC(=O)c1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzoic acid"
    assert "heptanoic" not in en
    assert "hexanoic" not in en
