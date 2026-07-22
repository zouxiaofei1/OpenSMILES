# IUPAC: P-64.1.1
# Layer: L2,L4,L5
"""Retained parent name acetophenone (benzene + one acetyl on ring carbon).

Parent = benzene with principal characteristic acetyl Ph–C(=O)–CH3.
Allow 0–2 extra ring halo / methyl / phenolic hydroxy; acetyl attach = locant 1.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted retained parent
    ("CC(=O)c1ccccc1", "acetophenone", "苯乙酮"),
    # positive: ring chloro
    ("CC(=O)c1ccc(Cl)cc1", "4-chloroacetophenone", "4-氯苯乙酮"),
    # positive: phenolic hydroxy
    ("CC(=O)c1ccc(O)cc1", "4-hydroxyacetophenone", "4-羟基苯乙酮"),
    # negative: benzaldehyde, acetone, benzene, benzoic acid, ethanol
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("c1ccccc1", "benzene", "苯"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("CCO", "ethanol", "乙醇"),
    # negative: propiophenone must NOT be acetophenone (keep chain ketone)
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_acetophenone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_acetophenone_not_chain_ketone() -> None:
    """Aromatic acetyl must not expand ring into chain ketone parent."""
    r = SMILESNNamer().name("CC(=O)c1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "acetophenone"
    assert "octan" not in en
    assert "one" not in en or en == "acetophenone"
