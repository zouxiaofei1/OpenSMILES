# IUPAC: P-65.1.1.1
# Layer: L2,L3,L5
"""Arene retained acid/aldehyde: ring amino/alkoxy/nitro as prefixes.

Benzoic acid (and benzaldehyde) keep COOH/CHO as principal FG; simple
ring amino, methoxy/ethoxy, and nitro are allowed prefixes alongside
existing halo/methyl/hydroxy (sub cap ≤3). Aligns with phenol/benzene.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: aminobenzoic
    ("Nc1ccccc1C(=O)O", "2-aminobenzoic acid", None),
    ("Nc1ccc(C(=O)O)cc1", "4-aminobenzoic acid", "4-氨基苯甲酸"),
    # positive: amino + methoxy + dimethyl
    ("NC1=CC(=C(C(=O)O)C=C1)OC", "4-amino-2-methoxybenzoic acid", "4-氨基-2-甲氧基苯甲酸"),
    ("NC1=C(C=C(C(=O)O)C=C1C)C", "4-amino-3,5-dimethylbenzoic acid", None),
    # positive: methoxy / nitro only
    ("COc1ccc(C(=O)O)cc1", "4-methoxybenzoic acid", "4-甲氧基苯甲酸"),
    ("O=[N+]([O-])c1ccc(C(=O)O)cc1", "4-nitrobenzoic acid", "4-硝基苯甲酸"),
    # optional benzaldehyde
    ("O=Cc1ccc(OC)cc1", "4-methoxybenzaldehyde", "4-甲氧基苯甲醛"),
    # negative near-miss: must not regress
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("OC(=O)c1ccc(O)cc1", "4-hydroxybenzoic acid", "4-羟基苯甲酸"),
    ("OC(=O)c1ccc(Cl)cc1", "4-chlorobenzoic acid", "4-氯苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_aminobenzoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_aminobenzoic_not_chain_acid() -> None:
    """Ring carboxyl + amino must not expand to aminoalkanoic acid."""
    r = SMILESNNamer().name("Nc1ccccc1C(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-aminobenzoic acid"
    assert "heptanoic" not in en
    assert "hexanoic" not in en
