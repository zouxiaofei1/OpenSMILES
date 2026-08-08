# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent 1H-indazole (IUPAC P-22.2.1 / P-25):
benzo[c]pyrazole / 1,2-diazaindene fused 6+5 aromatic system.
Unsubstituted + ≤2 ring monomethyl / monohalo. NH=1, N=2.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted
    ("c1ccc2[nH]ncc2c1", "1H-indazole", "1H-吲唑"),
    # positive: mono-methyl (position 3)
    ("Cc1n[nH]c2ccccc12", "3-methyl-1H-indazole", "3-甲基-1H-吲唑"),
    # positive: mono-halo
    ("Clc1ccc2[nH]ncc2c1", "5-chloro-1H-indazole", "5-氯-1H-吲唑"),
    # positive: diiodo (benchmark dual)
    ("IC1=NNC2=CC(I)=CC=C12", "3,6-diiodo-1H-indazole", "3,6-二碘-1H-吲唑"),
    # positive: mono carbonitrile (benchmark dual)
    ("N1N=CC2=CC=C(C=C12)C#N", "1H-indazole-6-carbonitrile", "1H-吲唑-6-甲腈"),
    # positive: mono carbaldehyde (benchmark dual)
    ("N1N=CC2=CC(=CC=C12)C=O", "1H-indazole-5-carbaldehyde", "1H-吲唑-5-甲醛"),
    # negative: must not regress indole / benzofuran / pyrazole / quinoline
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("c1ccc2occc2c1", "benzofuran", "苯并呋喃"),
    ("c1cn[nH]c1", "1H-pyrazole", "吡唑"),
    ("c1ccc2ncccc2c1", "quinoline", "喹啉"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2[nH]ncc2c1")
    assert r.success
    assert normalize_en(r.en) == "1h-indazole"
    assert "octane" not in normalize_en(r.en)


def test_diiodo_locants_3_6() -> None:
    """IC1=NNC2=CC(I)=CC=C12 is 3,6-diiodo (not 3,5)."""
    r = SMILESNNamer().name("IC1=NNC2=CC(I)=CC=C12")
    assert r.success
    en = normalize_en(r.en)
    assert en == "3,6-diiodo-1h-indazole"
    assert "3,5" not in en
