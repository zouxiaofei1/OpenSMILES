# IUPAC: P-66.5.1 / P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained pyridinecarbonitrile parent (pyridine-n-carbonitrile / 吡啶-n-甲腈).

IUPAC P-66.5.1 nitriles; P-22.2.1 retained pyridine parent with principal
nitrile suffix. N=1; CN ring-attach locant; ≤3 simple ring prefixes
(halo / methyl / methoxy / nitro). Prefer system name over nicotinonitrile.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted pyridine-2/3/4-carbonitrile
    ("N#Cc1ccccn1", "pyridine-2-carbonitrile", "吡啶-2-甲腈"),
    ("N#Cc1cccnc1", "pyridine-3-carbonitrile", "吡啶-3-甲腈"),
    ("N#Cc1ccncc1", "pyridine-4-carbonitrile", "吡啶-4-甲腈"),
    # positive: dimethyl (benchmark dual)
    (
        "CC=1C(=NC(=CC1)C)C#N",
        "3,6-dimethylpyridine-2-carbonitrile",
        "3,6-二甲基吡啶-2-甲腈",
    ),
    # positive: methoxy (system name, not 4-methoxynicotinonitrile)
    (
        "COC1=CC=NC=C1C#N",
        "4-methoxypyridine-3-carbonitrile",
        "4-甲氧基吡啶-3-甲腈",
    ),
    # positive: halo (N=1; lowest set of CN + halo locants)
    ("N#Cc1ccc(Cl)nc1", "2-chloropyridine-5-carbonitrile", "2-氯吡啶-5-甲腈"),
    ("N#Cc1cc(Cl)ccn1", "4-chloropyridine-2-carbonitrile", "4-氯吡啶-2-甲腈"),
    # positive: methyl + halo
    ("Cc1ncccc1C#N", "2-methylpyridine-3-carbonitrile", "2-甲基吡啶-3-甲腈"),
    # negative: must not steal benzonitrile / pyridine / acid / chain nitrile
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("O=C(O)c1ccncc1", "pyridine-4-carboxylic acid", "吡啶-4-甲酸"),
    ("CC#N", "acetonitrile", "乙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinecarbonitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_nitrile() -> None:
    """Pyridine nitrile must not expand ring into alkanenitrile parent."""
    r = SMILESNNamer().name("N#Cc1ccccn1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "pyridine-2-carbonitrile"
    assert "hexanenitrile" not in en
    assert "pentanenitrile" not in en
