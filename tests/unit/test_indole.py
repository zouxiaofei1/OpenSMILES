# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained parent 1H-indole (IUPAC P-22.2.1 / P-25):
benzo[b]pyrrole fused 6+5 aromatic system. Unsubstituted + at most one
ring monomethyl / monohalo (F/Cl/Br/I). NH = 1; standard locants 2–7.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    # positive: mono-methyl (pyrrole β = 3; α = 2)
    ("Cc1c[nH]c2ccccc12", "3-methyl-1H-indole", "3-甲基-1H-吲哚"),
    ("Cc1cc2ccccc2[nH]1", "2-methyl-1H-indole", "2-甲基-1H-吲哚"),
    # positive: mono-halo on benzene ring
    ("Brc1ccc2[nH]ccc2c1", "5-bromo-1H-indole", "5-溴-1H-吲哚"),
    ("Clc1ccc2[nH]ccc2c1", "5-chloro-1H-indole", "5-氯-1H-吲哚"),
    # negative: must not regress naphthalene / imidazole / benzene / pyridine
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("c1cnc[nH]1", "1H-imidazole", "咪唑"),
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_indole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2[nH]ccc2c1")
    assert r.success
    assert normalize_en(r.en) == "1h-indole"
    assert "octane" not in normalize_en(r.en)


def test_methyl_locant_is_3_not_2() -> None:
    """Cc1c[nH]c2ccccc12 is 3-methyl (β), not 2-methyl (α)."""
    r = SMILESNNamer().name("Cc1c[nH]c2ccccc12")
    assert r.success
    en = normalize_en(r.en)
    assert en == "3-methyl-1h-indole"
    assert "2-methyl" not in en
