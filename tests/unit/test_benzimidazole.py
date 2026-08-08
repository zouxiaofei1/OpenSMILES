# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained parent 1H-benzimidazole (IUPAC P-22.2.1 / P-25 / P-62.2.1):
fused aromatic 6C+5(C3N2); NH=1, N=3 (imidazole-type). Unsub / ≤2
halo+methyl+CF3; 2-primary amine (amino or 2-imino dual-NH tautomer)
→ benzimidazolamine (with ≤1 halo or CF3).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted 1H-benzimidazole
    ("c1ccc2[nH]cnc2c1", "1H-benzimidazole", "1H-苯并咪唑"),
    ("c1nc2ccccc2[nH]1", "1H-benzimidazole", "1H-苯并咪唑"),
    # positive: mono-halo / mono-methyl on ring
    ("Clc1ccc2[nH]cnc2c1", "5-chloro-1H-benzimidazole", "5-氯-1H-苯并咪唑"),
    ("Cc1nc2ccccc2[nH]1", "2-methyl-1H-benzimidazole", "2-甲基-1H-苯并咪唑"),
    # positive: unsubstituted 2-amine (amino tautomer)
    ("Nc1nc2ccccc2[nH]1", "1H-benzimidazol-2-amine", "1H-苯并咪唑-2-胺"),
    # positive: dual gold — 2-imino dual-NH tautomer
    ("N=c1[nH]c2ccccc2[nH]1", "1H-benzimidazol-2-amine", "1H-苯并咪唑-2-胺"),
    # positive: 2-amine + mono-halo
    (
        "Nc1nc2ccc(Br)cc2[nH]1",
        "6-bromo-1H-benzimidazol-2-amine",
        "6-溴-1H-苯并咪唑-2-胺",
    ),
    # negative: near neighbors must not regress
    ("c1ccc2[nH]ncc2c1", "1H-indazole", "1H-吲唑"),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "吲哚"),
    ("c1nc2ccccc2s1", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("c1c[nH]cn1", "1H-imidazole", "咪唑"),
    ("c1ccccc1", "benzene", "苯"),
    ("Nc1ccccc1", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzimidazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_hexane() -> None:
    r = SMILESNNamer().name("c1ccc2[nH]cnc2c1")
    assert r.success
    assert normalize_en(r.en) == "1h-benzimidazole"
    assert "hexane" not in normalize_en(r.en)


def test_gold_amine_not_hexane() -> None:
    r = SMILESNNamer().name("N=c1[nH]c2ccccc2[nH]1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1h-benzimidazol-2-amine"
    assert "hexane" not in en
    assert "methanamine" not in en


def test_amino_tautomer_not_methanamine() -> None:
    r = SMILESNNamer().name("Nc1nc2ccccc2[nH]1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1h-benzimidazol-2-amine"
    assert "methanamine" not in en


def test_chloro_locant_is_5() -> None:
    r = SMILESNNamer().name("Clc1ccc2[nH]cnc2c1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "5-chloro-1h-benzimidazole"
    assert "6-chloro" not in en


def test_amine_locant_is_2() -> None:
    r = SMILESNNamer().name("Nc1nc2ccccc2[nH]1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1h-benzimidazol-2-amine"
    assert "3-amine" not in en
