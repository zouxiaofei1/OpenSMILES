# IUPAC: P-66.1.1.1.3 / P-29.3
# Layer: L2,L3
"""Amide N-benzyl (and substituted benzyl) claim — not N-methyl.

Open-chain monoamide parent; N–CH2–Ph (Ph may carry simple leaves) is
N-benzyl / N-(…benzyl). N-phenyl and plain N-alkyl stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # N-benzyl formamide / acetamide
    ("O=CNCc1ccccc1", "N-benzylformamide", "N-苄基甲酰胺"),
    ("CC(=O)NCc1ccccc1", "N-benzylacetamide", "N-苄基乙酰胺"),
    # substituted benzyl
    (
        "O=CNCc1ccc(Cl)cc1",
        "N-(4-chlorobenzyl)formamide",
        "N-(4-氯苄基)甲酰胺",
    ),
    (
        "CC(=O)NCc1ccc(OC)cc1",
        "N-(4-methoxybenzyl)acetamide",
        "N-(4-甲氧基苄基)乙酰胺",
    ),
    (
        "O=CNCc1ccc(C)cc1",
        "N-(4-methylbenzyl)formamide",
        "N-(4-甲基苄基)甲酰胺",
    ),
    # negatives: N-phenyl, plain N-alkyl, primary amide
    ("O=CNc1ccccc1", "N-phenylformamide", "N-苯基甲酰胺"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("O=CNCCCC", "N-butylformamide", "N-丁基甲酰胺"),
    ("CC(=O)NCC", "N-ethylacetamide", "N-乙基乙酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_amide_n_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
