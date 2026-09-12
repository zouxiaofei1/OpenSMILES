# IUPAC: P-66.1.1.1.3 / P-29.3
# Layer: L2,L3
"""Amide N-benzyl (and substituted benzyl) claim — not N-methyl.

Open-chain monoamide parent; N–CH2–Ph (Ph may carry simple leaves) is
N-benzyl / N-(…benzyl). N-phenyl and plain N-alkyl stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # negatives: N-phenyl, plain N-alkyl, primary amide
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_amide_n_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
