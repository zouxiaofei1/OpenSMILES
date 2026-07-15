# IUPAC: P-63.1.1 / P-44
# Layer: L2
"""Reject methanol false parent for ring-carbon alcohols.

OH on a ring carbon with no open-chain arm must not collapse to C1 methanol.
Benzyl alcohol (exocyclic CH2OH) and true methanol must remain correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# Complex ring alcohols that collapsed to methanol
POS_NO_METHANOL = [
    "Nc1ccc2c(c1)CC(O)C2",  # aminoindanol-like
    "CC1CC(O)CC(C)(C)C1",  # substituted cyclohexanol not simple
]

NEG_CASES = [
    ("CO", "methanol", "甲醇"),
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)O", "propan-2-ol", "丙-2-醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("c1ccc(CO)cc1", "phenylmethanol", "苯基甲醇"),
]


@pytest.mark.parametrize("smiles", POS_NO_METHANOL)
def test_ring_alcohol_not_methanol(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methanol"
    assert not en.endswith("methanol") or "phenyl" in en


@pytest.mark.parametrize("smiles,en,zh", NEG_CASES)
def test_alcohol_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
