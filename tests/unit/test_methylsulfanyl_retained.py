# IUPAC: P-63.1.1 (substitutive nomenclature — retained prefix methylsulfanyl)
# Layer: L3
"""Retained methylsulfanyl (CH3S-) substituent leaf in RetainedBackend."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
# zh=None skips Chinese assertion.
CASES = [
    # Positive: simple CH3S- on primary amine chain
    ("CSCCCN", "3-methylsulfanylpropan-1-amine", None),
    # Positive: CH3S- on ethanol (C1 有可取代 H，2- 不可省略；匹配 methoxy 命名模式)
    ("CSCCO", "2-methylsulfanylethanol", None),
    # Positive: CH3S- on acetic acid
    ("CSCC(=O)O", "2-methylsulfanylacetic acid", None),
    # Negative: methoxy analogue — must remain unchanged
    ("COCCCN", "3-methoxypropan-1-amine", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_methylsulfanyl_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
