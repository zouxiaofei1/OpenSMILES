# IUPAC: P-65.1.1 / P-31.1 / P-93
# Layer: L2,L4,L5
"""Open-chain symmetric dialkyl alkenedioates (unsaturated diesters).

Mono-ene diesters of open-chain diacids: di{alkyl} (E/Z)-alk-n-enedioate /
(E/Z)-{烷}-n-烯二酸二{烷}酯. Saturated diesters must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: mono ester, saturated diacid, anhydride
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenedioate_diester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
