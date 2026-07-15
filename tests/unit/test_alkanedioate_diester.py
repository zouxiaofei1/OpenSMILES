# IUPAC: P-65.1.1 / P-65.6
# Layer: L2,L5
"""Open-chain symmetric dialkyl alkanedioates (diesters of diacids).

Symmetric saturated diesters: di{alkyl} {alkane}dioate / {二酸}二{烷}酯.
C2 retained oxalate/草酸; C3+ systematic propanedioate etc.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: symmetric open-chain saturated diesters
    ("CCOC(=O)C(=O)OCC", "diethyl oxalate", "草酸二乙酯"),
    ("COC(=O)C(=O)OC", "dimethyl oxalate", "草酸二甲酯"),
    ("CCOC(=O)CC(=O)OCC", "diethyl propanedioate", "丙二酸二乙酯"),
    ("CCOC(=O)CCC(=O)OCC", "diethyl butanedioate", "丁二酸二乙酯"),
    ("COC(=O)CCCC(=O)OC", "dimethyl pentanedioate", "戊二酸二甲酯"),
    ("CCOC(=O)CCCCC(=O)OCC", "diethyl hexanedioate", "己二酸二乙酯"),
    # negative: mono ester, diacid, anhydride (must not become diester)
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("O=C(O)C(=O)O", "oxalic acid", "草酸"),
    ("CC(=O)OC(=O)C", "acetic anhydride", "乙酸酐"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanedioate_diester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
