# IUPAC: P-65.1.1 / P-31.1 / P-93
# Layer: L2,L4,L5
"""Open-chain symmetric dialkyl alkenedioates (unsaturated diesters).

Mono-ene diesters of open-chain diacids: di{alkyl} (E/Z)-alk-n-enedioate /
(E/Z)-{烷}-n-烯二酸二{烷}酯. Saturated diesters must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: open-chain mono-ene symmetric diesters (E/Z + alkyl)
    (r"COC(=O)/C=C/C(=O)OC", "dimethyl (E)-but-2-enedioate", "(E)-丁-2-烯二酸二甲酯"),
    (r"COC(=O)/C=C\C(=O)OC", "dimethyl (Z)-but-2-enedioate", "(Z)-丁-2-烯二酸二甲酯"),
    (r"CCOC(=O)/C=C/C(=O)OCC", "diethyl (E)-but-2-enedioate", "(E)-丁-2-烯二酸二乙酯"),
    # positive: asymmetric mono-ene → lowest ene locant (P-31.1)
    (r"CCOC(=O)C/C=C/C(=O)OCC", "diethyl (E)-pent-2-enedioate", "(E)-戊-2-烯二酸二乙酯"),
    (r"COC(=O)CC/C=C/C(=O)OC", "dimethyl (E)-hex-2-enedioate", "(E)-己-2-烯二酸二甲酯"),
    # positive: saturated diesters must not regress
    ("CCOC(=O)C(=O)OCC", "diethyl oxalate", "草酸二乙酯"),
    ("CCOC(=O)CC(=O)OCC", "diethyl propanedioate", "丙二酸二乙酯"),
    # negative: mono ester, saturated diacid, anhydride
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("CC(=O)OC(=O)C", "acetic anhydride", "乙酸酐"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenedioate_diester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
