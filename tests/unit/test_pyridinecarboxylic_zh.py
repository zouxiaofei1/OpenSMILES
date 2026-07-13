# IUPAC: P-65.1.1
# Layer: L5
"""Chinese retained suffix for pyridinecarboxylic acids: 甲酸 (not 羧酸).

IUPAC P-65.1.1 retained pyridinecarboxylic acids; Chinese gold uses
吡啶-n-甲酸, aligned with benzoic→苯甲酸 and cycloalkanecarboxylic→环…甲酸.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: pyridine-2/3/4-carboxylic acid → 吡啶-n-甲酸
    ("OC(=O)c1ccccn1", "pyridine-2-carboxylic acid", "吡啶-2-甲酸"),
    ("OC(=O)c1cccnc1", "pyridine-3-carboxylic acid", "吡啶-3-甲酸"),
    ("OC(=O)c1ccncc1", "pyridine-4-carboxylic acid", "吡啶-4-甲酸"),
    # negative: carbocyclic carboxylic acids keep correct 甲酸 (no regression)
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinecarboxylic_zh(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
