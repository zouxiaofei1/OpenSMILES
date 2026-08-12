# Universal claimable-block class-level red bars.
# Goal: every selected parent accounts for every heavy atom via typed, ordered
# side-block naming (or fail/retry another parent). Task 1 establishes RED only.
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_universal_claimable_classes(smiles: str, en: str, zh: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)
