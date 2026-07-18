# Universal claimable-block class-level red bars.
# Goal: every selected parent accounts for every heavy atom via typed, ordered
# side-block naming (or fail/retry another parent). Task 1 establishes RED only.
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("O=C(NC1CCCCCC1)c1ccccc1", "N-cycloheptylbenzamide", "N-环庚基苯甲酰胺"),
    ("O=C(Nc1cc(C)cc(C)c1)c1ccccc1", "N-(3,5-dimethylphenyl)benzamide", "N-(3,5-二甲基苯基)苯甲酰胺"),
    ("COC(C)CC", "2-methoxybutane", "2-甲氧基丁烷"),
    ("CCCCC(CC(C)CC)C(=O)O", "2-(2-methylbutyl)hexanoic acid", "2-(2-甲基丁基)己酸"),
    ("c1ccc(cc1)CC(C)CC", "(2-methylbutyl)benzene", "(2-甲基丁基)苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_universal_claimable_classes(smiles: str, en: str, zh: str) -> None:
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)
