# IUPAC: P-25.4.1 / 表2.7
# Layer: e2e
"""phenanthrene 保留名 + 传统编号(PIN): 纯烃与甲基衍生物。"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("c1ccc2c(c1)ccc1ccccc12", "phenanthrene", "菲"),
    ("Cc1cccc2c1ccc1ccccc12", "1-methylphenanthrene", "1-甲基菲"),
    ("Cc1ccc2c(ccc3ccccc32)c1", "2-methylphenanthrene", "2-甲基菲"),
    ("Cc1ccc2ccc3ccccc3c2c1", "3-methylphenanthrene", "3-甲基菲"),
    ("Cc1cccc2ccc3ccccc3c12", "4-methylphenanthrene", "4-甲基菲"),
    ("Cc1cc2ccccc2c2ccccc12", "9-methylphenanthrene", "9-甲基菲"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_phenanthrene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
