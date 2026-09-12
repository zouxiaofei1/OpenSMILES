# IUPAC: P-25.3.3.3 / 表2.7
# Layer: e2e
"""pyrene 保留名 + 推荐编号(PIN): 纯烃与甲基衍生物。"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    ("c1cc2ccc3cccc4ccc(c1)c2c34", "pyrene", "芘"),
    ("Cc1ccc2ccc3cccc4ccc1c2c34", "1-methylpyrene", "1-甲基芘"),
    ("Cc1cc2ccc3cccc4ccc(c1)c2c34", "2-methylpyrene", "2-甲基芘"),
    ("Cc1cc2cccc3ccc4cccc1c4c32", "4-methylpyrene", "4-甲基芘"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyrene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
