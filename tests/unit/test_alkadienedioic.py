# IUPAC: P-65.1.1 / P-31.1 / P-72.2.2.1
# Layer: L2,L4,L5
"""Open-chain polyunsaturated dicarboxylic acids (alkadienedioic).

Exactly two carboxyls (acid or anion) + ≥2 open-chain C=C; parent covers both
carboxyl carbons and all double-bond carbons. English:
(2E,4E)-hexa-2,4-dienedioic acid / …dienedioate; Chinese: …-2,4-二烯二酸 / 根.
Reuses multi E/Z prefix and diacid anion conversion.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: open-chain dienedioic acids / anions with multi E/Z
    (
        r"O=C([O-])/C=C/C=C/C(=O)[O-]",
        "(2E,4E)-hexa-2,4-dienedioate",
        "(2E,4E)-己-2,4-二烯二酸根",
    ),
    (
        r"O=C([O-])/C=C\C=C/C(=O)[O-]",
        "(2Z,4Z)-hexa-2,4-dienedioate",
        "(2Z,4Z)-己-2,4-二烯二酸根",
    ),
    (
        r"O=C(O)/C=C/C=C/C(=O)O",
        "(2E,4E)-hexa-2,4-dienedioic acid",
        "(2E,4E)-己-2,4-二烯二酸",
    ),
    # negative: saturated diacid anion and mono-ene diacid must not regress
    (r"O=C([O-])CCC(=O)[O-]", "butanedioate", "丁二酸根"),
    (r"O=C(O)/C=C/C(=O)O", "(2E)-but-2-enedioic acid", "(2E)-丁-2-烯二酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkadienedioic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
