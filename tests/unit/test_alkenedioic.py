# IUPAC: P-44 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated dicarboxylic acids (alkenedioic acids).

Exactly two COOH + one C=C; parent chain through both carboxyl carbons and the
double bond. English: (E/Z)-alk-n-enedioic acid; Chinese: (E/Z)-{烷首}-n-烯二酸.
Stereo from RDKit BondStereo; no maleic/fumaric retained names.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: open-chain alkenedioic acids with E/Z
    (r"OC(=O)/C=C/C(=O)O", "(E)-but-2-enedioic acid", "(E)-丁-2-烯二酸"),
    (r"OC(=O)/C=C\C(=O)O", "(Z)-but-2-enedioic acid", "(Z)-丁-2-烯二酸"),
    (r"O=C(O)/C=C/CCC(=O)O", "(E)-hex-2-enedioic acid", None),
    # negative: saturated diacids, mono alkenoic acid, non-acids
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("OC(=O)CCC(=O)O", "butanedioic acid", "丁二酸"),
    ("OC(=O)C=C", "prop-2-enoic acid", "丙-2-烯酸"),
    ("CCO", "ethanol", "乙醇"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenedioic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
