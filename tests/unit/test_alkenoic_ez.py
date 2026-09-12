# IUPAC: P-93.4 / P-31.1
# Layer: L2,L5
"""Open-chain monounsaturated monocarboxylic acids with (E)/(Z) stereo.

When RDKit BondStereo is STEREOE/STEREOZ on the parent C=C, prefix (E)-/(Z)-.
No stereo marker → no descriptor (prop-2-enoic acid stays plain).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: mono-alkenoic acids with E/Z from BondStereo
    (r"CCCCCC/C=C/C(=O)O", "(2E)-non-2-enoic acid", None),
    (r"CC/C=C/CC(=O)O", "(3E)-hex-3-enoic acid", None),
    (r"CCCCC/C=C\CCC(=O)O", "(4Z)-dec-4-enoic acid", None),
    # negative: no forced stereo, diacid, non-acid
    ("OC(=O)C=C", "prop-2-enoic acid", "丙-2-烯酸"),
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
    (r"OC(=O)/C=C/C(=O)O", "(2E)-but-2-enedioic acid", "(2E)-丁-2-烯二酸"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenoic_ez(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
