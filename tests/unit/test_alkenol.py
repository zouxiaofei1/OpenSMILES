# IUPAC: P-63.1.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monoalcohols (alkenols).

Alcohol is the principal characteristic group; one non-aromatic C=C is inserted
as -ene_locant-en-OH_locant-ol / 首字-ene-烯-OH-醇. Parent chain covers the OH
carbon and both double-bond carbons. Numbering gives lowest OH locant first,
then lowest ene locant. No (E)/(Z) stereodescriptors this cycle.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenols (no E/Z)
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
    ("C=CCO", "prop-2-en-1-ol", "丙-2-烯-1-醇"),
    ("CC=CCO", "but-2-en-1-ol", "丁-2-烯-1-醇"),
    ("C=CCCCCO", "hex-5-en-1-ol", "己-5-烯-1-醇"),
    ("C=CC(C)O", "but-3-en-2-ol", "丁-3-烯-2-醇"),
    # negative: saturated alcohols, alkenoic acid, alkenal must not become alkenols
    ("CCO", "ethanol", "乙醇"),
    ("CC(O)C", "propan-2-ol", "丙-2-醇"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
