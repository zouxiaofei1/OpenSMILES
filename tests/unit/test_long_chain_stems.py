# IUPAC: P-14.2.1
# Layer: L5
"""Long-chain open parent stems C11–C35 (numerical multipliers + alkane stems).

Extends alkane / alcohol / acid / aldehyde / haloalkane naming past C10 so
L5 no longer returns empty for unsupported n_carbons>10. Chinese multi-char
stems (十一…三十五) must not be truncated via zh[0].
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: C11–C35 parents (alkane / aldehyde / halo / hydroxyacid / alkene / alcohol / acid)
    ("CCCCCCCCCCCCCCCCCCCCCCCCCC", "hexacosane", "二十六烷"),
    ("C(CCCCCCCCCC)=O", "undecanal", "十一醛"),
    ("CCCCCCCCCCCCCCI", "1-iodotetradecane", "1-碘十四烷"),
    ("CCCCCCCCCCC(O)C(=O)O", "2-hydroxydodecanoic acid", "2-羟基十二酸"),
    ("C=CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "pentatriacont-1-ene", "三十五-1-烯"),
    ("CCCCCCCCCCC", "undecane", "十一烷"),
    ("CCCCCCCCCCCCO", "dodecan-1-ol", "十二-1-醇"),
    ("CCCCCCCCCCC(=O)O", "undecanoic acid", "十一酸"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC=O", "triacontanal", "三十醛"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO", "tetratriacontan-1-ol", "三十四-1-醇"),
    # negative: C≤10 retained / boundary / known alkenol-as-alcohol must not regress
    ("CCO", "ethanol", "乙醇"),
    ("CCCCCCCCCC", "decane", "癸烷"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("C=CCCO", "butan-1-ol", "丁-1-醇"),
    ("CCCCCCCCCO", "nonan-1-ol", "壬-1-醇"),
    ("C=O", "formaldehyde", "甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_long_chain_stems(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
