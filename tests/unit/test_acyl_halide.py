# IUPAC: P-65.5
# Layer: L1,L2,L5
"""Acyl halides: alkanoyl chloride/bromide + benzoyl chloride/bromide (P-65.5).

Open-chain mono acyl bromides C2–C10 aligned with acyl chlorides; retained
isobutyryl bromide per benchmark gold; benzoyl bromide for unsubstituted
Ph–C(=O)Br. Negatives: acid, bromoalkane, aldehyde.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyl chloride regression
    ("CC(=O)Cl", "acetyl chloride", "乙酰氯"),
    ("ClC(=O)c1ccccc1", "benzoyl chloride", "苯甲酰氯"),
    ("c1ccc(C(=O)Cl)cc1", "benzoyl chloride", "苯甲酰氯"),
    # positive: open-chain mono acyl bromides
    ("CC(=O)Br", "acetyl bromide", "乙酰溴"),
    ("CCC(=O)Br", "propanoyl bromide", "丙酰溴"),
    ("CCCC(=O)Br", "butanoyl bromide", "丁酰溴"),
    # retained isobutyryl bromide (benchmark gold)
    ("C(C(C)C)(=O)Br", "isobutyryl bromide", "异丁酰溴"),
    ("CC(C)C(=O)Br", "isobutyryl bromide", "异丁酰溴"),
    # positive: benzoyl bromide
    ("c1ccc(C(=O)Br)cc1", "benzoyl bromide", "苯甲酰溴"),
    ("BrC(=O)c1ccccc1", "benzoyl bromide", "苯甲酰溴"),
    # negative: must not become acyl bromide
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCBr", "bromoethane", "溴乙烷"),
    ("CC=O", "acetaldehyde", "乙醛"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_acyl_halide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
