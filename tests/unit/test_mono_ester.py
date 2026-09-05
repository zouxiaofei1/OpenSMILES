# IUPAC: P-65.6
# Layer: L1,L2,L4,L5
"""Simple acyclic saturated monoesters (alkyl alkanoates).

Ester R–C(=O)–O–R': carbonyl C has =O and single-bond O with no H
(alkoxy oxygen, not OH). Functional class: alkyl alkanoate / 酸+烷词干+酯.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: simple saturated monoesters
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("COC(C)=O", "methyl acetate", "乙酸甲酯"),
    ("CCOC(=O)CC", "ethyl propanoate", "丙酸乙酯"),
    ("CCCC(=O)OC", "methyl butanoate", "丁酸甲酯"),
    ("CC(=O)OCCC", "propyl acetate", "乙酸丙酯"),
    ("CCOC(=O)CCC", "ethyl butanoate", "丁酸乙酯"),
    ("CCCCCC(=O)OC", "methyl hexanoate", "己酸甲酯"),
    # negative: acid / aldehyde / alkane / alcohol / ketone / alkyne / alkene
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCC", "propane", "丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_ester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
