# IUPAC: P-72.2.2.1 / P-65.1.1
# Layer: L0,L5
"""Alkali metal carboxylates: salt dissociation + functional class salt names.

L0 dissociates single alkali metal cation + one organic carboxylate fragment.
L5 prefixes metal (EN) or replaces 酸根 → 酸钠/钾/锂 (ZH).
Bare anions and neutral acids must not regress.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: alkali metal + carboxylate
    ("[Na+].[O-]C(=O)C", "sodium acetate", "乙酸钠"),
    ("[Na+].[O-]C(=O)c1ccccc1", "sodium benzoate", "苯甲酸钠"),
    ("FC1=C(C(=O)[O-])C=CC(=C1)F.[Na+]", "sodium 2,4-difluorobenzoate", "2,4-二氟苯甲酸钠"),
    ("O=C([O-])C.[K+]", "potassium acetate", "乙酸钾"),
    ("O=C[O-].[Na+]", "sodium formate", "甲酸钠"),
    ("[Li+].[O-]C(=O)C", "lithium acetate", "乙酸锂"),
    # positive: metal after R/S (and optional E/Z) on anion stem
    ("C[C@H](O)C(=O)[O-].[Na+]", "sodium (2S)-2-hydroxypropanoate", "(2S)-2-羟基丙酸钠"),
    ("C/C=C/[C@H](O)C(=O)[O-].[Na+]", "sodium (2S,3E)-2-hydroxypent-3-enoate", "(2S,3E)-2-羟基戊-3-烯酸钠"),
    # negative: bare anion / neutral acid / ester must not regress
    ("O=C([O-])C", "acetate", "乙酸根"),
    ("C[C@H](O)C(=O)[O-]", "(2S)-2-hydroxypropanoate", "(2S)-2-羟基丙酸根"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_metal_carboxylate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
