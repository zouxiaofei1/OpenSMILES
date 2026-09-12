# IUPAC: P-31.1 / P-65.1.1
# Layer: L2, L4, L5
"""Polyunsaturated acyclic mono-FG parents (dienoic / trienoic / …).

Before: only mono-ene+FG was detected; 2+ C=C + mono FG fell back to saturated
naming. L4/L5 already had ene_locants + polyalkenoic + E/Z; gap was L2.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# (smiles, expected_en, expected_zh_or_None)
CASES = [
    # dienoic acids
    ("C=CC=CC(=O)O", "penta-2,4-dienoic acid", "戊-2,4-二烯酸"),
    ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid", "己-2,4-二烯酸"),
    ("C=CC=CCC(=O)O", "hexa-3,5-dienoic acid", "己-3,5-二烯酸"),
    # mono-ene acid still works
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    # saturated acid not stolen
    ("CCCC(=O)O", "butanoic acid", "丁酸"),
    # alcohol polyene already had path; keep regression
    ("C=CC=CCO", "penta-2,4-dien-1-ol", "戊-2,4-二烯-1-醇"),
    # near-miss: do not claim ring C=C as open polyenoic
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyenoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_polyenoic_not_saturated_collapse() -> None:
    r = SMILESNNamer().name("C=CC=CC(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert "dienoic" in en
    assert en != "pentanoic acid"
