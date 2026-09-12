# IUPAC: P-62.2.1
# Layer: L1,L2,L4,L5
"""Simple acyclic primary monoamines (alkanamines).

Primary amine: N with exactly one carbon neighbor and ≥2 H; not amide nitrogen.
Parent chain through the carbon attached to N; kind=amine.
C1/C2 omit locant (methanamine/ethanamine); C≥3 use …-n-amine / …-n-胺.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: primary monoamines (straight + optional branched)
    ("CN", "methanamine", "甲胺"),
    ("CCN", "ethanamine", "乙胺"),
    ("CCCN", "propan-1-amine", "丙-1-胺"),
    ("CCCCN", "butan-1-amine", "丁-1-胺"),
    ("CCCCCN", "pentan-1-amine", "戊-1-胺"),
    ("CC(C)N", "propan-2-amine", "丙-2-胺"),
    # negative: alkane / alcohol / aldehyde / acid / ketone / ester / unsat
    ("CCC", "propane", "丙烷"),
    ("CCO", "ethanol", "乙醇"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC(C)=O", "propan-2-one", "丙-2-酮"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_amine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_amide_not_named_as_amine() -> None:
    """Amide nitrogen (N next to carbonyl) must not become *amine."""
    r = SMILESNNamer().name("CC(=O)N")
    assert r.success
    assert "amine" not in normalize_en(r.en)
