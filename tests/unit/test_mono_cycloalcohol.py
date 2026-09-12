# IUPAC: P-63.1.1 / P-22.1.1
# Layer: L2,L4,L5
"""Unsubstituted monocyclic monoalcohols (cycloalkanols).

Parent = saturated monocarbocycle with one ring-carbon OH.
Unsubstituted: omit locant (cyclohexanol not cyclohexan-1-ol).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted monocycloalkanols C3–C7
    ("OC1CC1", "cyclopropanol", "环丙醇"),
    ("OC1CCC1", "cyclobutanol", "环丁醇"),
    ("OC1CCCC1", "cyclopentanol", "环戊醇"),
    ("OC1CCCCC1", "cyclohexanol", "环己醇"),
    ("OC1CCCCCC1", "cycloheptanol", "环庚醇"),
    # negative: acyclic alcohol / unsubstituted ring / monoalkyl ring
    ("CCO", "ethanol", "乙醇"),
    ("CCCO", "propan-1-ol", "丙-1-醇"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("CC1CCCCC1", "methylcyclohexane", "甲基环己烷"),
    ("CCC", "propane", "丙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_mono_cycloalcohol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cyclohexanol_not_chain_alcohol() -> None:
    """Ring OH must not be named as acyclic alcohol (hexanol / nonan-*-ol)."""
    r = SMILESNNamer().name("OC1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "cyclohexanol"
    assert "hexanol" != en or en.startswith("cyclo")
    assert "nonan" not in en
