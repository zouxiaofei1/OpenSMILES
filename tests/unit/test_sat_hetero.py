# IUPAC: P-22.2.2
# Layer: L2,L3,L4,L5
"""Saturated monohetero parent scope (P-22.2.2) — negative guard only.

Positive sat-hetero cases (oxolane, piperidine, …) are not covered in this
file. The retained cases assert carbocycle / arene / open chain must not be
named as Hantzsch–Widman parents.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: carbocycle / arene / open chain must not break
    ("C1CCCCC1", "cyclohexane", "环己烷"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("CCCCC", "pentane", "戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_n_methylpiperidine_not_parent() -> None:
    """N-alkyl on ring N is out of scope (negative: must not claim piperidine)."""
    r = SMILESNNamer().name("CN1CCCCC1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    assert "piperidine" not in normalize_en(r.en)


def test_piperidine_carboxylic_not_oxolane() -> None:
    """Principal acid FG: must not invent sat-hetero parent name."""
    r = SMILESNNamer().name("O=C(O)C1CCCCN1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "oxolane" not in en
    assert "oxane" not in en
