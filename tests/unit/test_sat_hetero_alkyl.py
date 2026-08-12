# IUPAC: P-22.2.2 / P-14.3.4
# Layer: L2,L3,L4
"""Sat-hetero ring-alkyl scope (P-22.2.2) — negative guard only.

Positive ring-alkyl sat-hetero cases are not covered in this file. The retained
case asserts acetamide is not named as a sat-hetero parent.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: open-chain amide, not sat-hetero
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_alkyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_n_methyl_not_c_methyl_piperidine() -> None:
    """N-methylpiperidine is out of scope; must not invent a C-methyl name."""
    r = SMILESNNamer().name("CN1CCCCC1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "piperidine" not in en
    assert en != "4-methylpiperidine"
    assert en != "2-methylpiperidine"
    assert en != "3-methylpiperidine"
