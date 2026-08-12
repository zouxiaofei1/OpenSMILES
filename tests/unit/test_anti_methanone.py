# IUPAC: P-64.2.1 / P-44
# Layer: L2
"""Reject methanone false parent for ring-only ketone carbonyls.

When the ketone carbon has no open-chain carbon arms (only ring/aryl C
neighbors), open-chain kind=ketone must not collapse to C1 methanone.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

POS_NO_METHANONE = [
    "O[C@]1(C(=CC(C1)=O)C1=CC=CC=C1)C1=CC=CC=C1",
    "C(CCC)C=1OC2=C(C(C1)=O)C=C(C=C2)Cl",
]

NEG_CASES = [
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
    ("CCC(=O)CC", "pentan-3-one", "戊-3-酮"),
    ("CC(=O)C(C)C", "3-methylbutan-2-one", None),
]


@pytest.mark.parametrize("smiles", POS_NO_METHANONE)
def test_ring_ketone_not_methanone(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en != "methanone"
    assert "methanone" not in en


@pytest.mark.parametrize("smiles,en,zh", NEG_CASES)
def test_open_ketone_regressions(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
