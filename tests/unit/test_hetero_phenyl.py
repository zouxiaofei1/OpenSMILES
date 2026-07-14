# IUPAC: P-22.2.1 / P-29.3 / P-14.3.4
# Layer: L2,L3,L5
"""Phenyl on diazine / five-membered heteroarenes (plan 4.x)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    (
        "ClC1=NC=C(C(=N1)Cl)C1=CC=CC=C1",
        "2,4-dichloro-5-phenylpyrimidine",
        "2,4-二氯-5-苯基嘧啶",
    ),
    ("c1ccc(-c2ncccn2)cc1", "2-phenylpyrimidine", "2-苯基嘧啶"),
    ("C1(=CC=CC=C1)C=1SC=CC1", "2-phenylthiophene", "2-苯基噻吩"),
    ("c1ccc(-c2cccs2)cc1", "2-phenylthiophene", "2-苯基噻吩"),
    ("c1ccc(-c2ccco2)cc1", "2-phenylfuran", "2-苯基呋喃"),
    # regressions
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine", None),
    ("c1ccncc1", "pyridine", None),
    ("c1ccoc1", "furan", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hetero_phenyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
