# IUPAC: P-22.2.1
# Layer: L2,L4,L5
"""Retained monocyclic heteroarene parents (IUPAC P-22.2.1):
- 5-membered mono-hetero: furan / thiophene / 1H-pyrrole (unsub, monomethyl, monohalo)
- 6-membered diazines: pyrimidine / pyrazine / pyridazine (unsubstituted only)
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted 5-membered cores
    ("c1ccoc1", "furan", "呋喃"),
    ("c1ccsc1", "thiophene", "噻吩"),
    ("c1cc[nH]c1", "1H-pyrrole", "吡咯"),
    # positive: monomethyl / monohalo (hetero = 1 → lowest locant)
    ("Cc1ccco1", "2-methylfuran", "2-甲基呋喃"),
    ("Cc1cccs1", "2-methylthiophene", "2-甲基噻吩"),
    ("Clc1ccco1", "2-chlorofuran", "2-氯呋喃"),
    # positive: unsubstituted diazines (1,3 / 1,4 / 1,2)
    ("c1cncnc1", "pyrimidine", "嘧啶"),
    ("c1cnccn1", "pyrazine", "吡嗪"),
    ("c1ccnnc1", "pyridazine", "哒嗪"),
    # negative: must not regress carbocycles / pyridine / open chain
    ("c1ccccc1", "benzene", "苯"),
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Cc1ccccn1", "2-methylpyridine", "2-甲基吡啶"),
    ("CCCCC", "pentane", "戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_heteroarene5(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methylfuran_locant_is_2_not_5() -> None:
    """O fixed as 1; monomethyl must be 2- not 5-."""
    r = SMILESNNamer().name("Cc1ccco1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylfuran"
    assert "5-methyl" not in en


def test_pyrrole_not_alkane() -> None:
    r = SMILESNNamer().name("c1cc[nH]c1")
    assert r.success
    assert normalize_en(r.en) != "butane"
    assert "pyrrole" in normalize_en(r.en)
