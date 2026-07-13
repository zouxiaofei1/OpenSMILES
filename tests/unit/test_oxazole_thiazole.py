# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Retained 1,3-oxazole / 1,3-thiazole parents (IUPAC P-22.2.1 / P-14.3.4).

Monocyclic fully aromatic five-membered rings with N+O or N+S in 1,3 relation
(O/S fixed as locant 1, N as 3). Simple ring substituents: halo and straight
n-alkyl C1–C3, total ≤2. Also relax mono-hetero5 (furan/thiophene/pyrrole)
side chains from methyl-only to n-alkyl C1–C2 (ethyl).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted 1,3 cores
    ("c1cocn1", "1,3-oxazole", "恶唑"),
    ("c1cscn1", "1,3-thiazole", "噻唑"),
    # positive: dialkyl oxazole (benchmark dual gold)
    ("C(CC)C=1N=COC1CCC", "4,5-dipropyl-1,3-oxazole", "4,5-二丙基恶唑"),
    # positive: monomethyl / monohalo on azole (C2 between O/S and N)
    ("Cc1ncco1", "2-methyl-1,3-oxazole", "2-甲基恶唑"),
    ("Clc1nccs1", "2-chloro-1,3-thiazole", "2-氯噻唑"),
    ("Cc1ncoc1", "4-methyl-1,3-oxazole", "4-甲基恶唑"),
    # positive: mono-hetero5 ethyl extension (eval_en only for gold CCc1ccco1)
    ("CCc1ccco1", "2-ethylfuran", None),
    ("CCc1cccs1", "2-ethylthiophene", "2-乙基噻吩"),
    # negative near-miss: mono-hetero5 / diazole / sat / fused must not regress
    ("c1ccoc1", "furan", "呋喃"),
    ("c1ccsc1", "thiophene", "噻吩"),
    ("c1cn[nH]c1", "1H-pyrazole", "吡唑"),
    ("c1ccc2[nH]ncc2c1", "1H-indazole", "1H-吲唑"),
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_oxazole_thiazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_oxazole_not_isoxazole_or_chain() -> None:
    """1,3-oxazole must not collapse to ether/alkane fallback."""
    r = SMILESNNamer().name("c1cocn1")
    assert r.success
    en = normalize_en(r.en)
    assert "oxazole" in en
    assert "methoxy" not in en
    assert en != "ethane"


def test_dipropyl_oxazole_locants() -> None:
    """O=1,N=3 fixed; two propyls at 4 and 5."""
    r = SMILESNNamer().name("C(CC)C=1N=COC1CCC")
    assert r.success
    en = normalize_en(r.en)
    assert en == "4,5-dipropyl-1,3-oxazole"
    assert "4,5-dipropyl1,3" not in en


def test_ethylfuran_not_hexane() -> None:
    """Mono-hetero5 n-ethyl side must remain furan parent, not chain."""
    r = SMILESNNamer().name("CCc1ccco1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-ethylfuran"
    assert "hexane" not in en
