# IUPAC: P-61.5.1
# Layer: L1,L2,L3,L5
"""Nitro as prefix (nitrobenzene / nitrophenol): simple arene cases.

Parent = benzene or phenol; nitro is a non-senior prefix substituent.
Allow 1–2 nitros + optional halo/methyl (total subs ≤3 for benzene;
phenol: OH principal + nitro prefix).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: mono-nitro benzene
    ("[O-][N+](=O)c1ccccc1", "nitrobenzene", "硝基苯"),
    # positive: nitro + halo (lowest set of locants)
    ("O=[N+]([O-])c1cccc(I)c1", "1-iodo-3-nitrobenzene", None),
    # positive: 4-nitrophenol (OH principal, nitro prefix)
    ("Oc1ccc([N+]([O-])=O)cc1", "4-nitrophenol", "4-硝基苯酚"),
    # negative: bare / monohalo benzene, phenol, chain alcohol stay correct
    ("c1ccccc1", "benzene", "苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_nitrobenzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_nitromethane_not_nitrobenzene() -> None:
    """Aliphatic nitro must not be named as nitrobenzene."""
    r = SMILESNNamer().name("C[N+](=O)[O-]")
    en = normalize_en(r.en) if r.success else ""
    assert en != "nitrobenzene"
    assert "benzene" not in en
    assert "苯" not in (r.zh or "")
