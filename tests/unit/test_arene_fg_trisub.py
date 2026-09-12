# IUPAC: P-14.3.4
# Layer: L2
"""Arene FG retained parents: raise simple ring-sub cap 2→3.

IUPAC P-14.3.4 / P-63.1.4 / P-65.1.1.1 — phenol, aniline, benzoic acid
already keep 0–2 simple ring substituents (halo/methyl/nitro/alkoxy);
align with benzene's ≤3 simple ring-sub allowance so tri-substituted
retained parents stay on the FG parent (not chain / fail).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: tri-substituted retained FG parents
    ("CC1=C(C=C(C(=C1)C)C)O", "2,4,5-trimethylphenol", "2,4,5-三甲基苯酚"),
    ("ClC=1C(=C(N)C=CC1Cl)C", "3,4-dichloro-2-methylaniline", "3,4-二氯-2-甲基苯胺"),
    (
        "BrC=1C(=CC(=C(C(=O)O)C1)F)Cl",
        "5-bromo-4-chloro-2-fluorobenzoic acid",
        "5-溴-4-氯-2-氟苯甲酸",
    ),
    ("Oc1c(Br)cc(Br)cc1Br", "2,4,6-tribromophenol", None),
    # negative: must not regress existing correct names
    ("Cc1ccc(O)c(C)c1", "2,4-dimethylphenol", "2,4-二甲基苯酚"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("c1ccccc1", "benzene", "苯"),
    ("Clc1ccc(Cl)c(Cl)c1", "1,2,4-trichlorobenzene", "1,2,4-三氯苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_fg_trisub(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
