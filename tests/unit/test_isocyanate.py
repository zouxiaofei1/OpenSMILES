# IUPAC: P-61.9
# Layer: L1,L2,L3,L5
"""Isocyanate / isothiocyanate (R–N=C=X, X=O/S): functional class + arene prefix.

Alkyl mono: functional-class '{alkyl} isocyanate/isothiocyanate'.
Aryl: isocyanato/isothiocyanato prefixes on benzene (with halo/methyl).
Must not mis-detect as amide (N=C=O) or collapse isothiocyanate to alkane.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: alkyl isocyanates (functional class)
    ("CCCCN=C=O", "butyl isocyanate", "异氰酸丁酯"),
    ("CN=C=O", "methyl isocyanate", "异氰酸甲酯"),
    # positive: aryl isocyanate / isothiocyanate
    ("c1ccccc1N=C=O", "isocyanatobenzene", "异氰酸苯酯"),
    ("c1ccccc1N=C=S", "isothiocyanatobenzene", "异硫氰酸苯酯"),
    # positive: fused iso keeps relative locants (not mono-benzene omit)
    (
        "FC(C1=CC=C(C=C1)N=C=S)(F)F",
        "4-(trifluoromethyl)isothiocyanatobenzene",
        "4-三氟甲基异硫氰酸苯酯",
    ),
    ("Cc1ccc(N=C=S)cc1", "4-methylisothiocyanatobenzene", "4-甲基异硫氰酸苯酯"),
    # positive: multi-sub arene (prefix style; locants match gold EN)
    (
        "CC1=C(C=CC(=C1)C)N=C=S",
        "2,4-dimethylisothiocyanatobenzene",
        "2,4-二甲基异硫氰酸苯酯",
    ),
    (
        "ClC1=C(C=C(C=C1)C)N=C=O",
        "1-chloro-2-isocyanato-4-methylbenzene",
        "1-氯-2-异氰酸根合-4-甲基苯",
    ),
    # negative: true amide / nitrile / nitro must not regress
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CCCNC=O", "N-propylformamide", "N-丙基甲酰胺"),
    ("c1ccccc1C#N", "benzonitrile", "苯甲腈"),
    ("c1ccccc1[N+](=O)[O-]", "nitrobenzene", "硝基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_isocyanate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
