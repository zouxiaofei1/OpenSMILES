# IUPAC: P-65.6
# Layer: L2,L5
"""Ester alkoxy sides: benzyl / tert-butyl / isopropyl / phenyl (P-65.6).

O-CH2Ph is benzyl not methyl; O-C(CH3)3 is tert-butyl not ethyl.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: special alkoxy
    ("CC(=O)OCC1=CC=CC=C1", "benzyl acetate", "乙酸苄酯"),
    ("C(CCCCCCCCCCCCCCCCC)(=O)OCC1=CC=CC=C1", "benzyl octadecanoate", None),
    ("CC(=O)OC(C)(C)C", "tert-butyl acetate", "乙酸叔丁酯"),
    ("CC(=O)OC(C)C", "propan-2-yl acetate", "乙酸丙-2-基酯"),
    ("CC(=O)Oc1ccccc1", "phenyl acetate", "乙酸苯酯"),
    ("O=C(OCc1ccccc1)c1ccccc1", "benzyl benzoate", "苯甲酸苄酯"),
    # negative: linear alkyl esters must hold
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("CC(=O)OCCC", "propyl acetate", "乙酸丙酯"),
]

# Must NOT collapse to bare phenyl/benzyl acetate (substituted/heteroaryl)
NOT_BARE = [
    "CC(=O)Oc1ccc(Cl)cc1",
    "CC(=O)Oc1ncccc1",
    "CC(=O)Oc1ccccc1C",
    "CC(=O)OCC1=CC=C(Cl)C=C1",
    "CC(=O)OCc1ccccn1",
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ester_alkoxy_sides(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", NOT_BARE)
def test_not_bare_phenyl_or_benzyl(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en not in {
        normalize_en("phenyl acetate"),
        normalize_en("benzyl acetate"),
    }
