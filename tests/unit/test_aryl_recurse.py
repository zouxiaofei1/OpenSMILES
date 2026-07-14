# IUPAC: P-29.3 / P-14.3.1
# Layer: L2,L3
"""True recursive aryl substituent namer (depth can exceed 2).

Nested Ph / OPh on a Ph arm may themselves carry simple leaves (halo/Me/
alkoxy/nitro/CF3/OH/NH2) or further nested unsub Ph within max depth.
Bare biphenyl and depth-1/2 arms must stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # depth-3: nested Ph with simple leaves
    (
        "Clc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorophenyl)phenyl]ethanol",
        "2-[4-(4-氯苯基)苯基]乙醇",
    ),
    (
        "Cc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-methylphenyl)phenyl]ethanol",
        "2-[4-(4-甲基苯基)苯基]乙醇",
    ),
    (
        "COc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-methoxyphenyl)phenyl]ethanol",
        "2-[4-(4-甲氧基苯基)苯基]乙醇",
    ),
    (
        "O=[N+]([O-])c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-nitrophenyl)phenyl]ethanol",
        "2-[4-(4-硝基苯基)苯基]乙醇",
    ),
    (
        "Clc1cc(Cl)ccc1-c2ccc(CCO)cc2",
        "2-[4-(2,4-dichlorophenyl)phenyl]ethanol",
        "2-[4-(2,4-二氯苯基)苯基]乙醇",
    ),
    # nested phenoxy (unsub / mono-halo)
    (
        "c1ccc(Oc2ccc(CCO)cc2)cc1",
        "2-(4-phenoxyphenyl)ethanol",
        "2-(4-苯氧基苯基)乙醇",
    ),
    (
        "Clc1ccc(Oc2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorophenoxy)phenyl]ethanol",
        "2-[4-(4-氯苯氧基)苯基]乙醇",
    ),
    # prior depth-2 unsub nested Ph
    (
        "c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-(4-phenylphenyl)ethanol",
        "2-(4-苯基苯基)乙醇",
    ),
    # negatives
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_aryl_recurse(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
