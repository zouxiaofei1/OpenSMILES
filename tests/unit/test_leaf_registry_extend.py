# IUPAC: P-29.3
# Layer: L2
"""Registry leaves: n-alkyl C2–C4, methylthio, cyano; nested benzyl skeleton.

Ph arm on chain alcohol/amine/acid; Ar-CN must not become formonitrile parent.
Prior halo/MeO/nested Ph and bare arene alkyl parents stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # n-alkyl leaves
    ("CCc1ccc(CCO)cc1", "2-(4-ethylphenyl)ethanol", "2-(4-乙基苯基)乙醇"),
    ("CCCc1ccc(CCO)cc1", "2-(4-propylphenyl)ethanol", "2-(4-丙基苯基)乙醇"),
    ("CCCCc1ccc(CCO)cc1", "2-(4-butylphenyl)ethanol", "2-(4-丁基苯基)乙醇"),
    # methylthio
    (
        "CSc1ccc(CCO)cc1",
        "2-(4-methylsulfanylphenyl)ethanol",
        "2-(4-甲硫基苯基)乙醇",
    ),
    # cyano (alcohol parent, not formonitrile)
    ("N#Cc1ccc(CCO)cc1", "2-(4-cyanophenyl)ethanol", "2-(4-氰基苯基)乙醇"),
    ("N#Cc1ccc(CO)cc1", "(4-cyanophenyl)methanol", "(4-氰基苯基)甲醇"),
    # nested benzyl
    (
        "c1ccc(Cc2ccc(CCO)cc2)cc1",
        "2-(4-benzylphenyl)ethanol",
        "2-(4-苄基苯基)乙醇",
    ),
    (
        "Clc1ccc(Cc2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorobenzyl)phenyl]ethanol",
        "2-[4-(4-氯苄基)苯基]乙醇",
    ),
    # negatives
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("c1ccc(-c2ccc(CCO)cc2)cc1", "2-(4-phenylphenyl)ethanol", "2-(4-苯基苯基)乙醇"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_leaf_registry_extend(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
