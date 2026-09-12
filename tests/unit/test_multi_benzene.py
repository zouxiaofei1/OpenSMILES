# IUPAC: P-14.3.4 / P-22.1.3
# Layer: L2,L4,L5
"""Multi-substituted benzene: ring halo and/or methyl (C1), 2–3 substituents.

Lowest set of locants (ring rotation). Exactly two methyls → retained xylene
(en); Chinese remains systematic dimethylbenzene form.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: di/tri halo and/or methyl on benzene
    ("Clc1ccc(Cl)cc1", "1,4-dichlorobenzene", "1,4-二氯苯"),
    ("Clc1cccc(Cl)c1", "1,3-dichlorobenzene", "1,3-二氯苯"),
    ("Fc1ccc(F)cc1", "1,4-difluorobenzene", "1,4-二氟苯"),
    ("Fc1ccc(Br)cc1", "1-bromo-4-fluorobenzene", "1-溴-4-氟苯"),
    ("Cc1cc(C)cc(C)c1", "1,3,5-trimethylbenzene", "1,3,5-三甲基苯"),
    ("Cc1ccc(Cl)cc1", "1-chloro-4-methylbenzene", "1-氯-4-甲基苯"),
    ("Clc1ccc(Cl)c(Cl)c1", "1,2,4-trichlorobenzene", "1,2,4-三氯苯"),
    # negative: mono / unsubstituted / non-aromatic keep existing names
    ("c1ccccc1", "benzene", "苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("C1CCCCC1", "cyclohexane", "环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_multi_benzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C=Cc1ccccc1", "C#Cc1ccccc1"])
def test_unsaturated_sidechain_not_ethylbenzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != "ethylbenzene"
