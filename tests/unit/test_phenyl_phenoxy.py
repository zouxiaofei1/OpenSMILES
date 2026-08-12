# IUPAC: P-29.3 / P-14.3.1
# Layer: L2,L3
"""Unsubstituted phenyl / phenoxy (depth-1 aryl) on benzene parents.

Parent = benzene; ring–O–Ph → phenoxy; ring–Ph → phenyl.
Ph may carry at most one ring halo (locant from attachment = 1).
Must not break anisole / ethoxybenzene / bare benzene.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: biphenyl as benzene + phenyl
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    # negative: retained anisole / ethoxy / bare benzene
    ("CCOc1ccccc1", "ethoxybenzene", None),
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_phenyl_phenoxy(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
