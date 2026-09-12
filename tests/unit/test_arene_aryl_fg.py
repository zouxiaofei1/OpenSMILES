# IUPAC: P-29.3 / P-14.3.1 / P-63.1.4 / P-62.2.1 / P-66.6.1 / P-65.1.1.1
# Layer: L2,L3
"""Depth-1 phenyl / phenoxy on retained arene FG parents.

Parent = phenol / aniline / benzaldehyde / benzoic acid (FG ring = 1).
Ring–O–Ph → phenoxy; ring–Ph → phenyl. Ph: unsub or mono-halo only.
Must not regress bare FG parents, anisole, or phenoxybenzene.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: phenyl on benzaldehyde
    ("O=Cc1ccc(-c2ccccc2)cc1", "4-phenylbenzaldehyde", "4-苯基苯甲醛"),
    # negatives: bare FG parents and existing benzene aryl / alkoxy
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_aryl_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
