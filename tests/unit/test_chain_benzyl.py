# IUPAC: P-29.3 / P-63.1 / P-65.1 / P-62.2
# Layer: L2,L3,L5
"""Benzyl / phenyl on chain acid / amine parents (plan 3.4)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # diacid + C-benzyl / phenyl (open chain no longer walks into arene)
    (
        "O=C(O)CC(Cc1ccccc1)C(=O)O",
        "2-benzylbutanedioic acid",
        "2-苄基丁二酸",
    ),
    # phenylmethanamine system name
    ("NCc1ccccc1", "phenylmethanamine", "苯基甲胺"),
    ("c1ccc(CCN)cc1", "phenylethanamine", "苯基乙胺"),
    # aminoalcohol N-benzyl
    (
        "ClC=1C=C(CNCCO)C=CC1Cl",
        "2-[(3,4-dichlorobenzyl)amino]ethanol",
        "2-[(3,4-二氯苄基)氨基]乙醇",
    ),
    # regressions
    ("c1ccc(OCCN)cc1", "phenoxyethanamine", None),
    ("CC(=O)O", "acetic acid", None),
    ("NCCO", "aminoethanol", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_chain_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
