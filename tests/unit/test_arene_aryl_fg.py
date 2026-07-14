# IUPAC: P-29.3 / P-14.3.1 / P-63.1.4 / P-62.2.1 / P-66.6.1 / P-65.1.1.1
# Layer: L2,L3
"""Depth-1 phenyl / phenoxy on retained arene FG parents.

Parent = phenol / aniline / benzaldehyde / benzoic acid (FG ring = 1).
Ring–O–Ph → phenoxy; ring–Ph → phenyl. Ph: unsub or mono-halo only.
Must not regress bare FG parents, anisole, or phenoxybenzene.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: phenoxy on benzaldehyde / benzoic / phenol / aniline
    ("O=Cc1ccc(Oc2ccccc2)cc1", "4-phenoxybenzaldehyde", "4-苯氧基苯甲醛"),
    ("O=C(O)c1ccc(Oc2ccccc2)cc1", "4-phenoxybenzoic acid", "4-苯氧基苯甲酸"),
    ("Oc1ccc(Oc2ccccc2)cc1", "4-phenoxyphenol", "4-苯氧基苯酚"),
    ("Nc1ccc(Oc2ccccc2)cc1", "4-phenoxyaniline", "4-苯氧基苯胺"),
    # positive: phenyl on benzaldehyde
    ("O=Cc1ccc(-c2ccccc2)cc1", "4-phenylbenzaldehyde", "4-苯基苯甲醛"),
    # positive: multi-sub gold aniline (halo + nitro + phenoxy)
    (
        "Nc1c([N+](=O)[O-])ccc(Oc2ccccc2)c1Cl",
        "2-chloro-6-nitro-3-phenoxyaniline",
        None,
    ),
    # positive: Ph written first — FG ring still parent
    ("c1ccc(Oc2ccc(O)cc2)cc1", "4-phenoxyphenol", "4-苯氧基苯酚"),
    # negatives: bare FG parents and existing benzene aryl / alkoxy
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_aryl_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_phenoxy_benzaldehyde_not_chain() -> None:
    r = SMILESNNamer().name("O=Cc1ccc(Oc2ccccc2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "heptanal" not in en
    assert "phenoxybenzaldehyde" in en
