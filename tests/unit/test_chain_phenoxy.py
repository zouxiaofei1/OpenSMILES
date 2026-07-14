# IUPAC: P-29.3 / P-62.2 / P-63.2.2 / P-65.1
# Layer: L2,L3,L5
"""Phenoxy (and multi-sub Ph) on chain amine / alcohol / acid parents.

Plan 3.1: recursive aryl side on non-arene parents (sec-amine + aryloxy-acid).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # primary already works — regression
    ("c1ccc(OCCN)cc1", "phenoxyethanamine", None),
    ("Clc1ccc(OCCN)cc1", "2-(4-chlorophenoxy)ethanamine", None),
    ("OCCOc1ccccc1", "phenoxyethanol", None),
    # sec amine with ring (was 1-chlorohexane)
    (
        "ClC=1C=C(OCCNC)C=CC1",
        "2-(3-chlorophenoxy)-N-methylethanamine",
        "2-(3-氯苯氧基)-N-甲基乙胺",
    ),
    # acid with chloro+methyl phenoxy (was propanoic acid)
    (
        "Cc1cc(Cl)ccc1OC(C)C(=O)O",
        "2-(4-chloro-2-methylphenoxy)propanoic acid",
        "2-(4-氯-2-甲基苯氧基)丙酸",
    ),
    # already worked
    (
        "CC(C)(Oc1ccccc1)C(=O)O",
        "2-methyl-2-phenoxypropanoic acid",
        None,
    ),
    # negatives
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", None),
    ("CCNC", "N-methylethanamine", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_chain_phenoxy(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
