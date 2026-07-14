# IUPAC: P-29.3
# Layer: L2,L3
"""Depth-2 hydroxy / amino leaves on Ph arms; aromatic OH/NH2 not chain poly-FG.

Chain alcohol/amine/acid stays parent; phenolic OH and aniline-NH2 are arm leaves.
True alkanediol / diamine / phenol / aniline must remain correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: hydroxyphenyl on chain parents
    ("Oc1ccc(CCO)cc1", "2-(4-hydroxyphenyl)ethanol", "2-(4-羟基苯基)乙醇"),
    ("Oc1cccc(CCO)c1", "2-(3-hydroxyphenyl)ethanol", "2-(3-羟基苯基)乙醇"),
    ("Oc1ccc(CO)cc1", "(4-hydroxyphenyl)methanol", "(4-羟基苯基)甲醇"),
    (
        "O=C(O)Cc1ccc(O)cc1",
        "2-(4-hydroxyphenyl)acetic acid",
        "2-(4-羟基苯基)乙酸",
    ),
    # positive: aminophenyl on chain parents
    ("Nc1ccc(CCO)cc1", "2-(4-aminophenyl)ethanol", "2-(4-氨基苯基)乙醇"),
    ("Nc1ccc(CO)cc1", "(4-aminophenyl)methanol", "(4-氨基苯基)甲醇"),
    (
        "O=C(O)Cc1ccc(N)cc1",
        "2-(4-aminophenyl)acetic acid",
        "2-(4-氨基苯基)乙酸",
    ),
    (
        "Nc1ccccc1CC(=O)O",
        "2-(2-aminophenyl)acetic acid",
        "2-(2-氨基苯基)乙酸",
    ),
    # positive: hydroxy + halo mix on arm
    (
        "Oc1ccc(CCO)c(Cl)c1",
        "2-(2-chloro-4-hydroxyphenyl)ethanol",
        "2-(2-氯-4-羟基苯基)乙醇",
    ),
    # negative: true polyols / phenol / aniline / prior depth-2
    ("OCCO", "ethane-1,2-diol", "乙烷-1,2-二醇"),
    ("NCCN", "ethane-1,2-diamine", "乙烷-1,2-二胺"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("Oc1ccc(CC)cc1", "4-ethylphenol", "4-乙基苯酚"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_aryl_depth2_oh_nh2(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
