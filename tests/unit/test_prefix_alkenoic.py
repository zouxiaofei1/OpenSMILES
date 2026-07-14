# IUPAC: P-65.1.2 / P-31.1 / P-93.4
# Layer: L2,L3,L5
"""Amino / oxo / hydroxy prefixes on alkenoic acids (P-65.1.2, P-31.1, P-93.4).

Carboxyl is the principal characteristic group (-enoic acid / -enoate).
Aliphatic amino, oxo (ketone), and hydroxy are prefixes — they must not force
a saturated alkanoic parent. Parent chain covers COOH + all C=C + prefix carbons.
E/Z multi-prefixes precede substituent prefixes. Anion → -oate / 酸根.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: aminoalkenoic acids
    (
        "NCC=CC(=O)O",
        "4-aminobut-2-enoic acid",
        "4-氨基丁-2-烯酸",
    ),
    (
        "N[C@H](C(=O)O)CC=C",
        "2-aminopent-4-enoic acid",
        "2-氨基戊-4-烯酸",
    ),
    (
        "C=C(Cl)C[C@H](N)C(=O)O",
        "2-amino-4-chloropent-4-enoic acid",
        "2-氨基-4-氯戊-4-烯酸",
    ),
    (
        r"NC/C=C/C(=O)O",
        "(E)-4-aminobut-2-enoic acid",
        "(E)-4-氨基丁-2-烯酸",
    ),
    # positive: oxoalkenoic (mono + multi-ene anion)
    (
        r"CC(=O)/C=C/C(=O)O",
        "(E)-4-oxopent-2-enoic acid",
        "(E)-4-氧代戊-2-烯酸",
    ),
    (
        r"CCCCC/C=C\C/C=C\C=C\C(=O)C/C=C\CCCC(=O)[O-]",
        "(5Z,9E,11Z,14Z)-8-oxoicosa-5,9,11,14-tetraenoate",
        "(5Z,9E,11Z,14Z)-8-氧代二十-5,9,11,14-四烯酸根",
    ),
    # positive: hydroxy + oxo multi-ene
    (
        r"O=C(O)C(=O)/C=C/C=C\O",
        "(3E,5Z)-6-hydroxy-2-oxohexa-3,5-dienoic acid",
        "(3E,5Z)-6-羟基-2-氧代己-3,5-二烯酸",
    ),
    # negatives: saturated amino/oxo acids, plain alkenoic, hydroxy alkenoic, alkenol
    ("NCCCC(=O)O", "4-aminobutanoic acid", "4-氨基丁酸"),
    ("CC(=O)CC(=O)O", "3-oxobutanoic acid", "3-氧代丁酸"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    (
        r"OC/C=C/C(=O)O",
        "(E)-4-hydroxybut-2-enoic acid",
        "(E)-4-羟基丁-2-烯酸",
    ),
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_prefix_alkenoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
