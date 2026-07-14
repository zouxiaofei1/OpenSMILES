# IUPAC: P-65.1.2 / P-31.1 / P-93.4
# Layer: L2,L3,L4,L5
"""Hydroxyalkenoic acids / hydroxyalkenoates (P-65.1.2, P-31.1, P-93.4).

Carboxyl is the principal characteristic group (suffix -enoic acid / -enoate);
aliphatic OH is a hydroxy prefix. Non-aromatic C=C enters the parent as enoic
(not saturated alkanoic). Parent chain covers COOH + all C=C carbons + OH carbon.
E/Z stereo prefixes precede substituent prefixes.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: mono-ene hydroxyalkenoic acids
    (
        r"OC/C=C/C(=O)O",
        "(E)-4-hydroxybut-2-enoic acid",
        "(E)-4-羟基丁-2-烯酸",
    ),
    (
        "OCCCCCCC=CC(=O)O",
        "9-hydroxynon-2-enoic acid",
        "9-羟基壬-2-烯酸",
    ),
    (
        "C=C(O)C(=O)O",
        "2-hydroxyprop-2-enoic acid",
        "2-羟基丙-2-烯酸",
    ),
    (
        r"CCCCCCCC[C@H](O)/C=C/CCCCCCC(=O)O",
        "(E,10S)-10-hydroxyoctadec-8-enoic acid",
        "(E,10S)-10-羟基十八-8-烯酸",
    ),
    (
        r"CCCCCCC(O)C/C=C\CCCCCCCC(=O)[O-]",
        "(Z)-12-hydroxyoctadec-9-enoate",
        "(Z)-12-羟基十八-9-烯酸根",
    ),
    (
        r"C[C@@H](O)CCCC/C=C/C(=O)O",
        "(E,8R)-8-hydroxynon-2-enoic acid",
        "(E,8R)-8-羟基壬-2-烯酸",
    ),
    # positive: multi-ene hydroxyalkenoic acid
    (
        r"CCCCC[C@H](O)/C=C/C=C\CCCCCCCC(=O)O",
        "(9Z,11E,13S)-13-hydroxyoctadeca-9,11-dienoic acid",
        "(9Z,11E,13S)-13-羟基十八-9,11-二烯酸",
    ),
    # negatives: plain alkenoic, saturated hydroxy acid, alkenol, carboxylate
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    (
        "CCCCCCCCCCC(O)C(=O)O",
        "2-hydroxydodecanoic acid",
        "2-羟基十二酸",
    ),
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
    ("CCCCCCCCCCCC(=O)[O-]", "dodecanoate", "十二酸根"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hydroxy_alkenoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
