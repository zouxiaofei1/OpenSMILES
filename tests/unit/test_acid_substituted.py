# 合并自 6 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_aminoalkanoic_acid.py: Open-chain saturated monoaminoalkanoic acids (carboxylic acid parent + amino).
test_hydroxyalkanoic_acid.py: Open-chain saturated monohydroxyalkanoic acids (carboxylic acid parent + hydroxy).
test_hydroxy_alkenoic.py: Hydroxyalkenoic acids / hydroxyalkenoates (P-65.1.2, P-31.1, P-93.4).
test_oxoalkanoic_acid.py: Open-chain saturated mono-oxoalkanoic acids (carboxylic acid parent + oxo).
test_prefix_alkanedioic_acid.py: Open-chain saturated substituted alkanedioic acids / dioates.
test_prefix_alkenoic.py: Amino / oxo / hydroxy prefixes on alkenoic acids (P-65.1.2, P-31.1, P-93.4).
"""
from __future__ import annotations

import pytest

from opensmiles.layer1.analyzer import analyze
from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# 合并自 test_aminoalkanoic_acid.py
# IUPAC: P-65.1.2
# Layer: L3,L5
#
# Open-chain saturated monoaminoalkanoic acids (carboxylic acid parent + amino).
#
# Carboxyl is principal characteristic group; primary amino is amino prefix
# with locant from carboxyl numbering (COOH carbon = 1).
# ==========================================================================
aminoalkanoic_acid__CASES = [
    # positive: α/β/ω-amino monoacids
    ("NCC(=O)O", "2-aminoacetic acid", "2-氨基乙酸"),
    ("CCCCC(N)C(=O)O", "2-aminohexanoic acid", None),
    ("NCCCC(=O)O", "4-aminobutanoic acid", "4-氨基丁酸"),
    ("OC(=O)C(O)C", "2-hydroxypropanoic acid", "2-羟基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", aminoalkanoic_acid__CASES)
def test_aminoalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_hydroxyalkanoic_acid.py
# IUPAC: P-65.1.2
# Layer: L2,L3,L5
#
# Open-chain saturated monohydroxyalkanoic acids (carboxylic acid parent + hydroxy).
#
# Carboxyl is principal characteristic group; non-carboxyl OH is hydroxy prefix
# with locant from carboxyl numbering (COOH carbon = 1).
# ==========================================================================
hydroxyalkanoic_acid__CASES = [
    ("OC(=O)CCO", "3-hydroxypropanoic acid", "3-羟基丙酸"),
    ("OC(=O)CC(O)C", "3-hydroxybutanoic acid", "3-羟基丁酸"),
    ("O=C(O)CCCCCO", "6-hydroxyhexanoic acid", "6-羟基己酸"),
    ("OC(=O)C(O)CC", "2-hydroxybutanoic acid", "2-羟基丁酸"),
]


@pytest.mark.parametrize("smiles,en,zh", hydroxyalkanoic_acid__CASES)
def test_hydroxyalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_hydroxy_alkenoic.py
# IUPAC: P-65.1.2 / P-31.1 / P-93.4
# Layer: L2,L3,L4,L5
#
# Hydroxyalkenoic acids / hydroxyalkenoates (P-65.1.2, P-31.1, P-93.4).
#
# Carboxyl is the principal characteristic group (suffix -enoic acid / -enoate);
# aliphatic OH is a hydroxy prefix. Non-aromatic C=C enters the parent as enoic
# (not saturated alkanoic). Parent chain covers COOH + all C=C carbons + OH carbon.
# E/Z stereo prefixes precede substituent prefixes.
# ==========================================================================
hydroxy_alkenoic__CASES = [
    # positive: mono-ene hydroxyalkenoic acids
    (
        r"OC/C=C/C(=O)O",
        "(2E)-4-hydroxybut-2-enoic acid",
        "(2E)-4-羟基丁-2-烯酸",
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
        "(8E,10S)-10-hydroxyoctadec-8-enoic acid",
        "(8E,10S)-10-羟基十八-8-烯酸",
    ),
    (
        r"CCCCCCC(O)C/C=C\CCCCCCCC(=O)[O-]",
        "(9Z)-12-hydroxyoctadec-9-enoate",
        "(9Z)-12-羟基十八-9-烯酸根",
    ),

    # positive: multi-ene hydroxyalkenoic acid
    (
        r"CCCCC[C@H](O)/C=C/C=C\CCCCCCCC(=O)O",
        "(9Z,11E,13S)-13-hydroxyoctadeca-9,11-dienoic acid",
        "(9Z,11E,13S)-13-羟基十八-9,11-二烯酸",
    ),

    (
        "CCCCCCCCCCC(O)C(=O)O",
        "2-hydroxydodecanoic acid",
        "2-羟基十二酸",
    ),
]


@pytest.mark.parametrize("smiles,en,zh", hydroxy_alkenoic__CASES)
def test_hydroxy_alkenoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_oxoalkanoic_acid.py
# IUPAC: P-65.1.2
# Layer: L3
#
# Open-chain saturated mono-oxoalkanoic acids (carboxylic acid parent + oxo).
#
# Carboxyl is principal characteristic group; ketone carbonyl is oxo/氧代 prefix
# with locant from carboxyl numbering (COOH carbon = 1).
# ==========================================================================
oxoalkanoic_acid__CASES = [
    ("OC(=O)CCC(=O)C", "4-oxopentanoic acid", "4-氧代戊酸"),
    ("CC(=O)CC(=O)O", "3-oxobutanoic acid", "3-氧代丁酸"),
    ("CCCC(=O)C(=O)O", "2-oxopentanoic acid", None),
]


@pytest.mark.parametrize("smiles,en,zh", oxoalkanoic_acid__CASES)
def test_oxoalkanoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_prefix_alkanedioic_acid.py
# IUPAC: P-65.1.2
# Layer: L2,L3,L4,L5
#
# Open-chain saturated substituted alkanedioic acids / dioates.
#
# Parent remains alkanedioic acid / …dioate (or oxalic); chain hydroxy /
# amino / oxo are prefixes (incl. multi-prefix); numbering from either
# carboxyl carbon by lowest set. Anion → …dioate / …二酸根.
# ==========================================================================
prefix_alkanedioic_acid__CASES = [
    # positive: hydroxy / amino / oxo prefixes on saturated diacids
    ("O=C(O)CC(O)CC(=O)O", "3-hydroxypentanedioic acid", "3-羟基戊二酸"),
    ("NC(CC(=O)O)CC(=O)O", "3-aminopentanedioic acid", "3-氨基戊二酸"),
    ("O=C(O)C(N)CC(=O)O", "2-aminobutanedioic acid", "2-氨基丁二酸"),
    ("O=C(O)C(O)C(=O)O", "2-hydroxypropanedioic acid", "2-羟基丙二酸"),
    ("O=C([O-])CCCC(=O)C(=O)[O-]", "2-oxohexanedioate", "2-氧代己二酸根"),
    ("O=C(O)C(=O)C(O)C(=O)O", "2-hydroxy-3-oxobutanedioic acid", "2-羟基-3-氧代丁二酸"),
    ("O=C(O)C(O)C(O)C(O)C(=O)O", "2,3,4-trihydroxypentanedioic acid", "2,3,4-三羟基戊二酸"),
    ("O=C(O)C(C)(C)C(O)C(=O)O", "3-hydroxy-2,2-dimethylbutanedioic acid", "3-羟基-2,2-二甲基丁二酸"),
    # negative: unsubstituted diacid, mono hydroxyacid, alkyl-diacid, monoacid
    ("O=C(O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("O=C(O)C(O)C", "2-hydroxypropanoic acid", "2-羟基丙酸"),
    ("O=C(O)C(C)C(=O)O", "2-methylpropanedioic acid", "2-甲基丙二酸"),
]


@pytest.mark.parametrize("smiles,en,zh", prefix_alkanedioic_acid__CASES)
def test_prefix_alkanedioic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_prefix_alkenoic.py
# IUPAC: P-65.1.2 / P-31.1 / P-93.4
# Layer: L2,L3,L5
#
# Amino / oxo / hydroxy prefixes on alkenoic acids (P-65.1.2, P-31.1, P-93.4).
#
# Carboxyl is the principal characteristic group (-enoic acid / -enoate).
# Aliphatic amino, oxo (ketone), and hydroxy are prefixes — they must not force
# a saturated alkanoic parent. Parent chain covers COOH + all C=C + prefix carbons.
# E/Z multi-prefixes precede substituent prefixes. Anion → -oate / 酸根.
# ==========================================================================
prefix_alkenoic__CASES = [
    # positive: aminoalkenoic acids
    (
        "NCC=CC(=O)O",
        "4-aminobut-2-enoic acid",
        "4-氨基丁-2-烯酸",
    ),
    (
        "N[C@H](C(=O)O)CC=C",
        "(2S)-2-aminopent-4-enoic acid",
        "(2S)-2-氨基戊-4-烯酸",
    ),
    (
        "C=C(Cl)C[C@H](N)C(=O)O",
        "(2S)-2-amino-4-chloropent-4-enoic acid",
        "(2S)-2-氨基-4-氯戊-4-烯酸",
    ),
    (
        r"NC/C=C/C(=O)O",
        "(2E)-4-aminobut-2-enoic acid",
        "(2E)-4-氨基丁-2-烯酸",
    ),
    # positive: oxoalkenoic (mono + multi-ene anion)
    (
        r"CC(=O)/C=C/C(=O)O",
        "(2E)-4-oxopent-2-enoic acid",
        "(2E)-4-氧代戊-2-烯酸",
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
]


@pytest.mark.parametrize("smiles,en,zh", prefix_alkenoic__CASES)
def test_prefix_alkenoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
