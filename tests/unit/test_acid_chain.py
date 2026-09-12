# 合并自 10 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_carboxylic_acid.py: Monocarboxylic acids: systematic oic acid + retained formic/acetic.
test_alkanedioic_acid.py: Unsubstituted open-chain saturated alkanedioic acids (exactly two COOH).
test_alkadienedioic.py: Open-chain polyunsaturated dicarboxylic acids (alkadienedioic).
test_alkenedioic.py: Open-chain monounsaturated dicarboxylic acids (alkenedioic acids).
test_alkenoic_acid.py: Open-chain monounsaturated monocarboxylic acids (alkenoic acids).
test_alkynoic.py: Open-chain monounsaturated monocarboxylic alkynoic acids / alkynoates.
test_polyenoic_acid.py: Polyunsaturated acyclic mono-FG parents (dienoic / trienoic / …).
test_tricarboxylic.py: Open-chain tricarboxylic acids and their fully deprotonated anions.
test_carboxylate_anion.py: Carboxylate anions: C(=O)[O-] → alkanoate / …酸根 (not aldehyde).
test_metal_carboxylate.py: Alkali metal carboxylates: salt dissociation + functional class salt names.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_mono_carboxylic_acid.py
# IUPAC: P-65.1.1
# Layer: L1,L2,L4,L5
#
# Monocarboxylic acids: systematic oic acid + retained formic/acetic.
#
# Carboxyl carbon is part of the parent chain (locant 1); suffix has no locant.
# ==========================================================================
mono_carboxylic_acid__CASES = [
    ("CCCCC(=O)O", "pentanoic acid", "戊酸"),
    ("CC(C)C(=O)O", "2-methylpropanoic acid", "2-甲基丙酸"),
    ("CC(C)CCC(=O)O", "4-methylpentanoic acid", "4-甲基戊酸"),
    ("CCCl", "chloroethane", "氯乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_carboxylic_acid__CASES)
def test_mono_carboxylic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanedioic_acid.py
# IUPAC: P-65.1.1
# Layer: L2,L5
#
# Unsubstituted open-chain saturated alkanedioic acids (exactly two COOH).
#
# Parent chain through both carboxyl carbons; name alkanedioic acid /
# {烷首}二酸. C2 retained: oxalic acid / 草酸. No locants for terminal diacids.
# ==========================================================================
alkanedioic_acid__CASES = [
    # positive: unsubstituted open-chain diacids
    ("OC(=O)C(=O)O", "oxalic acid", "草酸"),
    ("O=C(O)CCCCCCCC(=O)O", "nonanedioic acid", "壬二酸"),
    ("CCC(=O)O", "propanoic acid", "丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanedioic_acid__CASES)
def test_alkanedioic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkadienedioic.py
# IUPAC: P-65.1.1 / P-31.1 / P-72.2.2.1
# Layer: L2,L4,L5
#
# Open-chain polyunsaturated dicarboxylic acids (alkadienedioic).
#
# Exactly two carboxyls (acid or anion) + ≥2 open-chain C=C; parent covers both
# carboxyl carbons and all double-bond carbons. English:
# (2E,4E)-hexa-2,4-dienedioic acid / …dienedioate; Chinese: …-2,4-二烯二酸 / 根.
# Reuses multi E/Z prefix and diacid anion conversion.
# ==========================================================================
alkadienedioic__CASES = [
    # positive: open-chain dienedioic acids / anions with multi E/Z
    (
        r"O=C([O-])/C=C/C=C/C(=O)[O-]",
        "(2E,4E)-hexa-2,4-dienedioate",
        "(2E,4E)-己-2,4-二烯二酸根",
    ),
    (
        r"O=C([O-])/C=C\C=C/C(=O)[O-]",
        "(2Z,4Z)-hexa-2,4-dienedioate",
        "(2Z,4Z)-己-2,4-二烯二酸根",
    ),
    (
        r"O=C(O)/C=C/C=C/C(=O)O",
        "(2E,4E)-hexa-2,4-dienedioic acid",
        "(2E,4E)-己-2,4-二烯二酸",
    ),
    # negative: saturated diacid anion and mono-ene diacid must not regress
    (r"O=C([O-])CCC(=O)[O-]", "butanedioate", "丁二酸根"),
]


@pytest.mark.parametrize("smiles,en,zh", alkadienedioic__CASES)
def test_alkadienedioic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenedioic.py
# IUPAC: P-44 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated dicarboxylic acids (alkenedioic acids).
#
# Exactly two COOH + one C=C; parent chain through both carboxyl carbons and the
# double bond. English: (E/Z)-alk-n-enedioic acid; Chinese: (E/Z)-{烷首}-n-烯二酸.
# Stereo from RDKit BondStereo; no maleic/fumaric retained names.
# ==========================================================================
alkenedioic__CASES = [
    (r"OC(=O)/C=C\C(=O)O", "(2Z)-but-2-enedioic acid", "(2Z)-丁-2-烯二酸"),
    (r"O=C(O)/C=C/CCC(=O)O", "(2E)-hex-2-enedioic acid", None),
]


@pytest.mark.parametrize("smiles,en,zh", alkenedioic__CASES)
def test_alkenedioic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenoic_acid.py
# IUPAC: P-65.1.1 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated monocarboxylic acids (alkenoic acids).
#
# Carboxyl is principal characteristic group (COOH = locant 1); one non-aromatic
# C=C is expressed as -n-enoic / -n-烯酸 with the lower double-bond carbon locant.
# No (E)/(Z) stereodescriptors this cycle.
# ==========================================================================
alkenoic_acid__CASES = [
    ("C=CCC(=O)O", "but-3-enoic acid", "丁-3-烯酸"),
    ("CC=CC(=O)O", "but-2-enoic acid", "丁-2-烯酸"),
    ("CC(C)=CC(=O)O", "3-methylbut-2-enoic acid", "3-甲基丁-2-烯酸"),
    ("CCC=CCC(=O)O", "hex-3-enoic acid", "己-3-烯酸"),
    ("CCC=CC(=O)O", "pent-2-enoic acid", "戊-2-烯酸"),
    ("OC(=O)C(=O)C", "2-oxopropanoic acid", "2-氧代丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenoic_acid__CASES)
def test_alkenoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkynoic.py
# IUPAC: P-65.1.1 / P-31.1 / P-14.3.4
# Layer: L4,L5
#
# Open-chain monounsaturated monocarboxylic alkynoic acids / alkynoates.
#
# Carboxyl (or ester carbonyl) is principal FG (locant 1); one non-aromatic C≡C
# is expressed as -n-ynoic / -n-炔酸 (or -ynoate / 炔酸…酯). First scope: mono
# C≡C, no C=C, open-chain mono acid/ester.
#
# P-14.3.4 例外：prop-2-ynoic acid（丙-2-炔酸）虽无歧义也不省略位次，与
# prop-2-enoic acid 同列；propiolic acid（丙炔酸）已非保留名（P-65.1.1.2.4）。
# ==========================================================================
alkynoic__CASES = [
    # positive: acyclic mono-alkynoic acids
    ("C#CCC(=O)O", "but-3-ynoic acid", "丁-3-炔酸"),
    ("C#CC(=O)O", "prop-2-ynoic acid", "丙-2-炔酸"),
    ("C(CCCCCCC#C)(=O)O", "non-8-ynoic acid", "壬-8-炔酸"),
    ("C#CCC(=O)[O-]", "but-3-ynoate", "丁-3-炔酸根"),
    ("C#CC(=O)[O-]", "prop-2-ynoate", "丙-2-炔酸根"),
    # positive: alkynoates (esters)
    ("C(C#C)(=O)OCC", "ethyl prop-2-ynoate", "丙-2-炔酸乙酯"),
    ("C#CC(=O)OC", "methyl prop-2-ynoate", "丙-2-炔酸甲酯"),
    ("C/C=C/C(=O)O", "(2E)-but-2-enoic acid", "(2E)-丁-2-烯酸"),
    ("CC#C", "propyne", "丙炔"),
]


@pytest.mark.parametrize("smiles,en,zh", alkynoic__CASES)
def test_alkynoic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_polyenoic_acid.py
# IUPAC: P-31.1 / P-65.1.1
# Layer: L2, L4, L5
#
# Polyunsaturated acyclic mono-FG parents (dienoic / trienoic / …).
#
# Before: only mono-ene+FG was detected; 2+ C=C + mono FG fell back to saturated
# naming. L4/L5 already had ene_locants + polyalkenoic + E/Z; gap was L2.
# ==========================================================================
polyenoic_acid__CASES = [
    # dienoic acids
    ("C=CC=CC(=O)O", "penta-2,4-dienoic acid", "戊-2,4-二烯酸"),
    ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid", "己-2,4-二烯酸"),
    ("C=CC=CCC(=O)O", "hexa-3,5-dienoic acid", "己-3,5-二烯酸"),
    # mono-ene acid still works
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    # saturated acid not stolen
    ("CCCC(=O)O", "butanoic acid", "丁酸"),
    # alcohol polyene already had path; keep regression
    ("C=CC=CCO", "penta-2,4-dien-1-ol", "戊-2,4-二烯-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", polyenoic_acid__CASES)
def test_polyenoic_acid(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_polyenoic_not_saturated_collapse() -> None:
    r = SMILESNNamer().name("C=CC=CC(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert "dienoic" in en
    assert en != "pentanoic acid"


# ==========================================================================
# 合并自 test_tricarboxylic.py
# IUPAC: P-65.1.1 / P-72.2.2.1
# Layer: L2,L4,L5
#
# Open-chain tricarboxylic acids and their fully deprotonated anions.
# ==========================================================================
def test_tricarboxylate_does_not_fall_back_to_monoacid() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)[O-])CC(=O)[O-]")
    assert result.success
    assert "acetate" not in normalize_en(result.en)
    assert "hexanoate" not in normalize_en(result.en)


def test_partial_deprotonation_does_not_claim_polycarboxylate() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)O)CC(=O)O")
    assert "tricarboxylate" not in normalize_en(result.en)


# ==========================================================================
# 合并自 test_carboxylate_anion.py
# IUPAC: P-65.1.1 / P-72.2.2.1
# Layer: L1,L2,L5
#
# Carboxylate anions: C(=O)[O-] → alkanoate / …酸根 (not aldehyde).
#
# IUPAC P-65.1.1 parent acid chain numbering; P-72.2.2.1 anion suffix -ate / 酸根.
# ==========================================================================
carboxylate_anion__CASES = [
    # positive: bare carboxylate anions
    ("CCCCCCCCCCCCCC(=O)[O-]", "tetradecanoate", "十四酸根"),
    ("CC(C)CC(=O)[O-]", "3-methylbutanoate", "3-甲基丁酸根"),
    ("CC(=O)[O-]", "acetate", "乙酸根"),
    ("CCCCCCCCCCC(O)C(=O)[O-]", "2-hydroxydodecanoate", "2-羟基十二酸根"),
    ("NCCCC(=O)[O-]", "4-aminobutanoate", "4-氨基丁酸根"),
    ("C(=O)[O-]", "formate", "甲酸根"),
    ("CCCCCCCCCCCCCCCCCCCCCCCCCC", "hexacosane", "二十六烷"),
]


@pytest.mark.parametrize("smiles,en,zh", carboxylate_anion__CASES)
def test_carboxylate_anion(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_metal_carboxylate.py
# IUPAC: P-72.2.2.1 / P-65.1.1
# Layer: L0,L5
#
# Alkali metal carboxylates: salt dissociation + functional class salt names.
#
# L0 dissociates single alkali metal cation + one organic carboxylate fragment.
# L5 prefixes metal (EN) or replaces 酸根 → 酸钠/钾/锂 (ZH).
# Bare anions and neutral acids must not regress.
# ==========================================================================
metal_carboxylate__CASES = [
    ("[Na+].[O-]C(=O)c1ccccc1", "sodium benzoate", "苯甲酸钠"),
    ("FC1=C(C(=O)[O-])C=CC(=C1)F.[Na+]", "sodium 2,4-difluorobenzoate", "2,4-二氟苯甲酸钠"),
    ("O=C([O-])C.[K+]", "potassium acetate", "乙酸钾"),
    ("O=C[O-].[Na+]", "sodium formate", "甲酸钠"),
    # positive: metal after R/S (and optional E/Z) on anion stem
    ("C[C@H](O)C(=O)[O-].[Na+]", "sodium (2S)-2-hydroxypropanoate", "(2S)-2-羟基丙酸钠"),
    ("C/C=C/[C@H](O)C(=O)[O-].[Na+]", "sodium (2S,3E)-2-hydroxypent-3-enoate", "(2S,3E)-2-羟基戊-3-烯酸钠"),
    # negative: bare anion / neutral acid / ester must not regress
    ("O=C([O-])C", "acetate", "乙酸根"),
    ("C[C@H](O)C(=O)[O-]", "(2S)-2-hydroxypropanoate", "(2S)-2-羟基丙酸根"),
]


@pytest.mark.parametrize("smiles,en,zh", metal_carboxylate__CASES)
def test_metal_carboxylate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
