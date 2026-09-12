# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_monoalkyl_cycloalkane.py: Monosubstituted monocycloalkanes: methylcyclohexane etc.
test_polyalkyl_cycloalkane.py: Polyalkyl monocyclic cycloalkanes: ≥2 linear C1–C4 alkyls on one ring.
test_monohalo_cycloalkane.py: Monosubstituted monohalo monocycloalkanes: chlorocyclohexane etc.
test_sub_cyclo_fg.py: Saturated monocyclic mono-alcohol/ketone/amine with simple ring subs.
test_sub_cycloalkene_fg.py: Alkyl/halo monocyclic monoalkenes + cycloalkenol / cycloalkenone.
test_cycloalkane_exocyclic_fg.py: Monocyclic cycloalkane with one exocyclic carbonyl-class FG (P-65/P-66).
test_cycloalkene_exocyclic_fg.py: Monocyclic cycloalkene with one exocyclic carbonyl-class FG (acid/ester/aldehyde).
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_monoalkyl_cycloalkane.py
# IUPAC: P-14.3.4 / P-22.1.1
# Layer: L2,L3,L4,L5
#
# Monosubstituted monocycloalkanes: methylcyclohexane etc.
#
# Parent = saturated monocarbocycle; single linear alkyl side chain.
# Monosubstituted: omit locant (methylcyclohexane not 1-methyl…).
# ==========================================================================
monoalkyl_cycloalkane__CASES = [
    ("CC1CCCC1", "methylcyclopentane", "甲基环戊烷"),
    ("CC1CCC1", "methylcyclobutane", "甲基环丁烷"),
    ("CC1CC1", "methylcyclopropane", "甲基环丙烷"),
    ("CCC1CCCCC1", "ethylcyclohexane", "乙基环己烷"),
    ("C1CC1", "cyclopropane", "环丙烷"),
    ("CCCCCC", "hexane", "己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", monoalkyl_cycloalkane__CASES)
def test_monoalkyl_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_dimethyl_not_methylcyclohexane() -> None:
    """Disubstituted ring out of simple monoalkyl scope naming as methylcyclohexane."""
    r = SMILESNNamer().name("CC1CCCCC1C")
    assert normalize_en(r.en) != "methylcyclohexane"


# ==========================================================================
# 合并自 test_polyalkyl_cycloalkane.py
# IUPAC: P-14.3.4 / P-22.1.1
# Layer: L2,L4,L5
#
# Polyalkyl monocyclic cycloalkanes: ≥2 linear C1–C4 alkyls on one ring.
#
# Parent = saturated monocarbocycle; multiple linear alkyl side chains.
# Lowest set of locants (P-14.3.4); multiplicative prefixes (P-29.2).
# Monosubstituted still omits locant; unsubstituted and open-chain intact.
# ==========================================================================
polyalkyl_cycloalkane__CASES = [
    # positive: polyalkyl monocycloalkanes (keep locants + mult prefixes)
    ("CC1CCCCC1C", "1,2-dimethylcyclohexane", "1,2-二甲基环己烷"),
    ("CC1CCC(C)CC1", "1,4-dimethylcyclohexane", "1,4-二甲基环己烷"),
    ("CC1CCCC(C)C1", "1,3-dimethylcyclohexane", "1,3-二甲基环己烷"),
    ("CC1(C)CCCCC1", "1,1-dimethylcyclohexane", "1,1-二甲基环己烷"),
    ("CC1CC(C)CC(C)C1", "1,3,5-trimethylcyclohexane", "1,3,5-三甲基环己烷"),
    ("CCC1CCCCC1C", "1-ethyl-2-methylcyclohexane", "1-乙基-2-甲基环己烷"),
    ("CC1CCCC1C", "1,2-dimethylcyclopentane", "1,2-二甲基环戊烷"),
    ("CCC1CCCCC1CC", "1,2-diethylcyclohexane", "1,2-二乙基环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", polyalkyl_cycloalkane__CASES)
def test_polyalkyl_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_monohalo_cycloalkane.py
# IUPAC: P-61.3.1 / P-14.3.4 / P-22.1.1
# Layer: L2,L3,L5
#
# Monosubstituted monohalo monocycloalkanes: chlorocyclohexane etc.
#
# Parent = saturated monocarbocycle; single halogen on ring carbon.
# Monosubstituted: omit locant (chlorocyclohexane not 1-chlorocyclohexane).
# ==========================================================================
monohalo_cycloalkane__CASES = [
    ("BrC1CCCCC1", "bromocyclohexane", "溴环己烷"),
    ("FC1CCCCC1", "fluorocyclohexane", "氟环己烷"),
    ("IC1CCCCC1", "iodocyclohexane", "碘环己烷"),
    ("ClC1CCCC1", "chlorocyclopentane", "氯环戊烷"),
    ("ClC1CCC1", "chlorocyclobutane", "氯环丁烷"),
    ("ClC1CC1", "chlorocyclopropane", "氯环丙烷"),


    # negative: chain halo and other parents
    ("CCCCCCCl", "1-chlorohexane", "1-氯己烷"),
    ("ClCC", "chloroethane", "氯乙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", monohalo_cycloalkane__CASES)
def test_monohalo_cycloalkane(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_halo_name() -> None:
    r = SMILESNNamer().name("ClC1CCCCC1")
    assert normalize_en(r.en) != "1-chlorohexane"


# ==========================================================================
# 合并自 test_sub_cyclo_fg.py
# IUPAC: P-63.1.1 / P-64.2.1 / P-62.2.1
# Layer: L2,L4,L5
#
# Saturated monocyclic mono-alcohol/ketone/amine with simple ring subs.
#
# Parent = unfused C3–C10 sat carbocycle with exactly one ring-carbon FG
# (OH / ketone carbonyl / primary amine). Ring may carry halo and/or claimable
# n-alkyl (C1–C4) sides. With any ring sub, keep FG locant 1
# (e.g. 2-methylcyclohexan-1-ol); unsubstituted still omits locant
# (cyclohexanol / cyclohexanone / cyclohexanamine).
# ==========================================================================
sub_cyclo_fg__CASES = [
    ("CCC1CCCCC1O", "2-ethylcyclohexan-1-ol", "2-乙基环己-1-醇"),
    ("CC1CCC(O)CC1", "4-methylcyclohexan-1-ol", "4-甲基环己-1-醇"),
    ("CC1CCCC1O", "2-methylcyclopentan-1-ol", "2-甲基环戊-1-醇"),
    # positive: monohalo cycloalcohol
    ("ClC1CCCCC1O", "2-chlorocyclohexan-1-ol", "2-氯环己-1-醇"),
    ("BrC1CCCCC1O", "2-bromocyclohexan-1-ol", "2-溴环己-1-醇"),
    # positive: monoalkyl / monohalo cycloamine
    ("CC1CCCCC1N", "2-methylcyclohexan-1-amine", "2-甲基环己-1-胺"),
    ("CC1CCC(N)CC1", "4-methylcyclohexan-1-amine", "4-甲基环己-1-胺"),
    ("ClC1CCCCC1N", "2-chlorocyclohexan-1-amine", "2-氯环己-1-胺"),
]


@pytest.mark.parametrize("smiles,en,zh", sub_cyclo_fg__CASES)
def test_sub_cyclo_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methyl_cyclohexanol_not_ethane() -> None:
    """Ring monoalcohol + methyl must not collapse to ethane."""
    r = SMILESNNamer().name("CC1CCCCC1O")
    assert r.success
    en = normalize_en(r.en)
    assert en == "2-methylcyclohexan-1-ol"
    assert en != "ethane"
    assert "hexanol" not in en or en.startswith("2-methylcyclo")


# ==========================================================================
# 合并自 test_sub_cycloalkene_fg.py
# IUPAC: P-31.1 / P-22.1.1 / P-63.1.1 / P-64.2.1 / P-14.3.4
# Layer: L2,L4,L5
#
# Alkyl/halo monocyclic monoalkenes + cycloalkenol / cycloalkenone.
#
# 1) Parent cycloalkene with claimable n-alkyl / ring halo:
#    double bond fixed at 1; lowest sub set → 1-methylcyclohexene.
# 2) Ring mono-OH + endocyclic C=C: FG@1, lowest ene → cyclohex-2-en-1-ol.
# 3) Ring mono-ketone + endocyclic C=C: FG@1, lowest ene → cyclohex-2-en-1-one.
# Unsubstituted cyclohexene / sat mono cyclo FG / open alkenol·one stay intact.
# ==========================================================================
sub_cycloalkene_fg__CASES = [
    # positive: alkyl / halo cycloalkene
    ("CC1=CCCCC1", "1-methylcyclohexene", "1-甲基环己烯"),
    ("CCC1=CCCCC1", "1-ethylcyclohexene", "1-乙基环己烯"),
    ("CC1=CCCC1", "1-methylcyclopentene", "1-甲基环戊烯"),
    ("CC1CCC=CC1", "4-methylcyclohexene", "4-甲基环己烯"),
    ("ClC1=CCCCC1", "1-chlorocyclohexene", "1-氯环己烯"),
    ("CC1CCCCC1O", "2-methylcyclohexan-1-ol", "2-甲基环己-1-醇"),
    ("CC=CCO", "but-2-en-1-ol", "丁-2-烯-1-醇"),
    ("CC(=O)C=C", "but-3-en-2-one", "丁-3-烯-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", sub_cycloalkene_fg__CASES)
def test_sub_cycloalkene_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_methylcyclohexene_not_propene() -> None:
    """Alkyl cycloalkene must not collapse to open-chain alkene."""
    r = SMILESNNamer().name("CC1=CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-methylcyclohexene"
    assert en != "propene"
    assert "cyclohexene" in en


# ==========================================================================
# 合并自 test_cycloalkane_exocyclic_fg.py
# IUPAC: P-66.6.1 / P-66.5.1 / P-66.1.1 / P-65.6 / P-65.5
# Layer: L2,L4,L5
#
# Monocyclic cycloalkane with one exocyclic carbonyl-class FG (P-65/P-66).
#
# Parent = unfused saturated monocarbocycle C3–C10 with exactly one ring-attached
# exocyclic FG of type: carbaldehyde / carbonitrile / carboxamide /
# carboxylate (alkyl ester) / carbonyl halide.
# Ring may carry 0–1 simple sub: monohalo or mono C1–C4 n-alkyl; FG attach = 1.
# IUPAC PIN stems: cyclohexanecarbaldehyde / carbonitrile / carboxamide /
# carboxylate / carbonyl chloride (aligned with cycloalkanecarboxylic acid).
# ==========================================================================
cycloalkane_exocyclic_fg__CASES = [
    ("CC(N)=O", "acetamide", "乙酰胺"),


    # positive: ring carries another located prefix → FG locant 1 must be explicit
    ("O=C1CCC(C(=O)[O-])C1", "3-oxocyclopentane-1-carboxylate", "3-氧代环戊烷-1-羧酸根"),
    ("CC1CCCC1C(=O)O", "2-methylcyclopentane-1-carboxylic acid", "2-甲基环戊烷-1-羧酸"),
    ("OC1CCCC1C=O", "2-hydroxycyclopentane-1-carbaldehyde", "2-羟基环戊烷-1-甲醛"),
    # positive: no extra prefix → locant 1 elided (unchanged)
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷羧酸"),
    ("O=CC1CCCCC1", "cyclohexanecarbaldehyde", "环己烷甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkane_exocyclic_fg__CASES)
def test_cycloalkane_exocyclic_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_cycloalkene_exocyclic_fg.py
# IUPAC: P-31.1.2 / P-65.2.2.1 / P-66.6.1.1.3
# Layer: L4,L5
#
# Monocyclic cycloalkene with one exocyclic carbonyl-class FG (acid/ester/aldehyde).
#
# Parent = unfused unsaturated monocarbocycle carrying one exocyclic COOH / COOR / CHO
# (ring C=C must survive in the stem: single ene locant-1 is elided in EN but kept in ZH
# cyclohexene vs 环己-1-烯; non-1 / polyene locants explicit). FG locant becomes explicit
# once the ring is unsaturated because numbering is no longer unique.
# Saturated cycloalkanecarboxylic acid / open-chain / arene must not change.
# ==========================================================================
cycloalkene_exocyclic_fg__CASES = [
    # positive: single endocyclic C=C, ene locant 1 (EN elided), FG@1 explicit
    ("O=C(O)C1=CC(=O)[C@@H](O)[C@H](O)C1",
     "(4S,5R)-4,5-dihydroxy-3-oxocyclohexene-1-carboxylic acid",
     "(4S,5R)-4,5-二羟基-3-氧代环己-1-烯-1-羧酸"),
    # carboxylate anion (去质子酸根)
    ("O=C([O-])C1=CC(=O)[C@@H](O)[C@H](O)C1",
     "(4S,5R)-4,5-dihydroxy-3-oxocyclohexene-1-carboxylate",
     "(4S,5R)-4,5-二羟基-3-氧代环己-1-烯-1-羧酸根"),
    # positive: diene, polyene locants explicit
    ("N[C@@H]1C(C(=O)O)=CC=C[C@@H]1O",
     "(5S,6R)-6-amino-5-hydroxycyclohexa-1,3-diene-1-carboxylic acid",
     "(5S,6R)-6-氨基-5-羟基环己-1,3-二烯-1-羧酸"),
    ("N[C@H]1C(C(=O)O)=CC=C[C@@H]1O",
     "(5S,6S)-6-amino-5-hydroxycyclohexa-1,3-diene-1-carboxylic acid",
     "(5S,6S)-6-氨基-5-羟基环己-1,3-二烯-1-羧酸"),
    ("CC1=CC=CC(O)C1(O)C(=O)O",
     "1,6-dihydroxy-2-methylcyclohexa-2,4-diene-1-carboxylic acid",
     "1,6-二羟基-2-甲基环己-2,4-二烯-1-羧酸"),
    # positive: ester / aldehyde, ene locant non-1 / 1
    ("CCOC(=O)[C@@]1(c2ccccc2)CCC=C[C@@H]1N(C)C",
     "ethyl (1R,2S)-2-(dimethylamino)-1-phenylcyclohex-3-ene-1-carboxylate",
     None),  # zh 二甲氨基 vs 二甲基氨基 为独立胺名规则，不在 R1 范围
    ("C=C(C)C1CC=C(C=O)CC1",
     "4-prop-1-en-2-ylcyclohexene-1-carbaldehyde",
     "4-丙-1-烯-2-基环己-1-烯-1-甲醛"),
    ("OC(=O)C1CCCC1", "cyclopentanecarboxylic acid", "环戊烷羧酸"),
]


@pytest.mark.parametrize("smiles,en,zh", cycloalkene_exocyclic_fg__CASES)
def test_cycloalkene_exocyclic_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
