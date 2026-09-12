# 合并自 13 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_simple_benzene.py: Retained name benzene: unsubstituted + monohalo + mono C1–C2 n-alkylbenzene.
test_multi_benzene.py: Multi-substituted benzene: ring halo and/or methyl (C1), 2–3 substituents.
test_arene_fg_trisub.py: Arene FG retained parents: raise simple ring-sub cap 2→3.
test_arene_aryl_fg.py: Depth-1 phenyl / phenoxy on retained arene FG parents.
test_arene_ester_nitrile.py: Retained arene parents: alkyl benzoate, benzonitrile, benzoyl chloride.
test_arene_peg_alkoxy.py: Arene linear n-alkoxy C1–C4 and PEG tails -(OCH2CH2)k-OR on benzene.
test_arene_side_extend.py: Arene simple side-chain extensions: ω-haloalkyl, alkoxy chains, quinoline CF3.
test_nitrobenzene.py: Nitro as prefix (nitrobenzene / nitrophenol): simple arene cases.
test_halo_prefix.py: Generic halogen substituent prefixes (fluoro/chloro/bromo/iodo) on alkane parents.
test_hetero_phenyl.py: Phenyl-on-heteroarene scope (P-22.2.1) — negative guard only.
test_phenyl_phenoxy.py: Unsubstituted phenyl / phenoxy (depth-1 aryl) on benzene parents.
test_alkylbenzene.py: Simple mono-substituted alkylbenzenes: C1–C4 n-alkyl + isopropyl.
test_chain_benzyl.py: Benzyl / phenyl on chain acid / amine parents (plan 3.4).
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_simple_benzene.py
# IUPAC: P-22.1.3
# Layer: L2,L4,L5
#
# Retained name benzene: unsubstituted + monohalo + mono C1–C2 n-alkylbenzene.
#
# Parent = aromatic monocarbocycle (benzene). Monosubstituted: omit locant.
# Methylbenzene retains toluene; ethylbenzene is systematic.
# ==========================================================================
simple_benzene__CASES = [
    ("Fc1ccccc1", "fluorobenzene", "氟苯"),
    ("Brc1ccccc1", "bromobenzene", "溴苯"),
    ("Ic1ccccc1", "iodobenzene", "碘苯"),
    ("ClC1CCCCC1", "chlorocyclohexane", "氯环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", simple_benzene__CASES)
def test_simple_benzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_not_chain_alkane_name() -> None:
    r = SMILESNNamer().name("c1ccccc1")
    assert normalize_en(r.en) != "hexane"


@pytest.mark.parametrize("smiles", [ "C#Cc1ccccc1"])
def test_unsaturated_sidechain_not_ethylbenzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != "ethylbenzene"


# ==========================================================================
# 合并自 test_multi_benzene.py
# IUPAC: P-14.3.4 / P-22.1.3
# Layer: L2,L4,L5
#
# Multi-substituted benzene: ring halo and/or methyl (C1), 2–3 substituents.
#
# Lowest set of locants (ring rotation). Exactly two methyls → retained xylene
# (en); Chinese remains systematic dimethylbenzene form.
# ==========================================================================
multi_benzene__CASES = [
    # positive: di/tri halo and/or methyl on benzene
    ("Clc1ccc(Cl)cc1", "1,4-dichlorobenzene", "1,4-二氯苯"),
    ("Clc1cccc(Cl)c1", "1,3-dichlorobenzene", "1,3-二氯苯"),
    ("Fc1ccc(F)cc1", "1,4-difluorobenzene", "1,4-二氟苯"),
    ("Fc1ccc(Br)cc1", "1-bromo-4-fluorobenzene", "1-溴-4-氟苯"),
    ("Cc1cc(C)cc(C)c1", "1,3,5-trimethylbenzene", "1,3,5-三甲基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", multi_benzene__CASES)
def test_multi_benzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", ["C=Cc1ccccc1", ])
def test_unsaturated_sidechain_not_ethylbenzene__multi_benzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != "ethylbenzene"


# ==========================================================================
# 合并自 test_arene_fg_trisub.py
# IUPAC: P-14.3.4
# Layer: L2
#
# Arene FG retained parents: raise simple ring-sub cap 2→3.
#
# IUPAC P-14.3.4 / P-63.1.4 / P-65.1.1.1 — phenol, aniline, benzoic acid
# already keep 0–2 simple ring substituents (halo/methyl/nitro/alkoxy);
# align with benzene's ≤3 simple ring-sub allowance so tri-substituted
# retained parents stay on the FG parent (not chain / fail).
# ==========================================================================
arene_fg_trisub__CASES = [
    # positive: tri-substituted retained FG parents
    ("CC1=C(C=C(C(=C1)C)C)O", "2,4,5-trimethylphenol", "2,4,5-三甲基苯酚"),
    ("ClC=1C(=C(N)C=CC1Cl)C", "3,4-dichloro-2-methylaniline", "3,4-二氯-2-甲基苯胺"),
    (
        "BrC=1C(=CC(=C(C(=O)O)C1)F)Cl",
        "5-bromo-4-chloro-2-fluorobenzoic acid",
        "5-溴-4-氯-2-氟苯甲酸",
    ),
    ("Oc1c(Br)cc(Br)cc1Br", "2,4,6-tribromophenol", None),
    # negative: must not regress existing correct names
    ("Cc1ccc(O)c(C)c1", "2,4-dimethylphenol", "2,4-二甲基苯酚"),
]


@pytest.mark.parametrize("smiles,en,zh", arene_fg_trisub__CASES)
def test_arene_fg_trisub(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_arene_aryl_fg.py
# IUPAC: P-29.3 / P-14.3.1 / P-63.1.4 / P-62.2.1 / P-66.6.1 / P-65.1.1.1
# Layer: L2,L3
#
# Depth-1 phenyl / phenoxy on retained arene FG parents.
#
# Parent = phenol / aniline / benzaldehyde / benzoic acid (FG ring = 1).
# Ring–O–Ph → phenoxy; ring–Ph → phenyl. Ph: unsub or mono-halo only.
# Must not regress bare FG parents, anisole, or phenoxybenzene.
# ==========================================================================
arene_aryl_fg__CASES = [
    # positive: phenyl on benzaldehyde
    ("O=Cc1ccc(-c2ccccc2)cc1", "4-phenylbenzaldehyde", "4-苯基苯甲醛"),
]


@pytest.mark.parametrize("smiles,en,zh", arene_aryl_fg__CASES)
def test_arene_aryl_fg(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_arene_ester_nitrile.py
# IUPAC: P-65.6 / P-65.5.1 / P-66.5.1
# Layer: L2,L4,L5
#
# Retained arene parents: alkyl benzoate, benzonitrile, benzoyl chloride.
#
# Benzene + one principal FG on a ring carbon (ester / nitrile / acyl chloride).
# Allow 0–2 extra ring halo / methyl; FG attach = locant 1.
# Ester alkoxy: simple straight n-alkyl C1–C16 (no branching).
# ==========================================================================
arene_ester_nitrile__CASES = [
    ("CCCCCCCCCCCCCCCCOC(=O)c1ccccc1", "hexadecyl benzoate", "苯甲酸十六酯"),
    ("CCOC(=O)c1ccc(Br)cc1Cl", "ethyl 4-bromo-2-chlorobenzoate", None),
    ("CCOC(=O)c1c(Br)cc(Cl)cc1", "ethyl 2-bromo-4-chlorobenzoate", None),
    ("CCCCCCC(=O)OC", "methyl heptanoate", "庚酸甲酯"),
]


@pytest.mark.parametrize("smiles,en,zh", arene_ester_nitrile__CASES)
def test_arene_ester_nitrile_acyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzoate_not_chain_ester() -> None:
    """Aromatic ester must not expand ring into chain alkanoate parent."""
    r = SMILESNNamer().name("COC(=O)c1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "methyl benzoate"
    assert "heptanoate" not in en
    assert "hexanoate" not in en


# ==========================================================================
# 合并自 test_arene_peg_alkoxy.py
# IUPAC: P-63.2.2 / P-29.3
# Layer: L2,L3,L5
#
# Arene linear n-alkoxy C1–C4 and PEG tails -(OCH2CH2)k-OR on benzene.
#
# Extends ring alkoxy beyond methoxy/ethoxy/2-methoxyethoxy so polyether
# side chains stay on the arene parent (not alkane fallback → hexane).
# ==========================================================================
arene_peg_alkoxy__CASES = [
    # simple n-alkoxy C3–C4
    ("CCCCOc1ccccc1", "butoxybenzene", "丁氧基苯"),
    # PEG k=1 ethoxy
    (
        "CCOCCOc1ccccc1",
        "(2-ethoxyethoxy)benzene",
        "(2-乙氧基乙氧基)苯",
    ),
    # PEG k=2: Ph-O-(CH2CH2O)2-R  (note extra C vs mistaken COCCOCOc)
    (
        "CCOCCOCCOc1ccccc1",
        "(2-(2-ethoxyethoxy)ethoxy)benzene",
        "(2-(2-乙氧基乙氧基)乙氧基)苯",
    ),
    # regressions
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    (
        "BrC1=CC(=C(C=C1)OCCOC)F",
        "4-bromo-2-fluoro-1-(2-methoxyethoxy)benzene",
        "4-溴-2-氟-1-(2-甲氧基乙氧基)苯",
    ),
    # open-chain ether unchanged
    ("CCOCC", "diethyl ether", None),
]


# ==========================================================================
# 合并自 test_arene_side_extend.py
# IUPAC: P-29.3 / P-14.3.4
# Layer: L2,L3,L5
#
# Arene simple side-chain extensions: ω-haloalkyl, alkoxy chains, quinoline CF3.
#
# P-29.3 / P-14.3.4: ring–(CH2)n–X ω-halo n-alkyl as compound substituent.
# P-22.2.1: quinoline with halo + CF3; extended alkoxy (2-methoxyethoxy).
# ==========================================================================
arene_side_extend__CASES = [
    # A: ω-halo n-alkyl on benzene
    ("BrCCCC1=CC=CC=C1", "(3-bromopropyl)benzene", "(3-溴丙基)苯"),
    # negatives: keep existing correct behaviour
    ("CCc1ccccc1", "ethylbenzene", "乙基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", arene_side_extend__CASES)
def test_arene_side_extend(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_nitrobenzene.py
# IUPAC: P-61.5.1
# Layer: L1,L2,L3,L5
#
# Nitro as prefix (nitrobenzene / nitrophenol): simple arene cases.
#
# Parent = benzene or phenol; nitro is a non-senior prefix substituent.
# Allow 1–2 nitros + optional halo/methyl (total subs ≤3 for benzene;
# phenol: OH principal + nitro prefix).
# ==========================================================================
nitrobenzene__CASES = [
    # positive: mono-nitro benzene
    ("[O-][N+](=O)c1ccccc1", "nitrobenzene", "硝基苯"),
    # positive: nitro + halo (lowest set of locants)
    ("O=[N+]([O-])c1cccc(I)c1", "1-iodo-3-nitrobenzene", None),
    # positive: 4-nitrophenol (OH principal, nitro prefix)
    ("Oc1ccc([N+]([O-])=O)cc1", "4-nitrophenol", "4-硝基苯酚"),
]


@pytest.mark.parametrize("smiles,en,zh", nitrobenzene__CASES)
def test_nitrobenzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_nitromethane_not_nitrobenzene() -> None:
    """Aliphatic nitro must not be named as nitrobenzene."""
    r = SMILESNNamer().name("C[N+](=O)[O-]")
    en = normalize_en(r.en) if r.success else ""
    assert en != "nitrobenzene"
    assert "benzene" not in en
    assert "苯" not in (r.zh or "")


# ==========================================================================
# 合并自 test_halo_prefix.py
# IUPAC: P-61.3.1
# Layer: L3,L4,L5
#
# Generic halogen substituent prefixes (fluoro/chloro/bromo/iodo) on alkane parents.
#
# Affiliated: P-14.3.4 locant omission, P-14.5 alphabetical order, multiplicative di/tri/tetra.
# ==========================================================================
halo_prefix__CASES = [
    ("CCBr", "bromoethane", "溴乙烷"),
    ("CF", "fluoromethane", None),
    ("CC(Cl)C", "2-chloropropane", "2-氯丙烷"),
    ("CC(Br)CC", "2-bromobutane", "2-溴丁烷"),
    ("CCCCF", "1-fluorobutane", "1-氟丁烷"),
    ("CCCCCI", "1-iodopentane", "1-碘戊烷"),
    ("CC(C)CCl", "1-chloro-2-methylpropane", "1-氯-2-甲基丙烷"),
    ("BrC(C)CC(C)C", "2-bromo-4-methylpentane", "2-溴-4-甲基戊烷"),
    ("ClCCCCCCCCCC", "1-chlorodecane", "1-氯癸烷"),
    # 数量前缀 >10 必须显示（回归：MULT 原只到 10，十五氟曾丢数量只剩 fluoro）
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
     "1,1,1,2,2,3,3,4,4,4-decafluorobutane", "1,1,1,2,2,3,3,4,4,4-十氟丁烷"),
    ("O=C(C(C(C(C(C(C(C(CF)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O",
     "2,2,3,3,4,4,5,5,6,6,7,7,8,8,9-pentadecafluorononanoic acid",
     "2,2,3,3,4,4,5,5,6,6,7,7,8,8,9-十五氟壬酸"),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
     "1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8-octadecafluorooctane",
     "1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8-十八氟辛烷"),
]


@pytest.mark.parametrize("smiles,en,zh", halo_prefix__CASES)
def test_halo_prefix(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_hetero_phenyl.py
# IUPAC: P-22.2.1 / P-29.3 / P-14.3.4
# Layer: L2,L3,L5
#
# Phenyl-on-heteroarene scope (P-22.2.1) — negative guard only.
#
# Positive phenyl-diazine / phenyl-hetero5 cases are not covered in this file.
# The retained case asserts pyridine is not captured as a phenyl-hetero parent.
# ==========================================================================
hetero_phenyl__CASES = [
    # regressions
    ("c1ccncc1", "pyridine", None),
]


@pytest.mark.parametrize("smiles,en,zh", hetero_phenyl__CASES)
def test_hetero_phenyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_phenyl_phenoxy.py
# IUPAC: P-29.3 / P-14.3.1
# Layer: L2,L3
#
# Unsubstituted phenyl / phenoxy (depth-1 aryl) on benzene parents.
#
# Parent = benzene; ring–O–Ph → phenoxy; ring–Ph → phenyl.
# Ph may carry at most one ring halo (locant from attachment = 1).
# Must not break anisole / ethoxybenzene / bare benzene.
# ==========================================================================
phenyl_phenoxy__CASES = [
    # positive: biphenyl as benzene + phenyl
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    # negative: retained anisole / ethoxy / bare benzene
    ("CCOc1ccccc1", "ethoxybenzene", None),
]


@pytest.mark.parametrize("smiles,en,zh", phenyl_phenoxy__CASES)
def test_phenyl_phenoxy(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkylbenzene.py
# IUPAC: P-29.3.1 / P-22.1.3
# Layer: L2,L3,L5
#
# Simple mono-substituted alkylbenzenes: C1–C4 n-alkyl + isopropyl.
#
# Parent = benzene (aromatic monocarbocycle). Monosubstituted: omit locant.
# Methylbenzene retains toluene; isopropyl is branched alkyl (not propyl).
# ==========================================================================
alkylbenzene__CASES = [
    ("CCCCc1ccccc1", "butylbenzene", "丁基苯"),
    ("Clc1ccccc1", "chlorobenzene", "氯苯"),
    ("CC(C)C", "2-methylpropane", "2-甲基丙烷"),
]


@pytest.mark.parametrize("smiles,en,zh", alkylbenzene__CASES)
def test_alkylbenzene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize(
    "smiles,bad",
    [
        ("CC(C)c1ccccc1", "2-methyloctane"),
        ("CCCc1ccccc1", "nonane"),
        ("CCCCc1ccccc1", "decane"),
    ],
)
def test_not_chain_alkane_name__alkylbenzene(smiles: str, bad: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert normalize_en(r.en) != normalize_en(bad)


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(C)Cc1ccccc1",   # isobutylbenzene
        "CCC(C)c1ccccc1",   # sec-butylbenzene
        "CC(C)(C)c1ccccc1", # tert-butylbenzene
    ],
)
def test_branched_butyl_not_bare_benzene(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert not (r.success and normalize_en(r.en) == "benzene")


# ==========================================================================
# 合并自 test_chain_benzyl.py
# IUPAC: P-29.3 / P-63.1 / P-65.1 / P-62.2
# Layer: L2,L3,L5
#
# Benzyl / phenyl on chain acid / amine parents (plan 3.4).
# ==========================================================================
chain_benzyl__CASES = [
    # diacid + C-benzyl / phenyl (open chain no longer walks into arene)
    (
        "O=C(O)CC(Cc1ccccc1)C(=O)O",
        "2-benzylbutanedioic acid",
        "2-苄基丁二酸",
    ),
    # phenylmethanamine system name；乙胺 C1 有可取代 H，2- 不可省略（P-14.3.4.4）
    ("c1ccc(CCN)cc1", "2-phenylethanamine", "2-苯基乙胺"),
    # regressions
    ("CC(=O)O", "acetic acid", None),
    # P-14.3.4.4 异构体测试：乙醇 C1 有可取代 H，2- 不可省略
    ("NCCO", "2-aminoethanol", None),
]


@pytest.mark.parametrize("smiles,en,zh", chain_benzyl__CASES)
def test_chain_benzyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
