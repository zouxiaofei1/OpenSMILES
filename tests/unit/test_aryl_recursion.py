# 合并自 8 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_aryl_depth2.py: Depth-2 simple leaves on Ph/OPh/Bn arms (alkoxy / nitro / CF3).
test_aryl_depth2_nested.py: Depth-2: propoxy/butoxy leaves; nested unsub Ph on Ph arm; phenol vs chain amine.
test_aryl_depth2_oh_nh2.py: Depth-2 hydroxy / amino leaves on Ph arms; aromatic OH/NH2 not chain poly-FG.
test_aryl_recurse.py: True recursive aryl substituent namer (depth can exceed 2).
test_recursive_aryl.py: Recursive depth-1 aryl: multi-halo Ph/OPh, benzyl, benzyloxy; pyridine parent.
test_aryl_alkenoic.py: Aryl-substituted open-chain alkenoic acids / alkenoates.
test_naphthyl_side.py: Unsubstituted naphthalen-n-yl side chains on open FG parents.
test_pyridinyl_side.py: Unsubstituted pyridin-n-yl side chains on open FG parents.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_aryl_depth2.py
# IUPAC: P-29.3
# Layer: L2,L3
#
# Depth-2 simple leaves on Ph/OPh/Bn arms (alkoxy / nitro / CF3).
#
# Parent remains chain alcohol/acid/amine; substituted phenyl is the recursive arm.
# Halo/methyl depth-1 and bare arene methoxy parents must stay correct.
# ==========================================================================
aryl_depth2__CASES = [
    # positive: methoxy / ethoxy / nitro / CF3 on phenylethanol arm
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("CCOc1ccc(CCO)cc1", "2-(4-ethoxyphenyl)ethanol", "2-(4-乙氧基苯基)乙醇"),
    ("O=[N+]([O-])c1ccc(CCO)cc1", "2-(4-nitrophenyl)ethanol", "2-(4-硝基苯基)乙醇"),
    (
        "FC(F)(F)c1ccc(CCO)cc1",
        "2-[4-(trifluoromethyl)phenyl]ethanol",
        "2-(4-三氟甲基苯基)乙醇",
    ),
    # positive: acid / amine parents
    (
        "O=C(O)Cc1ccc(OC)cc1",
        "2-(4-methoxyphenyl)acetic acid",
        "2-(4-甲氧基苯基)乙酸",
    ),
    (
        "NCCc1ccc(C(F)(F)F)cc1",
        "2-[4-(trifluoromethyl)phenyl]ethanamine",
        "2-(4-三氟甲基苯基)乙胺",
    ),
    # positive: meta methoxy + benzyl-style methanol
    ("COc1cccc(CO)c1", "(3-methoxyphenyl)methanol", "(3-甲氧基苯基)甲醇"),
    # negative: depth-1 halo arm + bare arene alkoxy (not recursive)
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
    ("COc1ccc(CC)cc1", "1-ethyl-4-methoxybenzene", "1-乙基-4-甲氧基苯"),
]


# ==========================================================================
# 合并自 test_aryl_depth2_nested.py
# IUPAC: P-29.3
# Layer: L2,L3
#
# Depth-2: propoxy/butoxy leaves; nested unsub Ph on Ph arm; phenol vs chain amine.
#
# Parent stays chain alcohol/amine; nested Ph is a recursive leaf (not Ph-on-Ph arm rejection).
# Bare arene propoxy and biphenyl-as-benzene must stay correct.
# ==========================================================================
aryl_depth2_nested__CASES = [
    # propoxy / butoxy leaves on phenylethanol
    ("CCCOc1ccc(CCO)cc1", "2-(4-propoxyphenyl)ethanol", "2-(4-丙氧基苯基)乙醇"),
    ("CCCCOc1ccc(CCO)cc1", "2-(4-butoxyphenyl)ethanol", "2-(4-丁氧基苯基)乙醇"),
    # nested unsubstituted phenyl on arm Ph
    (
        "c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-(4-phenylphenyl)ethanol",
        "2-(4-苯基苯基)乙醇",
    ),
    (
        "c1ccc(-c2cccc(CCO)c2)cc1",
        "2-(3-phenylphenyl)ethanol",
        "2-(3-苯基苯基)乙醇",
    ),
    # phenolic OH must not steal aliphatic amine parent
    (
        "Oc1ccc(CCN)cc1",
        "2-(4-hydroxyphenyl)ethanamine",
        "2-(4-羟基苯基)乙胺",
    ),
    # negatives: arene propoxy parent; biphenyl as benzene+phenyl; prior depth-2
    ("CCCOc1ccccc1", "propoxybenzene", "丙氧基苯"),
    ("CCCOc1ccc(CC)cc1", "1-ethyl-4-propoxybenzene", "1-乙基-4-丙氧基苯"),
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("Oc1ccc(CCO)cc1", "2-(4-hydroxyphenyl)ethanol", "2-(4-羟基苯基)乙醇"),
]


# ==========================================================================
# 合并自 test_aryl_depth2_oh_nh2.py
# IUPAC: P-29.3
# Layer: L2,L3
#
# Depth-2 hydroxy / amino leaves on Ph arms; aromatic OH/NH2 not chain poly-FG.
#
# Chain alcohol/amine/acid stays parent; phenolic OH and aniline-NH2 are arm leaves.
# True alkanediol / diamine / phenol / aniline must remain correct.
# ==========================================================================
aryl_depth2_oh_nh2__CASES = [
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


# ==========================================================================
# 合并自 test_aryl_recurse.py
# IUPAC: P-29.3 / P-14.3.1
# Layer: L2,L3
#
# True recursive aryl substituent namer (depth can exceed 2).
#
# Nested Ph / OPh on a Ph arm may themselves carry simple leaves (halo/Me/
# alkoxy/nitro/CF3/OH/NH2) or further nested unsub Ph within max depth.
# Bare biphenyl and depth-1/2 arms must stay correct.
# ==========================================================================
aryl_recurse__CASES = [
    # depth-3: nested Ph with simple leaves
    (
        "Clc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorophenyl)phenyl]ethanol",
        "2-[4-(4-氯苯基)苯基]乙醇",
    ),
    (
        "Cc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-methylphenyl)phenyl]ethanol",
        "2-[4-(4-甲基苯基)苯基]乙醇",
    ),
    (
        "COc1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-methoxyphenyl)phenyl]ethanol",
        "2-[4-(4-甲氧基苯基)苯基]乙醇",
    ),
    (
        "O=[N+]([O-])c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-[4-(4-nitrophenyl)phenyl]ethanol",
        "2-[4-(4-硝基苯基)苯基]乙醇",
    ),
    (
        "Clc1cc(Cl)ccc1-c2ccc(CCO)cc2",
        "2-[4-(2,4-dichlorophenyl)phenyl]ethanol",
        "2-[4-(2,4-二氯苯基)苯基]乙醇",
    ),
    # nested phenoxy (unsub / mono-halo)
    (
        "c1ccc(Oc2ccc(CCO)cc2)cc1",
        "2-(4-phenoxyphenyl)ethanol",
        "2-(4-苯氧基苯基)乙醇",
    ),
    (
        "Clc1ccc(Oc2ccc(CCO)cc2)cc1",
        "2-[4-(4-chlorophenoxy)phenyl]ethanol",
        "2-[4-(4-氯苯氧基)苯基]乙醇",
    ),
    # prior depth-2 unsub nested Ph
    (
        "c1ccc(-c2ccc(CCO)cc2)cc1",
        "2-(4-phenylphenyl)ethanol",
        "2-(4-苯基苯基)乙醇",
    ),
    # negatives
    ("c1ccc(-c2ccccc2)cc1", "phenylbenzene", "苯基苯"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("COc1ccc(CCO)cc1", "2-(4-methoxyphenyl)ethanol", "2-(4-甲氧基苯基)乙醇"),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
]


# ==========================================================================
# 合并自 test_recursive_aryl.py
# IUPAC: P-29.3 / P-14.3.1 / P-22.2.1
# Layer: L2,L3
#
# Recursive depth-1 aryl: multi-halo Ph/OPh, benzyl, benzyloxy; pyridine parent.
#
# - Ph/OPh: up to 3 terminal ring-halos; lowest locants from attach=1 (2- not 6-).
# - benzyl = parent–CH2–Ph; benzyloxy = parent–O–CH2–Ph (not ethoxy).
# - pyridine (unfused C5N) may carry phenyl/phenoxy/benzyl/benzyloxy.
# Negatives: anisole / ethoxybenzene / phenoxybenzene must stay correct.
# ==========================================================================
recursive_aryl__CASES = [
    # multi-halo phenyl as substituent on FG parent (attach=1 → lowest set)
    (
        "O=Cc1ccc(-c2ccc(Cl)c(Cl)c2)cc1",
        "4-(3,4-dichlorophenyl)benzaldehyde",
        "4-(3,4-二氯苯基)苯甲醛",
    ),
    # ortho-halo phenyl lowest locant is 2- (not 6-) when Ph is the arm
    (
        "O=Cc1ccc(-c2ccccc2Cl)cc1",
        "4-(2-chlorophenyl)benzaldehyde",
        "4-(2-氯苯基)苯甲醛",
    ),
    # benzylbenzene
    ("c1ccc(Cc2ccccc2)cc1", "benzylbenzene", "苄基苯"),
    # pyridine + phenyl
    ("n1ccc(-c2ccccc2)cc1", "4-phenylpyridine", "4-苯基吡啶"),
    # pyridine + benzyl
    ("n1ccc(Cc2ccccc2)cc1", "4-benzylpyridine", "4-苄基吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", recursive_aryl__CASES)
def test_recursive_aryl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_phenylpyridine_not_nonane() -> None:
    r = SMILESNNamer().name("c1ccc(-c2ccncc2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "nonane" not in en
    assert "phenylpyridine" in en


def test_dichloro_n_phenyl_not_unsub() -> None:
    """Multi-halo N-Ph must not be treated as unsubstituted N-phenylamide."""
    r = SMILESNNamer().name("O=C(Nc1cc(Cl)cc(Cl)c1)C")
    assert r.success
    en = normalize_en(r.en)
    assert "N-phenylacetamide" not in en
    assert "N-phenyl" not in en


def test_mixed_halo_phenyl_alpha_order() -> None:
    """Mixed halos: cite by English name order (bromo before chloro); the
    alpha-earlier substituent gets the lower locant (P-14.4)."""
    r = SMILESNNamer().name("O=Cc1ccc(-c2cc(Cl)cc(Br)c2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "3-bromo-5-chlorophenyl" in en
    assert "5-bromo-3-chlorophenyl" not in en


# ==========================================================================
# 合并自 test_aryl_alkenoic.py
# IUPAC: P-31.1 / P-65.1.1
# Layer: L2
#
# Aryl-substituted open-chain alkenoic acids / alkenoates.
#
# IUPAC P-31.1 / P-65.1.1 / P-65.6: the parent is the open-chain unsaturation +
# principal FG (COOH/ester). An aromatic ring is a substituent, not a reason to
# reject alkenoic/alkenoate selection. Gate is parent atoms not in ring, not
# molecule-level has_ring.
# ==========================================================================
aryl_alkenoic__CASES = [
    # positive: aryl-substituted open-chain alkenoic acids / alkenoates
    (
        r"OC(=O)/C=C/C1=CC=CC=C1",
        "(E)-3-phenylprop-2-enoic acid",
        "(E)-3-苯基丙-2-烯酸",
    ),
    (
        r"O=C(O)/C=C/c1ccc(C(F)(F)F)cc1",
        "(E)-3-[4-(trifluoromethyl)phenyl]prop-2-enoic acid",
        "(E)-3-(4-三氟甲基苯基)丙-2-烯酸",
    ),
    (
        r"COC(=O)/C=C/c1ccc(O)c(OC)c1",
        "methyl (E)-3-(4-hydroxy-3-methoxyphenyl)prop-2-enoate",
        "(E)-3-(4-羟基-3-甲氧基苯基)丙-2-烯酸甲酯",
    ),
    (
        r"COC(=O)/C=C/c1ccc(OC)c(OC)c1",
        "methyl (E)-3-(3,4-dimethoxyphenyl)prop-2-enoate",
        "(E)-3-(3,4-二甲氧基苯基)丙-2-烯酸甲酯",
    ),
    # positive: no-ring control remains alkenoic
    (r"O=C(O)/C=C/C", "(E)-but-2-enoic acid", "(E)-丁-2-烯酸"),
    # negative: saturated / diacid / ring-acid must not become enoic
    ("O=C(O)CCc1ccccc1", "3-phenylpropanoic acid", "3-苯基丙酸"),
    ("O=C(O)CC(=O)O", "propanedioic acid", "丙二酸"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


# ==========================================================================
# 合并自 test_naphthyl_side.py
# IUPAC: P-29.3 / P-25
# Layer: L2,L3
#
# Unsubstituted naphthalen-n-yl side chains on open FG parents.
# ==========================================================================
naphthyl_side__CASES = [
    (
        "OCCc1ccc2ccccc2c1",
        "2-(naphthalen-2-yl)ethanol",
        "2-(萘-2-基)乙醇",
    ),
    (
        "OCCc1cccc2ccccc12",
        "2-(naphthalen-1-yl)ethanol",
        "2-(萘-1-基)乙醇",
    ),
    (
        "N#CCc1ccc2ccccc2c1",
        "2-(naphthalen-2-yl)acetonitrile",
        "2-(萘-2-基)乙腈",
    ),
    (
        "NCCc1cccc2ccccc12",
        "2-(naphthalen-1-yl)ethanamine",
        "2-(萘-1-基)乙胺",
    ),
    (
        "O=C(C)c1ccc2ccccc2c1",
        "1-(naphthalen-2-yl)ethanone",
        "1-(萘-2-基)乙酮",
    ),
    # negatives
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("Cc1ccc2ccccc2c1", "2-methylnaphthalene", "2-甲基萘"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("OCCc1ccncc1", "2-(pyridin-4-yl)ethanol", "2-(吡啶-4-基)乙醇"),
]


# ==========================================================================
# 合并自 test_pyridinyl_side.py
# IUPAC: P-29.3 / P-44
# Layer: L2,L3
#
# Unsubstituted pyridin-n-yl side chains on open FG parents.
#
# Chain alcohol/amine/nitrile/ketone claim pyridin-2/3/4-yl; pyridine as
# parent with alkyl sides unchanged.
# ==========================================================================
pyridinyl_side__CASES = [
    (
        "OCCc1ccncc1",
        "2-(pyridin-4-yl)ethanol",
        "2-(吡啶-4-基)乙醇",
    ),
    (
        "NCCc1ccncc1",
        "2-(pyridin-4-yl)ethanamine",
        "2-(吡啶-4-基)乙胺",
    ),
    (
        "N#CCc1ccncc1",
        "2-(pyridin-4-yl)acetonitrile",
        "2-(吡啶-4-基)乙腈",
    ),
    (
        "N#CCc1ccccn1",
        "2-(pyridin-2-yl)acetonitrile",
        "2-(吡啶-2-基)乙腈",
    ),
    (
        "O=C(C)c1ccncc1",
        "1-(pyridin-4-yl)ethanone",
        "1-(吡啶-4-基)乙酮",
    ),
    (
        "OCCc1ccccn1",
        "2-(pyridin-2-yl)ethanol",
        "2-(吡啶-2-基)乙醇",
    ),
    # simple leaves on pyridine
    (
        "BrCC(=O)C1=NC=CC(=C1)OC",
        "2-bromo-1-(4-methoxypyridin-2-yl)ethanone",
        "2-溴-1-(4-甲氧基吡啶-2-基)乙酮",
    ),
    (
        "O=C(C)c1nc(OC)ccc1",
        "1-(6-methoxypyridin-2-yl)ethanone",
        "1-(6-甲氧基吡啶-2-基)乙酮",
    ),
    (
        "N#CCc1cc(Cl)ccn1",
        "2-(4-chloropyridin-2-yl)acetonitrile",
        "2-(4-氯吡啶-2-基)乙腈",
    ),
    # negatives: pyridine parent, phenyl ethanol
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Cc1ccncc1", "4-methylpyridine", "4-甲基吡啶"),
    ("Clc1ccc(CCO)cc1", "2-(4-chlorophenyl)ethanol", "2-(4-氯苯基)乙醇"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine", "4-苯基吡啶"),
]
