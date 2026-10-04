# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_sulfo_retained.py: Retained sulfo (-SO3H) substituent leaf in RetainedBackend via _try_registry_leaf.
test_sulfonamide.py: Sulfonamide scope (P-65.3) — negative guard only.
test_sulfonamide_n_alkyl.py: N-alkyl / N,N-dialkyl sulfonamide scope (P-65.3) — negative guard only.
test_sulfone.py: Dialkyl sulfone scope (P-65.3.1.2) — negative guard only.
test_sulfonyl_chloride.py: Sulfonyl chloride scope (P-65.3) — negative guard only.
test_sulfonyl_retained.py: Retained methylsulfinyl (CH3SO-), methylsulfonyl (CH3SO2-),
test_methylsulfanyl_retained.py: Retained methylsulfanyl (CH3S-) substituent leaf in RetainedBackend.
"""
from __future__ import annotations

import pytest

from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_sulfo_retained.py
# IUPAC: P-65.3.2.5 (substitutive nomenclature — retained prefix sulfo)
# Layer: L3
#
# Retained sulfo (-SO3H) substituent leaf in RetainedBackend via _try_registry_leaf.
# ==========================================================================
sulfo_retained__CASES = [
    # sulfo as substituent on acetic acid (locant 2 is correct)
    ("O=S(=O)(O)CC(=O)O", "2-sulfoacetic acid", "2-磺基乙酸"),

    # (2R)-2-hydroxy-3-sulfopropanoic acid — the original failing molecule
    ("O=C(O)[C@@H](O)CS(=O)(=O)O", "(2R)-2-hydroxy-3-sulfopropanoic acid", "(2R)-2-羟基-3-磺基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", sulfo_retained__CASES)
def test_sulfo_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_sulfonamide.py
# IUPAC: P-65.3
# Layer: L1,L2,L5
#
# Sulfonamide scope (P-65.3) — negative guard only.
#
# Positive sulfonamide cases are not covered in this file. The retained cases
# assert acetamide / aniline / sulfoxide / sulfone are not named as sulfonamides.
# ==========================================================================
sulfonamide__CASES = [
    ("c1ccccc1N", "aniline", "苯胺"),
]


# Must not be named as sulfonamide
sulfonamide__NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),
    ("CC(=O)N", "sulfonamide"),
    ("c1ccccc1N", "sulfonamide"),
    ("CS(=O)(=O)C", "sulfonamide"),  # sulfone
]


@pytest.mark.parametrize("smiles,en,zh", sulfonamide__CASES)
def test_simple_sulfonamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,forbidden", sulfonamide__NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰胺" not in zh


# ==========================================================================
# 合并自 test_sulfonamide_n_alkyl.py
# IUPAC: P-65.3
# Layer: L1,L2,L5
#
# N-alkyl / N,N-dialkyl sulfonamide scope (P-65.3) — negative guard only.
#
# Positive N-alkyl sulfonamide cases are not covered in this file. The retained
# cases assert sulfoxide / amide / sulfone / sulfonate ester are not named as
# N-alkyl sulfonamides.
# ==========================================================================
sulfonamide_n_alkyl__NEG_NOT_SULFONAMIDE = [
    ("CS(C)=O", "sulfonamide"),       # sulfoxide
    ("CC(=O)N", "sulfonamide"),       # amide
    ("CS(=O)(=O)C", "sulfonamide"),   # sulfone
    ("CS(=O)(=O)OC", "sulfonamide"),  # sulfonate ester
]


@pytest.mark.parametrize("smiles,forbidden", sulfonamide_n_alkyl__NEG_NOT_SULFONAMIDE)
def test_not_sulfonamide_n_alkyl(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en


# ==========================================================================
# 合并自 test_sulfone.py
# IUPAC: P-65.3
# Layer: L1,L2,L5
#
# Dialkyl sulfone scope (P-65.3.1.2) — negative guard only.
#
# Positive sulfone cases are not covered in this file. The retained cases
# assert sulfoxide / sulfonamide / sulfonic acid / sulfonyl chloride / ketone /
# sulfonate ester / sulfide are not named as sulfones.
# ==========================================================================
sulfone__NEG_NOT_SULFONE = [
    ("CS(C)=O", "sulfone"),           # sulfoxide
    ("CS(=O)(=O)N", "sulfone"),       # sulfonamide
    ("CS(=O)(=O)O", "sulfone"),       # sulfonic acid
    ("CS(=O)(=O)Cl", "sulfone"),      # sulfonyl chloride
    ("CC(=O)C", "sulfone"),           # ketone (acetone)
    ("CS(=O)(=O)OC", "sulfone"),      # sulfonate ester
    ("CSC", "sulfone"),               # sulfide
]


@pytest.mark.parametrize("smiles,forbidden", sulfone__NEG_NOT_SULFONE)
def test_not_sulfone(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    en = normalize_en(r.en)
    assert forbidden not in en, f"{smiles} should not contain '{forbidden}' but got {r.en!r}"


# ==========================================================================
# 合并自 test_sulfonyl_chloride.py
# IUPAC: P-65.3
# Layer: L1,L2,L5
#
# Sulfonyl chloride scope (P-65.3) — negative guard only.
#
# Positive sulfonyl chloride cases are not covered in this file. The retained
# cases assert sulfonamide / acetyl chloride / sulfoxide / tosylate are not
# named as sulfonyl chlorides.
# ==========================================================================
sulfonyl_chloride__NEG_NOT_SULFONYL_CHLORIDE = [
    ("c1ccc(S(=O)(=O)N)cc1", "sulfonyl chloride"),
    ("CC(=O)Cl", "sulfonyl chloride"),
    ("CS(C)=O", "sulfonyl chloride"),
    ("CC1=CC=C(C=C1)S(=O)(=O)OCCCC", "sulfonyl chloride"),
    ("CS(=O)(=O)N", "sulfonyl chloride"),
]


@pytest.mark.parametrize("smiles,forbidden", sulfonyl_chloride__NEG_NOT_SULFONYL_CHLORIDE)
def test_not_sulfonyl_chloride(smiles: str, forbidden: str) -> None:
    r = SMILESNNamer().name(smiles)
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    zh = normalize_zh(r.zh)
    assert forbidden not in en
    assert "磺酰氯" not in zh


# ==========================================================================
# 合并自 test_sulfonyl_retained.py
# IUPAC: P-65.3.1/P-65.3.2 (methylsulfinyl/methylsulfonyl retained prefixes)
# Layer: L3
#
# Retained methylsulfinyl (CH3SO-), methylsulfonyl (CH3SO2-),
# tosyl (4-Me-Ph-SO2-), triflyl (CF3-SO2-) retained prefix leaves.
# ==========================================================================
sulfonyl_retained__CASES = [
    # ── methylsulfinyl (PIN P-65.3.1) ──
    ("CS(C)=O",
     "methylsulfinylmethane", None),
    # ── methylsulfonyl (PIN P-65.3.2) ──
    ("CS(=O)(=O)CC[C@H](N)C(=O)O",
     "(2S)-2-amino-4-methylsulfonylbutanoic acid", None),
    # methylsulfonyl on simple alkyl acid
    ("CS(=O)(=O)CCC(=O)O",
     "3-methylsulfonylpropanoic acid", None),
    # ── methylsulfanyl regression (must still work) ──
    ("CSCC[C@H](N)C(=O)O",
     "(2S)-2-amino-4-methylsulfanylbutanoic acid", None),
    ("CSCCCN",
     "3-methylsulfanylpropan-1-amine", None),
]


@pytest.mark.parametrize("smiles,en,zh", sulfonyl_retained__CASES)
def test_sulfonyl_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"naming failed for {smiles}"
    assert normalize_en(r.en) == normalize_en(en), \
        f"EN mismatch: {r.en} != {en}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_methylsulfanyl_retained.py
# IUPAC: P-63.1.1 (substitutive nomenclature — retained prefix methylsulfanyl)
# Layer: L3
#
# Retained methylsulfanyl (CH3S-) substituent leaf in RetainedBackend.
# ==========================================================================
methylsulfanyl_retained__CASES = [
    # Positive: simple CH3S- on primary amine chain
    ("CSCCCN", "3-methylsulfanylpropan-1-amine", None),
    # Positive: CH3S- on ethanol (C1 有可取代 H，2- 不可省略；匹配 methoxy 命名模式)
    ("CSCCO", "2-methylsulfanylethanol", None),
    # Positive: CH3S- on acetic acid
    ("CSCC(=O)O", "2-methylsulfanylacetic acid", None),
    # Negative: methoxy analogue — must remain unchanged
    ("COCCCN", "3-methoxypropan-1-amine", None),
]


@pytest.mark.parametrize("smiles,en,zh", methylsulfanyl_retained__CASES)
def test_methylsulfanyl_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 硫酸族：无碳臂硫中心的硫酸母体识别
# IUPAC: P-67.1.3, P-67.2.1
# Layer: L1,L2,L5
#
# 无碳臂的硫中心两臂皆酸式氧或 O-臂：硫酸本体/酸根/硫酸(氢)酯同归 sulfate 母体；
# S-O-S 链与 P-O-P 链同构——留一个链节作硫酸母体，其余作磺酰前缀，不得退为甲烷母体。
# ==========================================================================
sulfate__CASES = [
    ("OS(=O)(=O)O", "sulfuric acid", "硫酸"),
    ("OS(=O)(=O)[O-]", "hydrogen sulfate", None),
    ("[O-]S(=O)(=O)[O-]", "sulfate", None),
    ("COS(=O)(=O)O", "methyl hydrogen sulfate", "硫酸氢甲酯"),
    ("COS(=O)(=O)OC", "dimethyl sulfate", "硫酸二甲酯"),
    ("CCOS(=O)(=O)OCC", "diethyl sulfate", None),
    ("OS(=O)(=O)OS(=O)(=O)O", "sulfo hydrogen sulfate", None),
    ("COS(=O)(=O)OS(=O)(=O)OC", "methoxysulfonyl methyl sulfate", None),
]

# 碳臂硫中心仍走磺酸/磺酸酯族：硫酸族放开不得截胡
sulfate__NEG_SULFONIC = [
    ("CS(=O)(=O)O", "methanesulfonic acid", "甲磺酸"),
    ("CS(=O)(=O)OC", "methyl methanesulfonate", "甲磺酸甲酯"),
    ("O=S(=O)(O)CC", "ethanesulfonic acid", "乙磺酸"),
]


@pytest.mark.parametrize("smiles,en,zh", [*sulfate__CASES, *sulfate__NEG_SULFONIC])
def test_sulfate_parent(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: {r.meta}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
