# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_mono_ester.py: Simple acyclic saturated monoesters (alkyl alkanoates).
test_alkenoate.py: Open-chain monounsaturated monoesters (alkenoates).
test_alkanedioate_diester.py: Open-chain symmetric dialkyl alkanedioates (diesters of diacids).
test_alkenedioate_diester.py: Open-chain symmetric dialkyl alkenedioates (unsaturated diesters).
test_ester_alkoxy_sides.py: Ester alkoxy sides: benzyl / tert-butyl / isopropyl / phenyl (P-65.6).
test_complex_benzoate.py: Complex O-alkyl benzoate parent hit (Ph–C(=O)–O–R).
test_principal_benzoate.py:
"""
from __future__ import annotations

import pytest

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.kind_registry import pack_parent_stem
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer2.principal_parent import rule_driven_parent_candidates
from namepredict.namer import SMILESNNamer, try_candidate
from namepredict.tools.re import normalize_en, normalize_zh
from rdkit import Chem

# ==========================================================================
# 合并自 test_mono_ester.py
# IUPAC: P-65.6
# Layer: L1,L2,L4,L5
#
# Simple acyclic saturated monoesters (alkyl alkanoates).
#
# Ester R–C(=O)–O–R': carbonyl C has =O and single-bond O with no H
# (alkoxy oxygen, not OH). Functional class: alkyl alkanoate / 酸+烷词干+酯.
# ==========================================================================
mono_ester__CASES = [
    ("COC(C)=O", "methyl acetate", "乙酸甲酯"),
    ("CCOC(=O)CC", "ethyl propanoate", "丙酸乙酯"),
    ("CCCC(=O)OC", "methyl butanoate", "丁酸甲酯"),
    ("CC(=O)OCCC", "propyl acetate", "乙酸丙酯"),
    ("CCOC(=O)CCC", "ethyl butanoate", "丁酸乙酯"),
    ("CCCCCC(=O)OC", "methyl hexanoate", "己酸甲酯"),
]


@pytest.mark.parametrize("smiles,en,zh", mono_ester__CASES)
def test_mono_ester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenoate.py
# IUPAC: P-65.6 / P-31.1
# Layer: L2,L4,L5
#
# Open-chain monounsaturated monoesters (alkenoates).
#
# Ester is the principal characteristic group (carbonyl C = locant 1); one
# non-aromatic C=C on the acyl chain is expressed as -n-enoate / -n-烯酸…酯.
# Alkoxy limited to unsubstituted C1–C4; no (E)/(Z) this cycle.
# ==========================================================================
alkenoate__CASES = [
    ("C=CC(=O)OCC", "ethyl prop-2-enoate", "丙-2-烯酸乙酯"),
    ("CC=CC(=O)OC", "methyl but-2-enoate", "丁-2-烯酸甲酯"),
    ("C=CCC(=O)OC", "methyl but-3-enoate", "丁-3-烯酸甲酯"),
    ("C=CC(=O)OCCC", "propyl prop-2-enoate", "丙-2-烯酸丙酯"),
    ("CCC=CC(=O)OCC", "ethyl pent-2-enoate", "戊-2-烯酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenoate__CASES)
def test_alkenoate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkanedioate_diester.py
# IUPAC: P-65.1.1 / P-65.6
# Layer: L2,L5
#
# Open-chain symmetric dialkyl alkanedioates (diesters of diacids).
#
# Symmetric saturated diesters: di{alkyl} {alkane}dioate / {二酸}二{烷}酯.
# C2 retained oxalate/草酸; C3+ systematic propanedioate etc.
# ==========================================================================
alkanedioate_diester__CASES = [
    ("O=C(O)C(=O)O", "oxalic acid", "草酸"),
]


@pytest.mark.parametrize("smiles,en,zh", alkanedioate_diester__CASES)
def test_alkanedioate_diester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_alkenedioate_diester.py
# IUPAC: P-65.1.1 / P-31.1 / P-93
# Layer: L2,L4,L5
#
# Open-chain symmetric dialkyl alkenedioates (unsaturated diesters).
#
# Mono-ene diesters of open-chain diacids: di{alkyl} (E/Z)-alk-n-enedioate /
# (E/Z)-{烷}-n-烯二酸二{烷}酯. Saturated diesters must not regress.
# ==========================================================================
alkenedioate_diester__CASES = [
    ("OC(=O)CC(=O)O", "propanedioic acid", "丙二酸"),
]


@pytest.mark.parametrize("smiles,en,zh", alkenedioate_diester__CASES)
def test_alkenedioate_diester(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_ester_alkoxy_sides.py
# IUPAC: P-65.6
# Layer: L2,L5
#
# Ester alkoxy sides: benzyl / tert-butyl / isopropyl / phenyl (P-65.6).
#
# O-CH2Ph is benzyl not methyl; O-C(CH3)3 is tert-butyl not ethyl.
# ==========================================================================
ester_alkoxy_sides__CASES = [
    # positive: special alkoxy
    ("CC(=O)OCC1=CC=CC=C1", "benzyl acetate", "乙酸苄酯"),
    ("C(CCCCCCCCCCCCCCCCC)(=O)OCC1=CC=CC=C1", "benzyl octadecanoate", None),
    ("CC(=O)OC(C)(C)C", "tert-butyl acetate", "乙酸叔丁酯"),
    ("CC(=O)OC(C)C", "propan-2-yl acetate", "乙酸丙-2-基酯"),
    ("CC(=O)Oc1ccccc1", "phenyl acetate", "乙酸苯酯"),
    ("O=C(OCc1ccccc1)c1ccccc1", "benzyl benzoate", "苯甲酸苄酯"),
    # negative: linear alkyl esters must hold
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
    ("CC(=O)OCCC", "propyl acetate", "乙酸丙酯"),
]

# Must NOT collapse to bare phenyl/benzyl acetate (substituted/heteroaryl)
ester_alkoxy_sides__NOT_BARE = [
    "CC(=O)Oc1ccc(Cl)cc1",
    "CC(=O)Oc1ncccc1",
    "CC(=O)Oc1ccccc1C",
    "CC(=O)OCC1=CC=C(Cl)C=C1",
    "CC(=O)OCc1ccccn1",
]


@pytest.mark.parametrize("smiles", ester_alkoxy_sides__NOT_BARE)
def test_not_bare_phenyl_or_benzyl(smiles: str) -> None:
    r = SMILESNNamer().name(smiles)
    en = normalize_en(r.en or "")
    assert en not in {
        normalize_en("phenyl acetate"),
        normalize_en("benzyl acetate"),
    }


# ==========================================================================
# 合并自 test_complex_benzoate.py
# IUPAC: P-65.6 / P-65.1.1.1
# Layer: L2, L5
#
# Complex O-alkyl benzoate parent hit (Ph–C(=O)–O–R).
#
# When R is not simple linear/special alkoxy, still select benzoate parent instead
# of collapsing to ethane/benzene. Alcohol-side radical naming is deferred
# (alkoxy_complex); L5 emits bare benzoate / 苯甲酸酯 as infrastructure.
# ==========================================================================
complex_benzoate__COMPLEX = (
    "COC[C@]12CN(C)C3[C@@H]4[C@H](OC)[C@H]1[C@@]3([C@@H](OC)C[C@H]2O)"
    "[C@@H]1C[C@@]2(O)[C@H](OC(=O)c3ccccc3)[C@@H]1[C@]4(O)[C@@H](O)[C@@H]2OC"
)

complex_benzoate__CASES = [
    # complex O-alkyl: parent hit (not dual-complete vs gold)
    (complex_benzoate__COMPLEX, "benzoate", "苯甲酸酯"),
    # simple regressions
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "苯甲酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", complex_benzoate__CASES)
def test_complex_benzoate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_complex_not_ethane_collapse() -> None:
    r = SMILESNNamer().name(complex_benzoate__COMPLEX)
    assert r.success
    en = normalize_en(r.en)
    assert en != "ethane"
    assert "benzoate" in en


# ==========================================================================
# 合并自 test_principal_benzoate.py
# IUPAC: P-65.6.1 / P-22.1.3
# Layer: L2
#
# 验证 principal 规则管线对「苯环上直接连酯基」产出 benzoate 保留母体，
# 不依赖经典生产者通道（_try_arene_other_fg / _benzoate_parent）。
# ==========================================================================
principal_benzoate__CASES = [
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "苯甲酸乙酯"),
    ("CCCCCCOC(=O)c1ccccc1", "hexyl benzoate", "苯甲酸己酯"),
    ("CCCCCCCCOC(=O)c1ccccc1", "octyl benzoate", "苯甲酸辛酯"),
    ("CCCCCCCCCCCCCCCCOC(=O)c1ccccc1", "hexadecyl benzoate", "苯甲酸十六酯"),
]

# 负例：近邻但不应被 benzoate 规则误伤
principal_benzoate__NEGATIVE = [
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),  # 苯甲酸走 ACID 保留名
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),    # 开链酯走开链母体
]


def principal_benzoate___principal_name(smiles: str):
    """Only the rule-driven principal pipeline + L3-L5, no legacy producers."""
    mol = preprocess(smiles)
    if mol is None:
        return None
    info = analyze(mol)
    for parent in rule_driven_parent_candidates(info):
        packed = pack_parent_stem(parent, info["mol"])
        hit = try_candidate(info, packed, depth=0)
        if hit is not None and hit.success:
            return hit
    return None


def test_principal_benzoate_cases():
    for smiles, en, zh in principal_benzoate__CASES:
        r = principal_benzoate___principal_name(smiles)
        assert r is not None, f"principal 管线未能产出候选: {smiles}"
        assert r.success, f"装配失败 {smiles}: {r.meta}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"


def test_principal_benzoate_negative_untouched():
    for smiles, en, zh in principal_benzoate__NEGATIVE:
        r = principal_benzoate___principal_name(smiles)
        assert r is not None, f"principal 管线未产出: {smiles}"
        assert normalize_en(r.en) == normalize_en(en), f"EN {smiles}: got {r.en!r} want {en!r}"
        assert normalize_zh(r.zh) == normalize_zh(zh), f"ZH {smiles}: got {r.zh!r} want {zh!r}"
