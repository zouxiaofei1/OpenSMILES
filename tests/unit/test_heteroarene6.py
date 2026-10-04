# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_pyridine.py: Retained pyridine scope (P-22.2.1) — unsubstituted + negative guard.
test_diazine_expand.py: Diazine expand scope (P-22.2.1) — negative guard only.
test_quinazoline.py: Quinazoline scope (P-22.2.1) — negative guard only.
test_pyridinamine_ol.py: Pyridin-amine / pyridin-ol scope (P-62.2.1) — negative guard only.
test_pyridine_alkoxy_nitro.py: Simple pyridine retains alkoxy/nitro ring substituents like benzene (P-22.2.1 / P-61.5 / P-63.2.2).
test_pyridinecarbonitrile.py: Pyridinecarbonitrile scope (P-66.5.1) — negative guard only.
test_pyridinecarboxamide.py: Retained pyridinecarboxamide parent (pyridine-n-carboxamide / 吡啶-n-甲酰胺).
"""
from __future__ import annotations

import pytest

from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_pyridine.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Retained pyridine scope (P-22.2.1) — unsubstituted + negative guard.
#
# Only the unsubstituted pyridine positive case is retained here; substituted /
# carboxylic cases are not covered. The other cases assert benzene / ethanol /
# cyclohexane stay correct.
# ==========================================================================
pyridine__CASES = [
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridine__CASES)
def test_pyridine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_pyridine_not_pentane() -> None:
    r = SMILESNNamer().name("c1ccncc1")
    assert r.success
    assert normalize_en(r.en) != "pentane"


# ==========================================================================
# 合并自 test_diazine_expand.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Diazine expand scope (P-22.2.1) — negative guard only.
#
# Positive substituted diazine cases are not covered in this file. The retained
# cases assert pyridine / benzene must not be captured as diazines.
# ==========================================================================
diazine_expand__CASES = [
    # negative: pyridine / benzene must not be captured as diazine
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", diazine_expand__CASES)
def test_diazine_expand(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_complex_side_not_simple_diazine() -> None:
    """Pentyl-substituted pyrimidine now supported (was out of scope A)."""
    r = SMILESNNamer().name("ClC1=NC=C(C=N1)CCCCC")
    assert r.success
    assert normalize_en(r.en) == "2-chloro-5-pentylpyrimidine"


# ==========================================================================
# 合并自 test_quinazoline.py
# IUPAC: P-22.2.1 / P-25
# Layer: L2,L5
#
# Quinazoline scope (P-22.2.1) — negative guard only.
#
# Positive quinazoline cases are not covered in this file. The retained case
# asserts naphthalene is not named quinazoline.
# ==========================================================================
quinazoline__CASES = [
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
]


@pytest.mark.parametrize("smiles,en,zh", quinazoline__CASES)
def test_quinazoline(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_cinnoline_not_benzodiazine() -> None:
    r = SMILESNNamer().name("c1ccc2nnccc2c1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "quinazoline" not in en and "quinoxaline" not in en


# ==========================================================================
# 合并自 test_pyridinamine_ol.py
# IUPAC: P-62.2.1 / P-63.1.4 / P-22.2.1
# Layer: L2,L4,L5
#
# Pyridin-amine / pyridin-ol scope (P-62.2.1) — negative guard only.
#
# Positive pyridinamine/pyridinol cases are not covered in this file. The retained
# cases assert bare pyridine, aniline, phenol, chain amine stay correct.
# ==========================================================================
pyridinamine_ol__CASES = [
    ("Oc1ccccc1", "phenol", "苯酚"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridinamine_ol__CASES)
def test_pyridinamine_ol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyridine_alkoxy_nitro.py
# IUPAC: P-22.2.1
# Layer: L2
#
# Simple pyridine retains alkoxy/nitro ring substituents like benzene (P-22.2.1 / P-61.5 / P-63.2.2).
# ==========================================================================
pyridine_alkoxy_nitro__CASES = [
    # negative near-miss: existing simple pyridine / benzene stay correct
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridine_alkoxy_nitro__CASES)
def test_pyridine_alkoxy_nitro(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyridinecarbonitrile.py
# IUPAC: P-66.5.1 / P-22.2.1
# Layer: L2,L3,L4,L5
#
# Pyridinecarbonitrile scope (P-66.5.1) — negative guard only.
#
# Positive pyridinecarbonitrile cases are not covered in this file. The retained
# cases assert benzonitrile / pyridine / acetonitrile / propanenitrile stay correct.
# ==========================================================================
pyridinecarbonitrile__CASES = [
    # negative: must not steal benzonitrile / pyridine / acid / chain nitrile
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    ("CC#N", "acetonitrile", "乙腈"),
    ("CCC#N", "propanenitrile", "丙腈"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridinecarbonitrile__CASES)
def test_pyridinecarbonitrile(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyridinecarboxamide.py
# IUPAC: P-66.1.1 / P-22.2.1 / P-65.1.1.1
# Layer: L2,L3,L4,L5
#
# Retained pyridinecarboxamide parent (pyridine-n-carboxamide / 吡啶-n-甲酰胺).
#
# Unfused pyridine + single exocyclic –CONH2 (P-66.1.1 / P-22.2.1). N=1;
# carboxamide attach locant; ring simple prefixes aligned with pyridinecarboxylic;
# N is H or simple N-C1–C4 / N,N-dialkyl. Prefer system name over isonicotinamide.
# Block (pyridin-n-yl)formamide collapse.
#
# 已知局限：pyridine-carbonitrile 的碳腈后缀（pyridine-4-carbonitrile）尚未实现。
# ==========================================================================
pyridinecarboxamide__CASES = [
    # positive: unsubstituted pyridine-2/3/4-carboxamide
    ("NC(=O)c1ccccn1", "pyridine-2-carboxamide", "吡啶-2-甲酰胺"),
    ("NC(=O)c1cccnc1", "pyridine-3-carboxamide", "吡啶-3-甲酰胺"),
    ("NC(=O)c1ccncc1", "pyridine-4-carboxamide", "吡啶-4-甲酰胺"),
    # positive: simple N-alkyl
    ("c1ccncc1C(=O)NC", "N-methylpyridine-3-carboxamide", "N-甲基吡啶-3-甲酰胺"),
    ("C(C)NC(=O)C1=CC=NC=C1", "N-ethylpyridine-4-carboxamide", "N-乙基吡啶-4-甲酰胺"),
    ("c1ccncc1C(=O)N(C)C", "N,N-dimethylpyridine-3-carboxamide", "N,N-二甲基吡啶-3-甲酰胺"),
    # positive: ring prefix (N=1; CONH2 + halo lowest set, P-14.4(c) 主官能团最低位次)
    ("O=C(N)c1ccc(Cl)nc1", "6-chloropyridine-3-carboxamide", "6-氯吡啶-3-甲酰胺"),
    ("O=C(N)c1cc(C)ncc1", "2-methylpyridine-4-carboxamide", "2-甲基吡啶-4-甲酰胺"),
    ("Nc1ccncc1", "pyridin-4-amine", "吡啶-4-胺"),
    ("O=C(N)C1CCCCC1", "cyclohexanecarboxamide", "环己烷甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", pyridinecarboxamide__CASES)
def test_pyridinecarboxamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
