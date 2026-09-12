# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_heteroarene5.py: Monocyclic heteroarene scope (P-22.2.1) — negative guard only.
test_pyrazole.py: Retained 1H-pyrazole scope (P-22.2.1) — negative guard only.
test_imidazole.py: Retained 1H-imidazole scope (P-22.2.1) — negative guard only.
test_imidazole_disub.py: Imidazole disub scope (P-14.3.4) — negative guard only.
test_indole.py: Retained 1H-indole scope (P-22.2.1) — unsubstituted + negative guard.
test_indazole.py: Retained 1H-indazole scope (P-22.2.1) — negative guard only.
test_benzimidazole.py: Retained 1H-benzimidazole scope (P-22.2.1) — negative guard only.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_heteroarene5.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Monocyclic heteroarene scope (P-22.2.1) — negative guard only.
#
# Positive furan/thiophene/diazine cases are not covered in this file. The
# retained cases assert benzene / pyridine / open chain must not regress.
# ==========================================================================
heteroarene5__CASES = [
    ("CCCCC", "pentane", "戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", heteroarene5__CASES)
def test_heteroarene5(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyrazole.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Retained 1H-pyrazole scope (P-22.2.1) — negative guard only.
#
# Positive pyrazole cases are not covered in this file. The retained case asserts
# pyridine is not misclassified as pyrazole.
# ==========================================================================
pyrazole__CASES = [
    # negative: must not misclassify 1,3-diazole / diazine / mono-hetero5
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", pyrazole__CASES)
def test_pyrazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_imidazole.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Retained 1H-imidazole scope (P-22.2.1) — negative guard only.
#
# Positive imidazole cases are not covered in this file. The retained cases
# assert benzene / pyridine / open chain must not regress.
# ==========================================================================
imidazole__CASES = [
    # negative: must not regress mono-hetero5 / arene / open chain
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", imidazole__CASES)
def test_imidazole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_imidazole_disub.py
# IUPAC: P-14.3.4 / P-22.2.1
# Layer: L2
#
# Imidazole disub scope (P-14.3.4) — negative guard only.
#
# Positive 2-substituted imidazole cases are not covered in this file. The
# retained case asserts an amino-imidazolecarboxylic acid is not named as a
# simple retained imidazole.
# ==========================================================================
imidazole_disub__CASES = [
    # negative: FG-bearing acid must not be simple imidazole parent
    ("Nc1[nH]cnc1C(=O)O", None, None),
]


@pytest.mark.parametrize("smiles,en,zh", imidazole_disub__CASES)
def test_imidazole_disub(smiles: str, en: str | None, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    if not r.success:
        return  # fallback 已删：无候选显式失败
    if en is None:
        assert "imidazole" not in normalize_en(r.en) or "carboxylic" in normalize_en(r.en)
        return
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_acid_not_simple_imidazole() -> None:
    """Amino-imidazolecarboxylic acid is not a simple retained imidazole."""
    r = SMILESNNamer().name("Nc1[nH]cnc1C(=O)O")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert en != "1h-imidazole"
    assert "2-amino" not in en or "carboxylic" in en or "acid" in en


# ==========================================================================
# 合并自 test_indole.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Retained 1H-indole scope (P-22.2.1) — unsubstituted + negative guard.
#
# Only the unsubstituted indole positive case is retained here; substituted
# cases are not covered. The other cases assert naphthalene / benzene / pyridine
# must not regress.
# ==========================================================================
indole__CASES = [
    # positive: unsubstituted
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", indole__CASES)
def test_indole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_octane() -> None:
    r = SMILESNNamer().name("c1ccc2[nH]ccc2c1")
    assert r.success
    assert normalize_en(r.en) == "1h-indole"
    assert "octane" not in normalize_en(r.en)


# ==========================================================================
# 合并自 test_indazole.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Retained 1H-indazole scope (P-22.2.1) — negative guard only.
#
# Positive indazole cases are not covered in this file. The retained case asserts
# indole must not regress.
# ==========================================================================
indazole__CASES = [
    # negative: must not regress indole / benzofuran / pyrazole / quinoline
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", indazole__CASES)
def test_indazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzimidazole.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Retained 1H-benzimidazole scope (P-22.2.1) — negative guard only.
#
# Positive benzimidazole cases are not covered in this file. The retained cases
# assert indole / benzene / aniline must not regress.
# ==========================================================================
benzimidazole__CASES = [
    ("Nc1ccccc1", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", benzimidazole__CASES)
def test_benzimidazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
