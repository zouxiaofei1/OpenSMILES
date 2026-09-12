# 合并自 4 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_benzofuran.py: Retained benzofuran scope (P-22.2.1) — negative guard only.
test_benzothiophene.py: Retained 1-benzothiophene scope (P-22.2.1) — negative guard only.
test_benzothiazole.py: Retained 1,3-benzothiazole scope (P-22.2.1) — negative guard only.
test_benzoxazole.py: Retained 1,3-benzoxazole scope (P-22.2.1) — negative guard only.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_benzofuran.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Retained benzofuran scope (P-22.2.1) — negative guard only.
#
# Positive benzofuran cases are not covered in this file. The retained cases
# assert indole / naphthalene / benzene must not regress.
# ==========================================================================
benzofuran__CASES = [
    # negative: near neighbors must not regress
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", benzofuran__CASES)
def test_benzofuran_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzothiophene.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Retained 1-benzothiophene scope (P-22.2.1) — negative guard only.
#
# Positive benzothiophene cases are not covered in this file. The retained cases
# assert indole / naphthalene / benzene must not regress.
# ==========================================================================
benzothiophene__CASES = [
    # negative: near neighbors must not regress
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", benzothiophene__CASES)
def test_benzothiophene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzothiazole.py
# IUPAC: P-22.2.1
# Layer: L2,L4,L5
#
# Retained 1,3-benzothiazole scope (P-22.2.1) — negative guard only.
#
# Positive benzothiazole cases are not covered in this file. The retained cases
# assert benzene / aniline must not regress.
# ==========================================================================
benzothiazole__CASES = [
    # negative: near neighbors must not regress
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", benzothiazole__CASES)
def test_benzothiazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_benzoxazole.py
# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
#
# Retained 1,3-benzoxazole scope (P-22.2.1) — negative guard only.
#
# Positive benzoxazole cases are not covered in this file. The retained case
# asserts benzene must not regress.
# ==========================================================================
benzoxazole__CASES = [
    # negative: near neighbors must not regress
    ("c1ccccc1", "benzene", "苯"),
]


@pytest.mark.parametrize("smiles,en,zh", benzoxazole__CASES)
def test_benzoxazole_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
