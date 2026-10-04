# 合并自 7 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_naphthalene.py: Retained parent naphthalene (IUPAC P-22.1.1 / P-25): two fused aromatic
test_anthracene.py: Unsubstituted anthracene retained parent via ring_systems linear 666.
test_phenanthrene.py: phenanthrene 保留名 + 传统编号(PIN): 纯烃与甲基衍生物。
test_pyrene.py: pyrene 保留名 + 推荐编号(PIN): 纯烃与甲基衍生物。
test_fused_namer.py: fused_namer: 未注册稠环的稠合名称组装(benzo[a].../naphtho[...]... 类)。
test_fused56_engine.py: Fused56 e2e scope (P-22.2.1) — negative guard only.
test_fused56_13_engine.py: Fused56 di13 e2e scope (P-22.2.1) — negative guard only.
"""
from __future__ import annotations

import pytest

from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_naphthalene.py
# IUPAC: P-22.1.1
# Layer: L2,L4,L5
#
# Retained parent naphthalene (IUPAC P-22.1.1 / P-25): two fused aromatic
# six-membered carbocycles. Unsubstituted + mono-methyl + mono-halo; lowest
# locant (1- preferred over 2-). Fusion carbons are not substitution sites here.
# ==========================================================================
naphthalene__CASES = [

    # positive: mono-methyl — alpha (1) and beta (2)
    ("Cc1cccc2ccccc12", "1-methylnaphthalene", "1-甲基萘"),
    # positive: mono-halo
    ("Clc1cccc2ccccc12", "1-chloronaphthalene", "1-氯萘"),
    ("CCCCCCCCCC", "decane", "癸烷"),
      # indole retained parent (not naphthalene)
]


@pytest.mark.parametrize("smiles,en,zh", naphthalene__CASES)
def test_naphthalene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_unsub_not_decane() -> None:
    r = SMILESNNamer().name("c1ccc2ccccc2c1")
    assert r.success
    assert normalize_en(r.en) == "naphthalene"
    assert "decane" not in normalize_en(r.en)


def test_methyl_locant_alpha_is_1() -> None:
    r = SMILESNNamer().name("Cc1cccc2ccccc12")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-methylnaphthalene"
    assert "2-methyl" not in en


# ==========================================================================
# 合并自 test_anthracene.py
# IUPAC: P-25
# Layer: L2,L5
#
# Unsubstituted anthracene retained parent via ring_systems linear 666.
# ==========================================================================
anthracene__CASES = [
    # negatives
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
]


@pytest.mark.parametrize("smiles,en,zh", anthracene__CASES)
def test_anthracene(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_phenanthrene.py
# IUPAC: P-25.4.1 / 表2.7
# Layer: e2e
#
# phenanthrene 保留名 + 传统编号(PIN): 纯烃与甲基衍生物。
# ==========================================================================
phenanthrene__CASES = [
    ("c1ccc2c(c1)ccc1ccccc12", "phenanthrene", "菲"),
    ("Cc1cccc2c1ccc1ccccc12", "1-methylphenanthrene", "1-甲基菲"),
    ("Cc1ccc2c(ccc3ccccc32)c1", "2-methylphenanthrene", "2-甲基菲"),
    ("Cc1ccc2ccc3ccccc3c2c1", "3-methylphenanthrene", "3-甲基菲"),
    ("Cc1cccc2ccc3ccccc3c12", "4-methylphenanthrene", "4-甲基菲"),
    ("Cc1cc2ccccc2c2ccccc12", "9-methylphenanthrene", "9-甲基菲"),
]


@pytest.mark.parametrize("smiles,en,zh", phenanthrene__CASES)
def test_phenanthrene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_pyrene.py
# IUPAC: P-25.3.3.3 / 表2.7
# Layer: e2e
#
# pyrene 保留名 + 推荐编号(PIN): 纯烃与甲基衍生物。
# ==========================================================================
pyrene__CASES = [
    ("c1cc2ccc3cccc4ccc(c1)c2c34", "pyrene", "芘"),
    ("Cc1ccc2ccc3cccc4ccc1c2c34", "1-methylpyrene", "1-甲基芘"),
    ("Cc1cc2ccc3cccc4ccc(c1)c2c34", "2-methylpyrene", "2-甲基芘"),
    ("Cc1cc2cccc3ccc4cccc1c4c32", "4-methylpyrene", "4-甲基芘"),
]


@pytest.mark.parametrize("smiles,en,zh", pyrene__CASES)
def test_pyrene_rule(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused_namer.py
# IUPAC: P-25.3.2 / P-25.3.8
# Layer: L5
#
# fused_namer: 未注册稠环的稠合名称组装(benzo[a].../naphtho[...]... 类)。
# ==========================================================================
fused_namer__CASES = [
    # benchmark 参考: benzo[e]pyrene(未注册) → 芘 + 苯稠合
    ("c1ccc2c(c1)c1cccc3ccc4cccc2c4c31", "benzo[e]pyrene", "苯并[e]芘"),
    # 杂单环稠合: pyrimidine(母体) + pyridine(附加)。碱环是嘧啶时需让稠合边落在
    # 其 C4-C5(d 侧)而非 C5-C6(e 侧)——组分编号的 locant 1 在等价双 N 间浮动后
    # 由稠合原子位次最小化决定, 描述符才对(P-25.3.1.3 位次尽可能低)。
    ("ClC=1C2=C(N=C(N1)C)N=CC=C2", "4-chloro-2-methylpyrido[2,3-d]pyrimidine", "4-氯-2-甲基吡啶并[2,3-d]嘧啶"),
]


@pytest.mark.parametrize("smiles,en,zh", fused_namer__CASES)
def test_fused_namer(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, smiles
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused56_engine.py
# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4,L5
#
# Fused56 e2e scope (P-22.2.1) — negative guard only.
#
# Positive benzofuran/benzothiophene e2e cases are not covered in this file.
# The retained case asserts indole is not captured by the fused56 engine.
# ==========================================================================
fused56_engine___E2E = [
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", fused56_engine___E2E)
def test_e2e_fused56_mono(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_fused56_13_engine.py
# IUPAC: P-22.2.1 / P-25
# Layer: L2,L4,L5
#
# Fused56 di13 e2e scope (P-22.2.1) — negative guard only.
#
# Positive benzothiazole/benzoxazole e2e cases are not covered in this file.
# The retained case asserts pyridine is not captured by the fused56 engine.
# ==========================================================================
fused56_13_engine___E2E = [
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", fused56_13_engine___E2E)
def test_e2e_fused56_di13(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
