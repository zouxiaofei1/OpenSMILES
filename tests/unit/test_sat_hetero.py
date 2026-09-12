# 合并自 3 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_sat_hetero.py: Saturated monohetero parent scope (P-22.2.2) — negative guard only.
test_sat_hetero_one.py: Sat-monohetero lactone / lactam scope (P-65.6.3.5.1) — negative guard only.
test_sat_hetero_alkyl.py: Sat-hetero ring-alkyl scope (P-22.2.2) — negative guard only.
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_sat_hetero.py
# IUPAC: P-22.2.2
# Layer: L2,L3,L4,L5
#
# Saturated monohetero parent scope (P-22.2.2) — negative guard only.
#
# Positive sat-hetero cases (oxolane, piperidine, …) are not covered in this
# file. The retained cases assert carbocycle / arene / open chain must not be
# named as Hantzsch–Widman parents.
# ==========================================================================
sat_hetero__CASES = [
    # negative: carbocycle / arene / open chain must not break
    ("C1CCCCC1", "cyclohexane", "环己烷"),
]


@pytest.mark.parametrize("smiles,en,zh", sat_hetero__CASES)
def test_sat_hetero(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_n_methylpiperidine() -> None:
    """N-alkyl on ring N：饱和杂环支持后命名 1-methylpiperidine（N 为 1 位）。"""
    r = SMILESNNamer().name("CN1CCCCC1")
    assert r.success
    assert normalize_en(r.en) == normalize_en("1-methylpiperidine")


def test_piperidine_carboxylic_not_oxolane() -> None:
    """Principal acid FG: must not invent sat-hetero parent name."""
    r = SMILESNNamer().name("O=C(O)C1CCCCN1")
    if not r.success:
        return  # fallback 已删：无候选显式失败
    en = normalize_en(r.en)
    assert "oxolane" not in en
    assert "oxane" not in en


# ==========================================================================
# 合并自 test_sat_hetero_one.py
# IUPAC: P-65.6.3.5.1 / P-66.1.5.1
# Layer: L2,L4,L5
#
# Sat-monohetero lactone / lactam scope (P-65.6.3.5.1) — negative guard only.
#
# Positive lactone/lactam cases are not covered in this file. The retained cases
# assert open-chain ester / amide keep their open-chain names.
# ==========================================================================
sat_hetero_one__CASES = [
    # negative: open-chain ester / amide keep open-chain names
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", sat_hetero_one__CASES)
def test_sat_hetero_one(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# ==========================================================================
# 合并自 test_sat_hetero_alkyl.py
# IUPAC: P-22.2.2 / P-14.3.4
# Layer: L2,L3,L4
#
# Sat-hetero ring-alkyl scope (P-22.2.2) — negative guard only.
#
# Positive ring-alkyl sat-hetero cases are not covered in this file. The retained
# case asserts acetamide is not named as a sat-hetero parent.
# ==========================================================================
sat_hetero_alkyl__CASES = [
    # negative: open-chain amide, not sat-hetero
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", sat_hetero_alkyl__CASES)
def test_sat_hetero_alkyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_n_methyl_not_c_methyl_piperidine() -> None:
    """N-甲基哌啶命名在 N 上（1-methylpiperidine），不误作 C-甲基。"""
    r = SMILESNNamer().name("CN1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1-methylpiperidine"
    assert en not in ("4-methylpiperidine", "2-methylpiperidine", "3-methylpiperidine")
