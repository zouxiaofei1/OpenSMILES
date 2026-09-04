# IUPAC: P-65.2.2 (环/稠环上的多 exocyclic -COOH → -Xcarboxylic acid)
# Layer: L4,L5
"""环烷/芳环上多个 exocyclic COOH 的系统命名。

Parent = carbocycle / benzene with multiplicity≥2 exocyclic COOH.
多羧酸用 -Xcarboxylic acid（非链式 Xanedioic 模板），必须带全部位次:
  cyclohexane-1,2-dicarboxylic acid / 环己烷-1,2-二羧酸
  benzene-1,4-dicarboxylic acid        / 苯-1,4-二羧酸
单羧酸（cyclohexanecarboxylic acid / benzoic acid）不得回归。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh")
CASES = [
    # 环己烷-1,2-二甲酸（用户示例，ortho）
    ("C1CCC(C(=O)O)C(C(=O)O)C1",
     "cyclohexane-1,2-dicarboxylic acid", "环己烷-1,2-二羧酸"),
    # 环己烷-1,4-二甲酸（para）
    ("O=C(O)C1CCC(C(=O)O)CC1",
     "cyclohexane-1,4-dicarboxylic acid", "环己烷-1,4-二羧酸"),
    # 苯-1,2-二甲酸（phthalic）
    ("C1=CC=C(C(=O)O)C(=C1)C(=O)O",
     "benzene-1,2-dicarboxylic acid", "苯-1,2-二羧酸"),
    # 苯-1,4-二甲酸（对苯二甲酸系统名）
    ("O=C(O)c1ccc(C(=O)O)cc1",
     "benzene-1,4-dicarboxylic acid", "苯-1,4-二羧酸"),
    # 苯-1,3,5-三甲酸（trimesic 系统名）
    ("O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1",
     "benzene-1,3,5-tricarboxylic acid", "苯-1,3,5-三羧酸"),
]

# 单羧酸对照：不得回归
NEG = [
    ("C1CCC(C(=O)O)CC1", "cyclohexanecarboxylic acid", "环己烷羧酸"),
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_ring_polycarboxylic(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh


@pytest.mark.parametrize("smiles,en,zh", NEG)
def test_monocarboxylic_no_regress(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, r.meta
    assert normalize_en(r.en) == normalize_en(en), r.en
    assert normalize_zh(r.zh) == normalize_zh(zh), r.zh
