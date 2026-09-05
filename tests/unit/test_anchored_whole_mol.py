"""整分子锚定键命中：顶层命名 *xxx 与取代基查表（resolve_name）一致。

之前顶层把 *xxx 当自由基母体硬算，暴露出与保留表不一致的编号/骨架错误
（quinolin-5-yl vs quinolin-2-yl、ethan-1-yl vs methoxy、FAIL 等）。
本测试要求：整分子 canonical == 锚定表键时，顶层命名直接返回保留名。
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.anchored_table import resolve_name


def _name(smiles: str, mode: str = "general"):
    return SMILESNNamer(name_mode=mode).name(smiles)


# ── registry 条目：顶层 *xxx == resolve_name（general/pin 双语）──

@pytest.mark.parametrize("smi,key", [
    ("*C(C)C", "isopropyl"),
    ("*C(C)(C)C", "tert-butyl"),
    ("*Cc1ccccc1", "benzyl"),
    ("*OC", "methoxy"),
    ("*S(C)(=O)=O", "methylsulfonyl"),
    ("*C#N", "cyano"),
    ("*O", "hydroxy"),
    ("*[O]", "oxidanyl"),
    ("*N", "amino"),
    ("*S", "sulfanyl"),
])
def test_top_level_star_matches_registry(smi, key):
    for mode in ("general", "pin"):
        r = _name(smi, mode)
        assert r.success, f"{smi} [{mode}] 命名失败: {r.meta.get('reason')}"
        en, zh = resolve_name(key, name_mode=mode)
        assert r.en == en, f"{smi} [{mode}] en: {r.en!r} != {en!r}"
        assert r.zh == zh, f"{smi} [{mode}] zh: {r.zh!r} != {zh!r}"


# ── 内联条目：顶层 *xxx == 内联 (en, zh) ──

@pytest.mark.parametrize("smi,en,zh", [
    ("*c1ccccc1", "phenyl", "苯基"),
    ("*C1CC1", "cyclopropyl", "环丙基"),
])
def test_top_level_star_matches_inline(smi, en, zh):
    r = _name(smi)
    assert r.success
    assert r.en == en
    assert r.zh == zh


# ── 命中路径标记 ──

def test_hit_sets_anchored_meta():
    r = _name("*C(C)C")
    assert (r.meta or {}).get("anchored") is True
    r2 = _name("*Cc1ccccc1", mode="pin")
    assert (r2.meta or {}).get("anchored") is True


# ── 普通分子不受影响（表键带 * 前缀，绝不误命中）──

def test_plain_molecules_not_hit():
    assert _name("c1ccccc1").en == "benzene"
    assert _name("CC(C)C").success
    assert _name("CC(C)C").en != "isopropyl"
    assert ( _name("CC(C)C").meta or {}).get("anchored") is not True


# ── 带附加取代基的自由基：整分子不在表，仍走管线 ──

def test_star_with_extra_substituent_not_hit():
    # 整分子 *c1ccc(Cl)cc1 不是单一锚定键 → 走管线，绝不误命中 phenyl
    r = _name("*c1ccc(Cl)cc1")
    assert (r.meta or {}).get("anchored") is not True
    assert r.en != "phenyl"


# ── *O/*[O]/*N 只服务整分子顶层，不进 L3 取代基查表 ──

def test_mononuclear_anchored_whole_mol_only():
    # 整分子 *O 顶层命中 anchored（hydroxy）
    r = _name("*O")
    assert r.en == "hydroxy"
    assert (r.meta or {}).get("anchored") is True


def test_hetero_anchored_not_used_as_substituent_lookup():
    # 多酮羰基氧不被 L3 误作羟基取代基（*O 不进取代基查表 + _ketone_fg_atoms 修复）
    assert _name("CC(=O)C(C)=O").en == "butane-2,3-dione"
    assert _name("CC(=O)CC(=O)C").en == "pentane-2,4-dione"
    # 真正醚氧取代基仍正常（*OC 在表，非 whole-only）
    assert _name("*OC").en == "methoxy"
