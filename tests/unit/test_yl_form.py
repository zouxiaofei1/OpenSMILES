"""Tests for yl_form FG suffix→prefix conversion (P-63.2.2 / P-63.2.1 / P-62.2)."""
import pytest
from namepredict.layer3.yl_form import yl_form


# ── alcohol → alkoxy (P-63.2.2) ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanol", "甲醇", 1, "methoxy", "甲氧基"),
    ("ethanol", "乙醇", 2, "ethoxy", "乙氧基"),
    ("propan-1-ol", "丙-1-醇", 1, "propoxy", "丙氧基"),
    ("propan-2-ol", "丙-2-醇", 2, "propan-2-yloxy", "丙-2-基氧基"),
])
def test_alcohol_to_alkoxy(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── thiol → alkylsulfanyl (P-63.2.1) ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanethiol", "甲硫醇", 1, "methylsulfanyl", "甲硫基"),
    ("ethanethiol", "乙硫醇", 2, "ethylsulfanyl", "乙硫基"),
])
def test_thiol_to_sulfanyl(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── primary amine → alkylamino (P-62.2) ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("methanamine", "甲胺", 1, "methylamino", "甲氨基"),
    ("ethanamine", "乙胺", 2, "ethylamino", "乙氨基"),
])
def test_primary_amine_to_amino(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"


# ── non-converting: fallback to -n-yl ──
@pytest.mark.parametrize("en,zh,k,exp_en,exp_zh", [
    ("N-methylethanamine", "N-甲基乙胺", 1, "N-methylethanamin-1-yl", "N-甲基乙胺-1-基"),
    ("benzene", "苯", 1, "phenyl", "苯基"),
    ("ethane", "乙烷", 1, "ethan-1-yl", "乙烷-1-基"),
])
def test_fallback_to_n_yl(en, zh, k, exp_en, exp_zh):
    got_en, got_zh, paren = yl_form(en, zh, k)
    assert got_en == exp_en, f"EN: expected {exp_en}, got {got_en}"
    assert got_zh == exp_zh, f"ZH: expected {exp_zh}, got {got_zh}"
