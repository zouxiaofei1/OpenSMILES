"""Tests for retained_substituents registry and resolve_name dispatch."""
from __future__ import annotations

import pytest

from namepredict.layer3.retained_substituents import (
    IupacLevel,
    get_retained,
    resolve_name,
)

# ── registry completeness ──

EXPECTED_KEYS = {
    "tert-butyl", "isopropyl", "isobutyl", "sec-butyl",
    "neopentyl", "isopentyl", "2-methylbutan-2-yl", "3-methylbut-2-enyl",
    "vinyl", "allyl", "isopropenyl",
    "phenyl", "benzyl", "methoxy", "methylsulfanyl",
}


def test_all_expected_keys_present():
    for key in EXPECTED_KEYS:
        assert get_retained(key) is not None, f"missing key: {key}"


def test_no_duplicate_entries():
    seen_en = {}
    for key in EXPECTED_KEYS:
        e = get_retained(key)
        seen_en.setdefault(e.en, []).append(key)
    dupes = {en: keys for en, keys in seen_en.items() if len(keys) > 1}
    assert not dupes, f"duplicate en names: {dupes}"


# ── resolve_name dispatch ──

def test_general_always_uses_retained_name():
    assert resolve_name("vinyl", name_mode="general") == ("vinyl", "乙烯基")
    assert resolve_name("tert-butyl", name_mode="general") == ("tert-butyl", "叔丁基")


def test_pin_uses_systematic_for_general_level():
    assert resolve_name("vinyl", name_mode="pin") == ("ethenyl", "乙烯基")
    assert resolve_name("allyl", name_mode="pin") == ("prop-2-en-1-yl", "丙-2-烯-1-基")
    assert resolve_name("isopropyl", name_mode="pin") == ("propan-2-yl", "丙-2-基")


def test_pin_uses_systematic_for_not_recommended():
    assert resolve_name("isobutyl", name_mode="pin") == ("2-methylpropyl", "2-甲基丙基")
    assert resolve_name("sec-butyl", name_mode="pin") == ("butan-2-yl", "丁-2-基")


def test_pin_keeps_retained_for_pin_level():
    assert resolve_name("tert-butyl", name_mode="pin") == ("tert-butyl", "叔丁基")
    assert resolve_name("phenyl", name_mode="pin") == ("phenyl", "苯基")
    assert resolve_name("methylsulfanyl", name_mode="pin") == ("methylsulfanyl", "甲硫基")


def test_system_names_for_pin_level_also_pin():
    """PIN-level entries: retained == systematic (same PIN name)."""
    assert resolve_name("2-methylbutan-2-yl", name_mode="pin") == ("2-methylbutan-2-yl", "2-甲基丁-2-基")
    assert resolve_name("methoxy", name_mode="pin") == ("methoxy", "甲氧基")


def test_default_mode_is_general():
    assert resolve_name("vinyl") == ("vinyl", "乙烯基")
    assert resolve_name("isopropyl") == ("isopropyl", "异丙基")


# ── IupacLevel sanity ──

def test_iupac_level_values():
    assert get_retained("tert-butyl").level == IupacLevel.PIN
    assert get_retained("isopropyl").level == IupacLevel.GENERAL
    assert get_retained("isobutyl").level == IupacLevel.NOT_RECOMMENDED
