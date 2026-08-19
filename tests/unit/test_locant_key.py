# IUPAC: P-25.3.3.1 稠合碳字母位次
# Layer: L4
"""locant_key: 数字/字母位次统一排序。"""
from __future__ import annotations

from namepredict.layer4.locant_key import locant_key, locant_str_sort


def test_locant_key_parses_letter_suffix():
    assert locant_key("4a") == (4, "a")
    assert locant_key("10") == (10, "")
    assert locant_key(4) == (4, "")


def test_locant_str_sort_orders_mixed():
    locs = locant_str_sort(["5", "10", "4a", 4, "3"])
    assert locs == ["3", 4, "4a", "5", "10"]


def test_locant_key_numeric_ordering():
    # "10" > "2"(数值而非字典序)
    assert locant_key("10") > locant_key("2")
