# IUPAC: P-22.2.3
# Layer: L2
"""Saturated monocycle skeletal replacement stems (EN xa…azacycloalkane / ZH x杂环y烷).

Pure functions only — no RDKit, no namer. Sites are pre-numbered (locant, Z).
"""
from __future__ import annotations

import pytest

from namepredict.layer2.scaffold.sat_hetero_stem import (
    sat_hetero_stem_en,
    sat_hetero_stem_zh,
)

# (ring_size, sites, expected_en, expected_zh)
CASES = [
    # single hetero
    (6, [(1, 16)], "1-thiacyclohexane", "硫杂环己烷"),
    (7, [(1, 7)], "1-azacycloheptane", "氮杂环庚烷"),
    (5, [(1, 8)], "1-oxacyclopentane", "氧杂环戊烷"),
    # multi same
    (5, [(1, 8), (3, 8)], "1,3-dioxacyclopentane", "1,3-二氧杂环戊烷"),
    (8, [(1, 7), (5, 7)], "1,5-diazacyclooctane", "1,5-二氮杂环辛烷"),
    (6, [(1, 8), (4, 8)], "1,4-dioxacyclohexane", "1,4-二氧杂环己烷"),
    # multi mixed (ZH compact form locked)
    (6, [(1, 8), (4, 16)], "1-oxa-4-thiacyclohexane", "1,4-氧硫杂环己烷"),
    (6, [(1, 8), (3, 7)], "1-oxa-3-azacyclohexane", "1,3-氧氮杂环己烷"),
]


@pytest.mark.parametrize("size,sites,en,zh", CASES)
def test_sat_hetero_stem_pairs(size, sites, en, zh) -> None:
    assert sat_hetero_stem_en(size, sites) == en
    assert sat_hetero_stem_zh(size, sites) == zh


def test_unknown_z_returns_none() -> None:
    assert sat_hetero_stem_en(6, [(1, 9)]) is None
    assert sat_hetero_stem_zh(6, [(1, 9)]) is None
    assert sat_hetero_stem_en(6, [(1, 8), (4, 99)]) is None
    assert sat_hetero_stem_zh(6, [(1, 8), (4, 99)]) is None


def test_invalid_size_or_empty_returns_none() -> None:
    assert sat_hetero_stem_en(2, [(1, 8)]) is None
    assert sat_hetero_stem_zh(13, [(1, 8)]) is None
    assert sat_hetero_stem_en(6, []) is None
    assert sat_hetero_stem_zh(6, []) is None
