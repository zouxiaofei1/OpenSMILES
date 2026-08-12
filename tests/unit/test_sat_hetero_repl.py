# IUPAC: P-22.2.3
# Layer: L2,L4,L5
"""Sat-monohetero replacement parents (P-22.2.3) — negative guard only.

Positive table-out cases are not covered in this file. The retained case
asserts a thiane ring is not collapsed to methane.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

def test_table_out_not_methane() -> None:
    r = SMILESNNamer().name("C1CCSCC1")
    assert r.success
    assert normalize_en(r.en) != normalize_en("methane")
