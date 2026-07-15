# IUPAC: P-22.2.3
# Layer: L2,L4,L5
"""Saturated monohetero replacement parents (table-out → x杂环y烷).

Retained `_MONO` hits (oxane / piperidine / oxolane / …) must keep existing
names. Table-out only uses skeletal replacement stems.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# positive: table-out sat monohetero
# negative: retained `_MONO` must not become 1-oxacyclohexane etc.
CASES = [
    # positive table-out
    ("C1CCSCC1", "1-thiacyclohexane", "硫杂环己烷"),
    ("C1CCCCNC1", "1-azacycloheptane", "氮杂环庚烷"),
    ("C1CSCCO1", "1-oxa-4-thiacyclohexane", "1,4-氧硫杂环己烷"),
    # negative retained
    ("C1CCOCC1", "oxane", "氧杂环己烷"),
    ("C1CCNCC1", "piperidine", "哌啶"),
    ("C1CCOC1", "oxolane", "氧杂环戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_repl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_table_out_not_methane() -> None:
    r = SMILESNNamer().name("C1CCSCC1")
    assert r.success
    assert normalize_en(r.en) != normalize_en("methane")
