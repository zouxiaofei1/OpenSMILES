# IUPAC: P-22.2.1
# Layer: L2,L3,L4,L5
"""Expand retained diazine parents (IUPAC P-22.2.1 / P-14.3.4 / P-62.2.1).

Scope A — simple substituted pyrimidine / pyrazine / pyridazine:
  ring ≤3 (or 4 all-methyl) simple subs: halo, methyl, methoxy (+ nitro if aligned).
Scope B — pyrimidinamine parent (pyrimidine + ring primary NH2 suffix).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: unsubstituted diazines (must not regress)
    ("c1cncnc1", "pyrimidine", "嘧啶"),
    ("c1cnccn1", "pyrazine", "吡嗪"),
    ("c1ccnnc1", "pyridazine", "哒嗪"),
    # scope B: pyrimidinamine + simple ring sub
    ("Cc1cnc(N)nc1", "5-methylpyrimidin-2-amine", "2-氨基-5-甲基嘧啶"),
    ("COC=1C=NC(=NC1)N", "5-methoxypyrimidin-2-amine", "2-氨基-5-甲氧基嘧啶"),
    # bare pyrimidinamine (synthetic)
    ("Nc1ncncc1", "pyrimidin-4-amine", "嘧啶-4-胺"),
    ("Nc1ncnc(C)c1", "6-methylpyrimidin-4-amine", "4-氨基-6-甲基嘧啶"),
    # scope A: simple substituted diazine (no amine)
    ("BrC1=NC(=NC=C1)OC", "4-bromo-2-methoxypyrimidine", "4-溴-2-甲氧基嘧啶"),
    ("CC1=NC(C)=C(C)N=C1C", "2,3,5,6-tetramethylpyrazine", "2,3,5,6-四甲基吡嗪"),
    ("Clc1ncncc1", "4-chloropyrimidine", "4-氯嘧啶"),
    ("Cc1cncnc1", "5-methylpyrimidine", "5-甲基嘧啶"),
    ("COc1ncncc1", "4-methoxypyrimidine", "4-甲氧基嘧啶"),
    # negative: pyridine / benzene must not be captured as diazine
    ("c1ccncc1", "pyridine", "吡啶"),
    ("c1ccccc1", "benzene", "苯"),
    ("Nc1ccccn1", "pyridin-2-amine", "吡啶-2-胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_diazine_expand(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_pyrimidinamine_locant_not_flipped() -> None:
    """N fixed 1,3; 2-amine + 5-methyl must not become 4-amine + 5-methyl."""
    r = SMILESNNamer().name("Cc1cnc(N)nc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "5-methylpyrimidin-2-amine"
    assert "pyrimidin-4" not in en


def test_complex_side_not_simple_diazine() -> None:
    """Pentyl side chain is out of scope A (left for later)."""
    r = SMILESNNamer().name("ClC1=NC=C(C=N1)CCCCC")
    assert r.success
    en = normalize_en(r.en)
    assert "pentylpyrimidine" not in en
    assert "chloropyrimidine" not in en or "pentyl" not in en
