# IUPAC: P-22.2.2 / P-14.3.4
# Layer: L2,L3,L4
"""Saturated monohetero parents: simple ring C-alkyl sides (P-22.2.2 / P-14.3.4).

Allow ring-C linear n-alkyl C1–C4 (multi) and retained branched (tert-butyl /
isopropyl / …) on piperidine / pyrrolidine / morpholine / oxolane / oxane, etc.
N-alkyl remains out of scope (_n_unsub). Monohalo on ring C stays allowed.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: mono-methyl regressions (hetero = 1 → fixed locants)
    ("CC1CCNCC1", "4-methylpiperidine", "4-甲基哌啶"),
    ("CC1CCCNC1", "3-methylpiperidine", "3-甲基哌啶"),
    ("CC1CCCCN1", "2-methylpiperidine", "2-甲基哌啶"),
    # positive: linear n-alkyl C2–C3 (was falling to alkane)
    ("CCC1CCNCC1", "4-ethylpiperidine", "4-乙基哌啶"),
    ("CCCC1CCNCC1", "4-propylpiperidine", "4-丙基哌啶"),
    ("CCCC1CCCNC1", "3-propylpiperidine", "3-丙基哌啶"),
    # positive: retained branched tert-butyl (lowest set → 3 on pyrrolidine)
    ("C(C)(C)(C)C1CNCC1", "3-tert-butylpyrrolidine", "3-叔丁基吡咯烷"),
    # positive: monohalo regression
    ("ClC1CCNCC1", "4-chloropiperidine", "4-氯哌啶"),
    # positive: multi-alkyl SMILES-invariant alpha (locant set tie → stem order)
    # ethyl < methyl alphabetically → lower locant on ethyl (P-14.5 / P-14.3.5)
    ("CC1CCC(CC)O1", "2-ethyl-5-methyloxolane", "2-乙基-5-甲基氧杂环戊烷"),
    ("CCC1CCC(C)O1", "2-ethyl-5-methyloxolane", "2-乙基-5-甲基氧杂环戊烷"),
    ("CC1CNCC(CC)C1", "3-ethyl-5-methylpiperidine", "3-乙基-5-甲基哌啶"),
    # negative: unsubstituted parent must stay bare name
    ("C1CCNCC1", "piperidine", "哌啶"),
    # negative: open-chain amide, not sat-hetero
    ("CC(=O)N", "acetamide", "乙酰胺"),
    # negative: plain mono-alkyl must keep lowest set (not alpha-flip)
    ("CC1CCCO1", "2-methyloxolane", "2-甲基氧杂环戊烷"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sat_hetero_alkyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_n_methyl_not_c_methyl_piperidine() -> None:
    """N-methylpiperidine is out of scope; must not invent a C-methyl name."""
    r = SMILESNNamer().name("CN1CCCCC1")
    assert r.success
    en = normalize_en(r.en)
    assert "piperidine" not in en
    assert en != "4-methylpiperidine"
    assert en != "2-methylpiperidine"
    assert en != "3-methylpiperidine"
