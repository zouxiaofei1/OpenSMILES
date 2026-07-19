# IUPAC: P-65.3.1/P-65.3.2 (methylsulfinyl/methylsulfonyl retained prefixes)
# Layer: L3
"""Retained methylsulfinyl (CH3SO-), methylsulfonyl (CH3SO2-),
tosyl (4-Me-Ph-SO2-), triflyl (CF3-SO2-) retained prefix leaves."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # ── methylsulfonyl (PIN P-65.3.2) ──
    ("CS(=O)(=O)CC[C@H](N)C(=O)O",
     "(2S)-2-amino-4-methylsulfonylbutanoic acid", None),
    # methylsulfonyl on simple alkyl acid
    ("CS(=O)(=O)CCC(=O)O",
     "3-methylsulfonylpropanoic acid", None),
    # ── methylsulfinyl (P-65.3.1) ──
    ("CS(=O)CC[C@H](N)C(=O)O",
     "(2S)-2-amino-4-methylsulfinylbutanoic acid", None),
    # ── methylsulfanyl regression (must still work) ──
    ("CSCC[C@H](N)C(=O)O",
     "(2S)-2-amino-4-methylsulfanylbutanoic acid", None),
    ("CSCCCN",
     "3-methylsulfanylpropan-1-amine", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sulfonyl_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"naming failed for {smiles}"
    assert normalize_en(r.en) == normalize_en(en), \
        f"EN mismatch: {r.en} != {en}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
