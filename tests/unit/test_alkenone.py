# IUPAC: P-64.2.1 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monoketones (alkenones / alkynones).

Ketone is the principal characteristic group; one non-aromatic C=C or C≡C is
expressed as -n-en-m-one / -n-烯-m-酮 or -n-yn-m-one / -n-炔-m-酮. Ketone locant
is minimized first (P-64.2.1); when tied, unsaturation takes the lower set
(P-31.1). E/Z via BondStereo when defined.
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenones
    ("C=CC(C)=O", "but-3-en-2-one", "丁-3-烯-2-酮"),
    ("CC=CC(C)=O", "pent-3-en-2-one", "戊-3-烯-2-酮"),
    ("C=CCC(C)=O", "pent-4-en-2-one", "戊-4-烯-2-酮"),
    ("CCC=CC(C)=O", "hex-3-en-2-one", "己-3-烯-2-酮"),
    # positive: alkynone
    ("C#CC(C)=O", "but-3-yn-2-one", "丁-3-炔-2-酮"),
    # E/Z mono-alkenone
    (r"C/C=C/C(C)=O", "(3E)-pent-3-en-2-one", "(3E)-戊-3-烯-2-酮"),
    # negatives: saturated ketone / ring / aromatic must not become alkenone
    ("CCC(C)=O", "butan-2-one", "丁-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenone(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
