# IUPAC: P-65.7.1
# Layer: L1,L2,L5
"""Unsubstituted symmetrical open-chain alkanoic anhydrides.

Symmetric R–C(=O)–O–C(=O)–R; both acyl chains saturated, acyclic, no other
main FG. EN: C2 acetic anhydride; C≥3 …oic anhydride (acid stem + anhydride).
ZH: 酸名 + 酐 (乙酸酐, 丙酸酐, …).
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # negative: acid / ester / acyl chloride / ketone must not become anhydride
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CC(=O)C", "propan-2-one", "丙-2-酮"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkanoic_anhydride(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
