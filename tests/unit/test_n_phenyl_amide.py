# IUPAC: P-66.1.1.1.3
# Layer: L2,L3,L5
"""Simple open-chain N-phenyl alkanamides (P-66.1.1.1.3).

Exactly one amide with unsubstituted N-phenyl (N-phenylalkanamide).
Parent remains formamide/acetamide/propanamide…; N-phenyl is a prefix.
No N-benzyl, N-heteroaryl, or ring-substituted phenyl this round.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: N-phenyl alkanamides
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("CCC(=O)Nc1ccccc1", "N-phenylpropanamide", "N-苯基丙酰胺"),
    ("O=CNc1ccccc1", "N-phenylformamide", "N-苯基甲酰胺"),
    # negative: N-alkyl / primary stay correct; C-phenyl benzamide not N-phenyl
    ("CC(=O)NC", "N-methylacetamide", "N-甲基乙酰胺"),
    ("CC(=O)N(C)C", "N,N-dimethylacetamide", "N,N-二甲基乙酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccc(C(=O)N)cc1", "benzamide", "苯甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_phenyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
