# IUPAC: P-66.1.1.1.3
# Layer: L1,L2,L3,L5
"""Open-chain simple N-substituted amides (N-alkyl / N,N-dialkyl alkanamides).

Exactly one amide: carbonyl C with =O and N (not acid/ester/acyl chloride).
N bears 0–2 unsubstituted straight alkyls C1–C4 (sec or tert amide).
Acyl is unsubstituted straight formyl/acetyl/propanoyl…
No other principal FG; no ring/aromatic/C=C.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: N-alkyl / N,N-dialkyl alkanamides
    ("CC(=O)NC", "N-methylacetamide", "N-甲基乙酰胺"),
    ("CC(=O)NCC", "N-ethylacetamide", "N-乙基乙酰胺"),
    ("CC(=O)N(C)C", "N,N-dimethylacetamide", "N,N-二甲基乙酰胺"),
    ("O=CNC", "N-methylformamide", "N-甲基甲酰胺"),
    ("CCC(=O)NC", "N-methylpropanamide", "N-甲基丙酰胺"),
    ("CC(=O)N(C)CC", "N-ethyl-N-methylacetamide", "N-乙基-N-甲基乙酰胺"),
    # negative: primary amide / acid / aldehyde / sec amine must stay correct
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CC(=O)O", "acetic acid", "乙酸"),
    ("CC=O", "acetaldehyde", "乙醛"),
    ("CCNCC", "N-ethylethanamine", "N-乙基乙胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_n_alkyl_amide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
