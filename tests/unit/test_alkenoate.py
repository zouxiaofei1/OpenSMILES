# IUPAC: P-65.6 / P-31.1
# Layer: L2,L4,L5
"""Open-chain monounsaturated monoesters (alkenoates).

Ester is the principal characteristic group (carbonyl C = locant 1); one
non-aromatic C=C on the acyl chain is expressed as -n-enoate / -n-烯酸…酯.
Alkoxy limited to unsubstituted C1–C4; no (E)/(Z) this cycle.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: acyclic mono-alkenoates, simple C1–C4 alkoxy (no E/Z)
    ("C=CC(=O)OC", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
    ("C=CC(=O)OCC", "ethyl prop-2-enoate", "丙-2-烯酸乙酯"),
    ("CC=CC(=O)OC", "methyl but-2-enoate", "丁-2-烯酸甲酯"),
    ("C=CCC(=O)OC", "methyl but-3-enoate", "丁-3-烯酸甲酯"),
    ("C=CC(=O)OCCC", "propyl prop-2-enoate", "丙-2-烯酸丙酯"),
    ("CCC=CC(=O)OCC", "ethyl pent-2-enoate", "戊-2-烯酸乙酯"),
    # negative: saturated esters, alkenoic acid, alkenenitrile must not break
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
    ("C=CC#N", "prop-2-enenitrile", "丙-2-烯腈"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenoate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
