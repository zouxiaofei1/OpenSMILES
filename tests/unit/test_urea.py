# IUPAC: P-66.1.6.1.1
# Layer: L1,L2,L5
"""Simple mono-urea functional parent (carbonic diamide retained name urea).

Scope (first cut):
- unsubstituted urea
- mono N-aryl urea (phenyl / p-tolyl / halo-phenyl)
- N,N-dimethyl-N'-aryl urea with numerical locants 1/3 (benchmark style)
Must not regress true amides, carbamates, isocyanates, or aniline.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("NC(=O)N", "urea", "脲"),
    # positive: mono N-aryl
    ("c1ccc(NC(=O)N)cc1", "phenylurea", "苯基脲"),
    ("C1(=CC=C(C=C1)NC(=O)N)C", "p-tolylurea", "对甲苯基脲"),
    ("N=C(O)Nc1ccc(Cl)c(Cl)c1", "(3,4-dichlorophenyl)urea", None),
    # positive: 1,1-dimethyl-3-aryl (benchmark numerical locants)
    (
        "CN(C)C(=O)Nc1ccc(Cl)cc1",
        "3-(4-chlorophenyl)-1,1-dimethylurea",
        None,
    ),
    # negative: true amide / carbamate / isocyanate / aniline
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("CCOC(=O)N", "ethyl carbamate", "氨基甲酸乙酯"),
    ("CN=C=O", "methyl isocyanate", "异氰酸甲酯"),
    ("c1ccccc1N", "aniline", "苯胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_urea(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
