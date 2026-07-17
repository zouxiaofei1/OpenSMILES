# IUPAC: P-66.1.1 / P-65.1.1.1
# Layer: L2,L3,L4,L5
"""Retained parent name benzamide (Ph–C(=O)–N; P-66.1.1 / 中文 6.5.8.1).

Single amide, single benzene, direct Ph–C(=O)–N. Ring ≤3 simple prefixes
(halo / n-alkyl·Me / OH / alkoxy / amino / nitro / CF3). N is H or simple
N-C1–C4 / N,N-dialkyl / N-phenyl (L3 N-extract). Prevent formamide collapse.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted retained parent
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("NC(=O)c1ccccc1", "benzamide", "苯甲酰胺"),
    # positive: ring simple prefixes (attach = 1)
    ("O=C(N)c1ccc(Cl)cc1", "4-chlorobenzamide", "4-氯苯甲酰胺"),
    ("O=C(N)c1ccc(O)cc1", "4-hydroxybenzamide", "4-羟基苯甲酰胺"),
    ("O=C(N)c1ccccc1C", "2-methylbenzamide", "2-甲基苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("O=C(N)c1ccc(N)cc1", "4-aminobenzamide", "4-氨基苯甲酰胺"),
    ("O=C(N)c1ccc([N+](=O)[O-])cc1", "4-nitrobenzamide", "4-硝基苯甲酰胺"),
    # positive: N-simple
    ("c1ccccc1C(=O)NC", "N-methylbenzamide", "N-甲基苯甲酰胺"),
    ("c1ccccc1C(=O)N(C)C", "N,N-dimethylbenzamide", "N,N-二甲基苯甲酰胺"),
    ("c1ccccc1C(=O)Nc2ccccc2", "N-phenylbenzamide", "N-苯基苯甲酰胺"),
    # negative: must not become benzamide / must keep existing names
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccccc1C(=O)O", "benzoic acid", "苯甲酸"),
    ("c1ccccc1C=O", "benzaldehyde", "苯甲醛"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide", "N-苯基乙酰胺"),
    ("O=C(N)C1CCCCC1", "cyclohexanecarboxamide", "环己烷甲酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_benzamide(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzamide_not_formamide_collapse() -> None:
    """Ar–CONH2 must not collapse to phenylformamide / formamide."""
    r = SMILESNNamer().name("c1ccccc1C(=O)N")
    assert r.success
    en = normalize_en(r.en)
    assert en == "benzamide"
    assert "formamide" not in en
    assert "phenyl" not in en
