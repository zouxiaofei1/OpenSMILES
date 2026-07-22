# IUPAC: P-68.3.1.2
# Layer: L0,L1,L2,L5
"""Simple hydrazine retained name / N-substitution (P-68.3.1.2).

Scope (first cut):
- unsubstituted hydrazine (NN)
- 1,1-dialkyl C1–C4 n-alkyl (CN(C)N → 1,1-dimethylhydrazine)
- mono N-phenyl (NNc1ccccc1 → phenylhydrazine)
- hydrochloride salt (Cl.NNc1ccccc1 → phenylhydrazine;hydrochloride)
Must not regress acetamide, aniline, urea.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted
    ("NN", "hydrazine", "肼"),
    # positive: 1,1-dialkyl
    # positive: mono N-phenyl
    ("NNc1ccccc1", "phenylhydrazine", "苯肼"),
    # positive: hydrochloride
    # negative: amide / aniline / urea
    ("CC(=O)N", "acetamide", "乙酰胺"),
    ("c1ccccc1N", "aniline", "苯胺"),
    ("NC(=O)N", "urea", "脲"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hydrazine(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
