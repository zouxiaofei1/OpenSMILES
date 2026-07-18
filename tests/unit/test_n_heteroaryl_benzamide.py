# IUPAC: P-66.1.1 / P-29 / P-22.2.1
# Layer: L2,L3,L5
"""N-heteroaryl benzamide via recursive substituent path."""
from __future__ import annotations

import pytest
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

MVP = [
    (
        "c1ccccc1C(=O)Nc2ncccc2",
        "N-(pyridin-2-yl)benzamide",
        "N-(吡啶-2-基)苯甲酰胺",
    ),
    (
        "c1ccccc1C(=O)Nc2nccs2",
        "N-(1,3-thiazol-2-yl)benzamide",
        "N-(1,3-噻唑-2-基)苯甲酰胺",
    ),
    (
        "O=C(Nc1nccs1)c1ccc(OC)cc1",
        "4-methoxy-N-(1,3-thiazol-2-yl)benzamide",
        "4-甲氧基-N-(1,3-噻唑-2-基)苯甲酰胺",
    ),
    (
        "c1ccccc1C(=O)Nc2nc(C)cs2",
        "N-(4-methyl-1,3-thiazol-2-yl)benzamide",
        "N-(4-甲基-1,3-噻唑-2-基)苯甲酰胺",
    ),
]

REGRESS = [
    ("c1ccccc1C(=O)N", "benzamide", "苯甲酰胺"),
    ("O=C(N)c1ccc(OC)cc1", "4-methoxybenzamide", "4-甲氧基苯甲酰胺"),
    ("c1ccccc1C(=O)Nc2ccccc2", "N-phenylbenzamide", "N-苯基苯甲酰胺"),
    ("c1ccccc1C(=O)NC", "N-methylbenzamide", "N-甲基苯甲酰胺"),
    ("CC(=O)N", "acetamide", "乙酰胺"),
]


@pytest.mark.parametrize("smiles,en,zh", MVP)
def test_n_heteroaryl_benzamide_mvp(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", REGRESS)
def test_n_heteroaryl_benzamide_regress(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
