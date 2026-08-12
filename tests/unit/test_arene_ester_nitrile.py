# IUPAC: P-65.6 / P-65.5.1 / P-66.5.1
# Layer: L2,L4,L5
"""Retained arene parents: alkyl benzoate, benzonitrile, benzoyl chloride.

Benzene + one principal FG on a ring carbon (ester / nitrile / acyl chloride).
Allow 0–2 extra ring halo / methyl; FG attach = locant 1.
Ester alkoxy: simple straight n-alkyl C1–C16 (no branching).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: alkyl benzoate (P-65.6), C1–C16 n-alkyl alcohol side
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "苯甲酸乙酯"),
    ("CCCCCCOC(=O)c1ccccc1", "hexyl benzoate", "苯甲酸己酯"),
    ("CCCCCCCCOC(=O)c1ccccc1", "octyl benzoate", "苯甲酸辛酯"),
    ("CCCCCCCCCCCCCCCCOC(=O)c1ccccc1", "hexadecyl benzoate", "苯甲酸十六酯"),
    ("CCOC(=O)c1ccc(Br)cc1Cl", "ethyl 4-bromo-2-chlorobenzoate", None),
    ("CCOC(=O)c1c(Br)cc(Cl)cc1", "ethyl 2-bromo-4-chlorobenzoate", None),
    # positive: benzonitrile (P-66.5.1)
    ("N#Cc1ccccc1", "benzonitrile", "苯甲腈"),
    # negative: must not steal existing retained / open-chain parents
    ("OC(=O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("O=Cc1ccccc1", "benzaldehyde", "苯甲醛"),
    ("CCCCCCC(=O)OC", "methyl heptanoate", "庚酸甲酯"),
    ("CC(=O)OCC", "ethyl acetate", "乙酸乙酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_arene_ester_nitrile_acyl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzoate_not_chain_ester() -> None:
    """Aromatic ester must not expand ring into chain alkanoate parent."""
    r = SMILESNNamer().name("COC(=O)c1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "methyl benzoate"
    assert "heptanoate" not in en
    assert "hexanoate" not in en
