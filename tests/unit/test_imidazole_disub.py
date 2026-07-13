# IUPAC: P-14.3.4 / P-22.2.1
# Layer: L2
"""Five-membered heteroarene simple-sub limit 1→2 (IUPAC P-14.3.4 / P-22.2.1).

Retained parents furan/thiophene/pyrrole/imidazole/pyrazole may carry up to two
ring simple substituents: two monomethyl, two monohalo, or one of each. Each
side chain must remain pure monomethyl.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: iodo + methyl imidazole (gold dual target)
    ("IC=1NC=C(N1)C", "2-iodo-4-methyl-1H-imidazole", "2-碘-4-甲基-1H-咪唑"),
    # positive: dimethyl imidazole
    ("Cc1nc(C)c[nH]1", "2,4-dimethyl-1H-imidazole", "2,4-二甲基-1H-咪唑"),
    # positive: dichloro imidazole
    ("Clc1nc(Cl)c[nH]1", "2,4-dichloro-1H-imidazole", "2,4-二氯-1H-咪唑"),
    # positive: mono still ok (no regression)
    ("Cc1cn[nH]c1", "4-methyl-1H-pyrazole", "4-甲基-1H-吡唑"),
    # negative: unsubstituted keeps short Chinese parent
    ("c1cnc[nH]1", "1H-imidazole", "咪唑"),
    # negative: FG-bearing acid must not be simple imidazole parent
    ("Nc1[nH]cnc1C(=O)O", None, None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_imidazole_disub(smiles: str, en: str | None, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    if en is None:
        assert "imidazole" not in normalize_en(r.en) or "carboxylic" in normalize_en(r.en)
        return
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_acid_not_simple_imidazole() -> None:
    """Amino-imidazolecarboxylic acid is not a simple retained imidazole."""
    r = SMILESNNamer().name("Nc1[nH]cnc1C(=O)O")
    assert r.success
    en = normalize_en(r.en)
    assert en != "1h-imidazole"
    assert "2-amino" not in en or "carboxylic" in en or "acid" in en
