# IUPAC: P-65.1.1
# Layer: L2,L3,L4,L5
"""Five-membered heteroarene carboxylic acids (IUPAC P-65.1.1 / P-22.2.1).

Retained / systematic monocyclic hetero5 parents with one ring-attached COOH:
furan / thiophene / 1H-pyrrole / 1H-imidazole / 1H-pyrazole.
Chinese suffix 甲酸 (aligned with pyridinecarboxylic / benzoic).
Allow ≤2 simple ring prefixes: halo / methyl / primary amino.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: unsubstituted mono-hetero5 carboxylic acids
    ("O=C(O)c1ccc[nH]1", "1H-pyrrole-2-carboxylic acid", "1H-吡咯-2-甲酸"),
    ("O=C(O)c1ccco1", "furan-2-carboxylic acid", "呋喃-2-甲酸"),
    ("O=C(O)c1cccs1", "thiophene-2-carboxylic acid", "噻吩-2-甲酸"),
    # positive: diazole carboxylic acids (NH=1, other N lowest)
    ("O=C(O)c1c[nH]cn1", "1H-imidazole-4-carboxylic acid", "1H-咪唑-4-甲酸"),
    ("O=C(O)c1cn[nH]c1", "1H-pyrazole-4-carboxylic acid", "1H-吡唑-4-甲酸"),
    # positive: mono methyl / halo / amino prefixes
    ("O=C(O)c1ccc(C)o1", "5-methylfuran-2-carboxylic acid", "5-甲基呋喃-2-甲酸"),
    ("Cc1ccc(C(=O)O)[nH]1", "5-methyl-1H-pyrrole-2-carboxylic acid", "5-甲基-1H-吡咯-2-甲酸"),
    ("O=C(O)c1occc1Cl", "3-chlorofuran-2-carboxylic acid", "3-氯呋喃-2-甲酸"),
    ("Nc1[nH]cnc1C(=O)O", "5-amino-1H-imidazole-4-carboxylic acid", "5-氨基-1H-咪唑-4-甲酸"),
    # negative: must not steal existing correct parents
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
    ("c1cc[nH]c1", "1H-pyrrole", "吡咯"),
    ("O=C(O)c1ccncc1", "pyridine-4-carboxylic acid", "吡啶-4-甲酸"),
    ("OC(=O)c1ccccn1", "pyridine-2-carboxylic acid", "吡啶-2-甲酸"),
    ("OC(=O)C1CCCCC1", "cyclohexanecarboxylic acid", "环己烷甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_hetero5_carboxylic(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_pyrrole_2_carboxylic_not_pentanoic() -> None:
    """Ring COOH on pyrrole must not collapse into open-chain pentanoic acid."""
    r = SMILESNNamer().name("O=C(O)c1ccc[nH]1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "1h-pyrrole-2-carboxylic acid"
    assert "pentanoic" not in en
