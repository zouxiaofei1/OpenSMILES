# IUPAC: P-65.3.2.5 (substitutive nomenclature — retained prefix sulfo)
# Layer: L3
"""Retained sulfo (-SO3H) substituent leaf in RetainedBackend via _try_registry_leaf."""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
# zh=None skips Chinese assertion.
CASES = [
    # sulfo as substituent on acetic acid (locant 2 is correct)
    ("O=S(=O)(O)CC(=O)O", "2-sulfoacetic acid", "2-磺基乙酸"),

    # (2R)-2-hydroxy-3-sulfopropanoic acid — the original failing molecule
    ("O=C(O)[C@@H](O)CS(=O)(=O)O", "(2R)-2-hydroxy-3-sulfopropanoic acid", "(2R)-2-羟基-3-磺基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_sulfo_retained(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
