# IUPAC: P-65.6 / P-65.1.1.1
# Layer: L2, L5
"""Complex O-alkyl benzoate parent hit (Ph–C(=O)–O–R).

When R is not simple linear/special alkoxy, still select benzoate parent instead
of collapsing to ethane/benzene. Alcohol-side radical naming is deferred
(alkoxy_complex); L5 emits bare benzoate / 苯甲酸酯 as infrastructure.
"""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.namer import SMILESNNamer

# chebi-258: aconitine-like polycycle benzoate (error-driven pick)
COMPLEX = (
    "COC[C@]12CN(C)C3[C@@H]4[C@H](OC)[C@H]1[C@@]3([C@@H](OC)C[C@H]2O)"
    "[C@@H]1C[C@@]2(O)[C@H](OC(=O)c3ccccc3)[C@@H]1[C@]4(O)[C@@H](O)[C@@H]2OC"
)

CASES = [
    # complex O-alkyl: parent hit (not dual-complete vs gold)
    (COMPLEX, "benzoate", "苯甲酸酯"),
    # simple regressions
    ("COC(=O)c1ccccc1", "methyl benzoate", "苯甲酸甲酯"),
    ("CCOC(=O)c1ccccc1", "ethyl benzoate", "苯甲酸乙酯"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_complex_benzoate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_complex_benzoate_parent_kind() -> None:
    mol = Chem.MolFromSmiles(COMPLEX)
    info = analyze(mol)
    p = select_parent(info)
    assert p.get("kind") == "benzoate"
    assert p.get("alkoxy_complex") is True


def test_complex_not_ethane_collapse() -> None:
    r = SMILESNNamer().name(COMPLEX)
    assert r.success
    en = normalize_en(r.en)
    assert en != "ethane"
    assert "benzoate" in en
