# IUPAC: P-31.1 / P-44
# Layer: L2,L4,L5
"""Merge open-chain mono-FG alken* ParentKind into saturated kinds.

Unsaturated information lives only on parent fields double_bond /
double_bonds (+ L4 ene_locant / ene_locants). Names must stay correct;
parent.kind must be the saturated FG name.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_selector import select_parent
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_kind", "expected_en", "expected_zh")
UNSAT_CASES = [
    ("C=CC(=O)O", "acid", "prop-2-enoic acid", "丙-2-烯酸"),
    ("C=CC=O", "aldehyde", "prop-2-enal", "丙-2-烯醛"),
    ("C=CCO", "alcohol", "prop-2-en-1-ol", "丙-2-烯-1-醇"),
    ("COC(=O)C=C", "ester", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
    ("C=CC#N", "nitrile", "prop-2-enenitrile", "丙-2-烯腈"),
    ("C=CC(N)=O", "amide", "prop-2-enamide", "丙-2-烯酰胺"),
    ("O=C(O)C=CC(=O)O", "acid", "but-2-enedioic acid", "丁-2-烯二酸"),
]
SAT_CASES = [
    ("CC(=O)O", "acid", "acetic acid", "乙酸"),
    ("CCO", "alcohol", "ethanol", "乙醇"),
]
NAME_CASES = UNSAT_CASES + SAT_CASES


def _parent(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return select_parent(analyze(mol))


@pytest.mark.parametrize("smiles,kind,en,zh", UNSAT_CASES)
def test_merge_alken_parent_kind_unsat(
    smiles: str, kind: str, en: str, zh: str | None,
) -> None:
    p = _parent(smiles)
    assert p.get("kind") == kind
    assert p.get("double_bond") or p.get("double_bonds")


@pytest.mark.parametrize("smiles,kind,en,zh", SAT_CASES)
def test_merge_alken_parent_kind_sat(
    smiles: str, kind: str, en: str, zh: str | None,
) -> None:
    p = _parent(smiles)
    assert p.get("kind") == kind
    assert not p.get("double_bond") and not p.get("double_bonds")
