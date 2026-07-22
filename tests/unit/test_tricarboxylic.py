# IUPAC: P-65.1.1 / P-72.2.2.1
# Layer: L2,L4,L5
"""Open-chain tricarboxylic acids and their fully deprotonated anions."""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.polycarboxylic import try_polycarboxylic_parent
from namepredict.layer4.numbering import _orient_polycarboxylic
from namepredict.namer import SMILESNNamer

CASES = [
    (
        "O=C(O)CC(C(=O)O)CC(=O)O",
        "propane-1,2,3-tricarboxylic acid",
        "丙烷-1,2,3-三羧酸",
    ),
    (
        "O=C([O-])CC(C(=O)[O-])CC(=O)[O-]",
        "propane-1,2,3-tricarboxylate",
        "丙烷-1,2,3-三羧酸根",
    ),
    (
        "O=C(O)C(O)(CC(=O)O)CC(=O)O",
        "2-hydroxypropane-1,2,3-tricarboxylic acid",
        "2-羟基丙烷-1,2,3-三羧酸",
    ),
    (
        "O=C(O)C(C(=O)O)C(C(=O)O)C(=O)O",
        "ethane-1,1,2,2-tetracarboxylic acid",
        "乙烷-1,1,2,2-四羧酸",
    ),
]


def test_tricarboxylate_does_not_fall_back_to_monoacid() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)[O-])CC(=O)[O-]")
    assert result.success
    assert "acetate" not in normalize_en(result.en)
    assert "hexanoate" not in normalize_en(result.en)


def test_partial_deprotonation_does_not_claim_polycarboxylate() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)O)CC(=O)O")
    assert "tricarboxylate" not in normalize_en(result.en)


@pytest.mark.parametrize("smiles,old_name", [
])
def test_unsupported_polyacids_close_ordinary_candidate_fallback(
    smiles: str, old_name: str,
) -> None:
    result = SMILESNNamer().name(smiles)
    assert not result.success
    assert result.meta["reason"] == "unsupported"
    assert normalize_en(result.en) != normalize_en(old_name)


def test_unsaturated_polyacid_parent_retains_all_core_unsaturation() -> None:
    info = analyze(Chem.MolFromSmiles(r"O=C([O-])/C=C\C(=C/C(=O)[O-])C(=O)[O-]"))
    parent = try_polycarboxylic_parent(info)
    assert parent is not None
    assert len(parent["double_bonds"]) == 2


def test_branched_carbon_skeleton_is_not_claimed() -> None:
    info = analyze(Chem.MolFromSmiles("CC(C(C(=O)O)C(C(=O)O)C(=O)O)C"))
    assert try_polycarboxylic_parent(info) is None


def test_polyacid_orientation_uses_lowest_carboxyl_locant_set() -> None:
    chain = [10, 11, 12, 13]
    parent = {"cooh_c_idxs": [11, 12, 13]}
    assert _orient_polycarboxylic(chain, parent, []) == list(reversed(chain))
