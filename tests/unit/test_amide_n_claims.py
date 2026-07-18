"""Amide N sides: no false simple n_alkyl/n_phenyl; shared claims name complex N."""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.alkenamide import _amide_n_meta
from namepredict.layer2.parent_selector import select_parent
from namepredict.namer import SMILESNNamer

CYCLOHEPTYL = "O=C(NC1CCCCCC1)c1ccccc1"
DIMETHYLPHENYL = "O=C(Nc1cc(C)cc(C)c1)c1ccccc1"


def test_cycloheptyl_not_simple_n_alkyl():
    mol = Chem.MolFromSmiles(CYCLOHEPTYL)
    info = analyze(mol)
    meta = _amide_n_meta(info)
    assert "n_alkyl_n" not in meta
    assert "n_alkyl_ns" not in meta
    assert meta.get("n_block") is True


def test_dimethylphenyl_not_simple_n_phenyl():
    mol = Chem.MolFromSmiles(DIMETHYLPHENYL)
    info = analyze(mol)
    meta = _amide_n_meta(info)
    assert "n_phenyl" not in meta
    assert meta.get("n_block") is True


@pytest.mark.parametrize(
    "smiles,en,zh",
    [
        (CYCLOHEPTYL, "N-cycloheptylbenzamide", "N-环庚基苯甲酰胺"),
        (DIMETHYLPHENYL, "N-(3,5-dimethylphenyl)benzamide", "N-(3,5-二甲基苯基)苯甲酰胺"),
    ],
)
def test_complex_n_final_names(smiles, en, zh):
    result = SMILESNNamer().name(smiles)
    assert result.success
    assert normalize_en(result.en) == normalize_en(en)
    assert normalize_zh(result.zh) == normalize_zh(zh)


def test_simple_n_methyl_still_works():
    result = SMILESNNamer().name("c1ccccc1C(=O)NC")
    assert result.success
    assert normalize_en(result.en) == normalize_en("N-methylbenzamide")


def test_simple_n_phenyl_still_works():
    result = SMILESNNamer().name("c1ccccc1C(=O)Nc2ccccc2")
    assert result.success
    assert normalize_en(result.en) == normalize_en("N-phenylbenzamide")


def test_parent_owned_excludes_n_side_atoms():
    mol = Chem.MolFromSmiles(CYCLOHEPTYL)
    info = analyze(mol)
    parent = select_parent(info)
    owned = parent["owned_atoms"]
    # cycloheptyl carbons must be outside ownership
    n_root = parent.get("n_block_root")
    assert n_root is not None
    assert n_root not in owned
