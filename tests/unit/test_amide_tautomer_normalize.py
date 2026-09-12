# Layer: L0(tautomer), L1, L2, L5
"""CheBI 酰胺烯醇互变异构归一化：C(O)=N → C(=O)-NH。

非芳香中性 C(OH)=N 写法（亚胺酸/烯醇）应被当作酮式酰胺解析与命名；
带电 N / O⁻ 阴离子 / 芳香烯醇（吡啶酮类）不归一化。
"""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.tools.re import normalize_en
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.namer import SMILESNNamer


def _canon(smiles: str) -> str:
    """返回 SMILES 对应结构的规范串（含立体）。"""
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)


def test_preprocess_normalizes_enol_to_amide() -> None:
    """preprocess 后烯醇与酮式 canonical 相同，且 L1 检为 amide 而非醇。"""
    mol = preprocess("CC(O)=NC")
    assert _canon("CC(=O)NC") == Chem.MolToSmiles(mol)
    info = analyze(mol)
    assert info["has_amide"]
    assert not (info.get("hydroxyls") or [])


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(O)=NC", "N-methylacetamide"),
        ("CCCCN=CO", "N-butylformamide"),
        ("CC(O)=NCCc1cnc[nH]1", "N-[2-(1H-imidazol-5-yl)ethyl]acetamide"),
    ],
)
def test_enol_named_as_amide(smiles: str, expected: str) -> None:
    """烯醇写法命出的酰胺名与其酮式/参考一致。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(expected)


def test_enol_equals_keto_naming() -> None:
    """同一化合物的烯醇与酮式写法得到相同名称。"""
    n = SMILESNNamer()
    for enol, keto in [
        ("CC(O)=NC", "CC(=O)NC"),
        ("N[C@@H](CCC(O)=N)C(=O)O", "N[C@@H](CCC(=O)N)C(=O)O"),
    ]:
        e = n.name(enol)
        k = n.name(keto)
        assert e.success and k.success
        assert normalize_en(e.en) == normalize_en(k.en)


def test_stereo_preserved_across_normalization() -> None:
    """归一化只翻目标键级，不破坏其它手性中心。"""
    enol = preprocess("CC(O)=N[C@@H]1CCC[C@H]1O")
    keto = preprocess("CC(=O)N[C@@H]1CCC[C@H]1O")
    assert Chem.MolToSmiles(enol) == Chem.MolToSmiles(keto)
    r = SMILESNNamer().name("CC(O)=N[C@@H]1CCC[C@H]1O")
    assert r.success
    assert normalize_en(r.en) == normalize_en("N-[(1R,2R)-2-hydroxycyclopentyl]acetamide")


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)NC",  # 已是酮式，不改
        "Oc1ccccn1",  # 芳香烯醇（2-羟基吡啶→吡啶酮类），不改
        "CN(C)CCCSc1ccccc1[NH+]=C(O)/C=C/c1ccccc1",  # 带电 N 亚胺鎓，不改
    ],
)
def test_guard_no_false_normalization(smiles: str) -> None:
    """护栏：非中性/芳香/已酮式结构 preprocess 前后结构不变。"""
    mol = preprocess(smiles)
    assert Chem.MolToSmiles(mol) == _canon(smiles)
