"""简单分子/单质查表：≤3 重原子的单元素与极简单分子取保留名（constants.SIMPLE_MOLECULES）。"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# 单元素分子：单质气体与臭氧（N₂/F₂/Cl₂/Br₂/I₂ 由 L2 杂原子链命名，见下）
element__CASES = [
    ("[H][H]", "hydrogen", "氢"),
    ("O=O", "oxygen", "氧"),
    ("[O-][O+]=O", "ozone", "臭氧"),
]

# 母体氢化物与卤化氢（P-21：保留名与二元名）
hydride__CASES = [
    ("O", "water", "水"),
    ("N", "ammonia", "氨"),
    ("S", "hydrogen sulfide", "硫化氢"),
    ("P", "phosphine", "膦"),
    ("[SiH4]", "silane", "硅烷"),
    ("[AsH3]", "arsine", "胂"),
    ("F", "hydrogen fluoride", "氟化氢"),
    ("Cl", "hydrogen chloride", "氯化氢"),
    ("Br", "hydrogen bromide", "溴化氢"),
    ("I", "hydrogen iodide", "碘化氢"),
]

# 简单氧化物 / 硫化物 / 氮化物 / 卤素化合物
oxide__CASES = [
    ("[C-]#[O+]", "carbon monoxide", "一氧化碳"),
    ("O=C=O", "carbon dioxide", "二氧化碳"),
    ("S=C=S", "carbon disulfide", "二硫化碳"),
    ("O=C=S", "carbonyl sulfide", "氧硫化碳"),
    ("O=S=O", "sulfur dioxide", "二氧化硫"),
    ("[N]=O", "nitric oxide", "一氧化氮"),
    ("[N-]=[N+]=O", "nitrous oxide", "一氧化二氮"),
    ("[O-][N+]=O", "nitrogen dioxide", "二氧化氮"),
    ("NO", "hydroxylamine", "羟胺"),
    ("OCl", "hypochlorous acid", "次氯酸"),
    ("NCl", "chloramine", "氯胺"),
]

# 单质金属（中性单原子）
metal__CASES = [
    ("[Fe]", "iron", "铁"),
    ("[Cu]", "copper", "铜"),
    ("[Zn]", "zinc", "锌"),
    ("[Ni]", "nickel", "镍"),
    ("[Ag]", "silver", "银"),
    ("[Au]", "gold", "金"),
    ("[Hg]", "mercury", "汞"),
    ("[Sn]", "tin", "锡"),
    ("[Pb]", "lead", "铅"),
    ("[Al]", "aluminium", "铝"),
    ("[W]", "tungsten", "钨"),
    ("[U]", "uranium", "铀"),
]

# 均一杂原子链（P-21.2.2）由 L2 命名，简单表不覆盖
heterane__CASES = [
    ("N#N", "diazyne"),
    ("NN", "diazane"),
    ("OO", "dioxidane"),
    ("FF", "difluorane"),
    ("ClCl", "dichlorane"),
    ("BrBr", "dibromane"),
    ("II", "diiodane"),
]

# 重原子 > 3 或非表物种不查表：仍由常规管线裁决
heavy_limit__CASES = [
    ("O=S(=O)=O", None),           # SO₃：4 重原子，表外，仍失败
    ("IC(I)I", "triiodomethane"),  # CHI₃：4 重原子，常规管线出三碘甲烷
    ("CO", "methanol"),            # 甲醇：2 重原子但非表物种，保留系统名
    ("CS", "methanethiol"),
]

# 多组分：简单物种作为片段参与拼装（水合物 / 单质混合）
component__CASES = [
    ("O.C1CO1", "oxirane; water", "环氧乙烷; 水"),
    ("[Fe].C1CO1", "iron; oxirane", "铁; 环氧乙烷"),
]


def _check(smiles: str, en: str, zh: str | None) -> None:
    """跑命名并比对 en/zh（样板同 test_multi_component）。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", element__CASES)
def test_element_molecule(smiles: str, en: str, zh: str | None) -> None:
    """单元素分子取单质名。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", hydride__CASES)
def test_hydride(smiles: str, en: str, zh: str | None) -> None:
    """母体氢化物/卤化氢取 P-21 保留名或二元名。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", oxide__CASES)
def test_simple_oxide(smiles: str, en: str, zh: str | None) -> None:
    """简单氧化物/硫化物/氮化物取惯用名。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en,zh", metal__CASES)
def test_element_metal(smiles: str, en: str, zh: str | None) -> None:
    """单质金属取元素名。"""
    _check(smiles, en, zh)


@pytest.mark.parametrize("smiles,en", heterane__CASES)
def test_heterane_chain_not_overridden(smiles: str, en: str) -> None:
    """均一杂原子链走 L2 系统名，简单表不覆盖。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,en", heavy_limit__CASES)
def test_heavy_atom_limit(smiles: str, en: str | None) -> None:
    """重原子超出 3 或非表物种不查表。"""
    r = SMILESNNamer().name(smiles)
    if en is None:
        assert not r.success
        return
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


@pytest.mark.parametrize("smiles,en,zh", component__CASES)
def test_simple_species_as_component(smiles: str, en: str, zh: str | None) -> None:
    """简单物种作为多组分片段时同样取表名。"""
    _check(smiles, en, zh)
