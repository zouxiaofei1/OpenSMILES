# Layer: L0(charge)
"""L0 酸性质子重定位/电荷归一化：负电荷写在弱酸位(酚氧/烯醇 O⁻)而羧酸质子化时，
preprocess 把质子搬到弱酸位使负电荷收敛到羧酸根 → carboxylate 后缀命名。

IUPAC：carboxylate 阴离子后缀 P-72.2.2.1 / 酸根；酚/烯醇以中性 hydroxy 前缀表达。
"""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.namer import SMILESNNamer


def _canon(smiles: str) -> str:
    """返回 SMILES 对应结构的规范串（含立体）。"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


# (输入, 期望 en, 等价格式：负电已在羧酸/强酸根 + 中性弱酸)
RELO_CASES = [
    (
        "O=C(O)CCCCCCCCCCCCCCCCCCc1ccc([O-])cc1",
        "19-(4-hydroxyphenyl)nonadecanoate",
        "[O-]C(=O)CCCCCCCCCCCCCCCCCCc1ccc(O)cc1",
    ),
    (
        "N[C@H](Cc1ccc([O-])cc1)C(=O)O",
        "(2R)-2-amino-3-(4-hydroxyphenyl)propanoate",
        "N[C@H](Cc1ccc(O)cc1)C(=O)[O-]",
    ),
    (
        "N[C@@H](Cc1ccc([O-])cc1)C(=O)O",
        "(2S)-2-amino-3-(4-hydroxyphenyl)propanoate",
        "N[C@@H](Cc1ccc(O)cc1)C(=O)[O-]",
    ),
    (
        "O=C(O)c1[nH]c(-c2c[nH]c3ccccc23)cc1-c1c([O-])[nH]c2ccccc12",
        "3-(2-hydroxy-1H-indol-3-yl)-5-(1H-indol-3-yl)-1H-pyrrole-2-carboxylate",
        "O=C([O-])c1[nH]c(-c2c[nH]c3ccccc23)cc1-c1c(O)[nH]c2ccccc12",
    ),
]


def test_normalize_matches_equivalent_form() -> None:
    """弱酸负电写法 preprocess 后与其『负电在强酸 + 中性弱酸』等价格式 canonical 相同。"""
    for smiles, _en, equiv in RELO_CASES:
        mol = preprocess(smiles)
        assert Chem.MolToSmiles(mol) == _canon(equiv)


def test_stereo_preserved_across_normalization() -> None:
    """电荷重定位只改远处 O 的 FormalCharge/H，不翻转 Cα 手性中心。"""
    r_pp = Chem.MolToSmiles(preprocess("N[C@H](Cc1ccc([O-])cc1)C(=O)O"))
    s_pp = Chem.MolToSmiles(preprocess("N[C@@H](Cc1ccc([O-])cc1)C(=O)O"))
    assert "[C@H]" in r_pp and "[C@@H]" in s_pp


@pytest.mark.parametrize("smiles,en,equiv", RELO_CASES)
def test_weak_anion_named_as_carboxylate(smiles: str, en: str, equiv: str) -> None:
    """负电在弱酸位(酚氧/烯醇)而羧酸质子化的输入，命出的名 = 负电在羧酸的等价格式名。"""
    n = SMILESNNamer()
    r = n.name(smiles)
    re = n.name(equiv)
    assert r.success and re.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_en(r.en) == normalize_en(re.en)


def test_chebi433_bilingual() -> None:
    """chebi-433：羧酸根中文词尾 '酸根'，酚以 '羟基' 前缀表达。"""
    r = SMILESNNamer().name("O=C(O)CCCCCCCCCCCCCCCCCCc1ccc([O-])cc1")
    assert r.success
    assert normalize_zh(r.zh) == normalize_zh("19-(4-羟基苯基)十九酸根")


@pytest.mark.parametrize(
    "smiles",
    [
        "[O-]c1ccccc1",  # 纯酚盐：无强酸供体，不改
        "CCCCCCCCCCCC(=O)[O-]",  # 纯羧酸根：无弱受体，不改
        "[O-]C(=O)CCCO",  # 羧酸根 + 中性醇：不反搬醇(弱供体不作供体)
        "CC(=O)O",  # 中性羧酸：无弱受体，不改
        "O=[N+]([O-])c1ccc(C(=O)O)cc1",  # 硝基 O⁻ 邻 +1 内平衡，被跳过
        "COc1ccc([O-])cc1",  # 弱受体无强供体，不改
        "CC(=O)[O-].CCO",  # 分离片段：羧酸根片段无弱受体，不改
        "CC(=O)OCCO",  # 酯+醇：无弱受体/强供体配对，不改
        "Oc1ccccc1",  # 中性酚：无阴离子位，不改
    ],
)
def test_guard_no_false_relocation(smiles: str) -> None:
    """护栏：无 (弱受体, 强酸供体) 配对的输入 preprocess 前后结构不变。"""
    mol = preprocess(smiles)
    assert Chem.MolToSmiles(mol) == _canon(smiles)
