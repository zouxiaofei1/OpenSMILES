# IUPAC: P-67.1.3
# Layer: L2,L3,L5
"""磷酸整名（kind=phosphate 真集成，P 中心无机功能母体 + O-侧递归命名）。

磷酸作为整分子主官能团母体 kind=phosphate：链=[P 单原子]、owned={P+4O}，
O–R 臂由 L3 iter_claims 作 o_side 取代基整体递归命名（任意复杂：支链/芳环/
含 OH/胺/手性/糖环），L5 phosphate_names 按 n_oh 插 hydrogen/dihydrogen 组装。

此前实现把磷酸做在 namer 顶层短路、O-侧退化为直链烷基表；本测试锁定真 kind
路径下的任意复杂度 O-侧。臂内更高 FG（COOH 等 → phosphonooxy 降级）与
C–P 膦酸 / P–O–P / P(III) 不在本轮，谓词返回 None 放行旧路径。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
# 中文仅对简单词干（甲/乙/磷酸…）断言；复杂/立体臂 zh 留 None。
POSITIVE = [
    # 磷酸与全同中性酯
    ("OP(=O)(O)O", "phosphoric acid", "磷酸"),
    ("P(=O)(OC)(OC)OC", "trimethyl phosphate", "磷酸三甲酯"),
    ("O=P(OCC)(OCC)OCC", "triethyl phosphate", "磷酸三乙酯"),
    # 单/二烷基（O-侧为直链）磷酸酯
    ("O=P(O)(O)OC", "methyl dihydrogen phosphate", "磷酸二氢甲酯"),
    ("O=P(O)(O)OCC", "ethyl dihydrogen phosphate", "磷酸二氢乙酯"),
    ("O=P(O)(OC)OC", "dimethyl hydrogen phosphate", "磷酸氢二甲酯"),
    ("O=P(O)(OCC)OCC", "diethyl hydrogen phosphate", "磷酸氢二乙酯"),
    # 磷酸碱金属盐（k=0）
    ("[K+].O=P([O-])(O)O", "potassium dihydrogen phosphate", "磷酸二氢钾"),
    ("[Na+].[Na+].[O-]P(=O)([O-])O", "disodium hydrogen phosphate", "磷酸氢二钠"),
    # O-侧为任意复杂中性自由基（递归命名，直链限制解除）
    ("O=P(O)(O)OC(C)C", "propan-2-yl dihydrogen phosphate", None),
    ("O=P(O)(O)Oc1ccccc1", "phenyl dihydrogen phosphate", None),
    ("O=P(O)(O)OC[C@H](O)CO", "(2R)-2,3-dihydroxypropyl dihydrogen phosphate", None),
    ("CNCCOP(=O)(O)O", "2-(methylamino)ethyl dihydrogen phosphate", None),
    # O-侧为芳/饱和杂环（磷酸 O 连杂环臂）
    ("O=P(O)(O)OC1CCCCO1", "oxan-2-yl dihydrogen phosphate", None),
    ("O=P(O)(O)OC1CCNCC1", "piperidin-4-yl dihydrogen phosphate", None),
    ("O=P(O)(O)Oc1ccco1", "furan-2-yl dihydrogen phosphate", None),
    ("O=P(O)(O)Oc1ccncc1", "pyridin-4-yl dihydrogen phosphate", None),
    ("O=P(O)(O)Oc1ncccc1", "pyridin-2-yl dihydrogen phosphate", None),
    ("O=P(O)(O)Oc1nccs1", "1,3-thiazol-2-yl dihydrogen phosphate", None),
    ("O=P(O)(O)OCc1ccncc1", "pyridin-4-ylmethyl dihydrogen phosphate", None),
]


@pytest.mark.parametrize("smiles,en,zh", POSITIVE)
def test_phosphate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: {r.meta}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 全名回归守卫：不含磷 / 含羧酸的分子不得被 phosphate 母体截胡。
GUARDS = [
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("[Na+].[O-]C(=O)C", "sodium acetate", "乙酸钠"),
    ("CC(=O)O", "acetic acid", "乙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", GUARDS)
def test_non_phosphate_guard(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 结构近邻负例：detector 必须返回 None（不命中磷酸整名），维持旧行为。
# P(III) / P–O–P 焦磷酸 / P–C 膦酸 / 臂内更高 FG（COOH → phosphonooxy 属 M2）。
EXCLUDED = [
    "OP(O)O",                 # 亚磷酸 P(III)，无 P=O
    "O=P(O)(O)OP(=O)(O)O",   # 焦磷酸 P–O–P
    "CCP(=O)(O)O",            # 膦酸 P–C
    "O=C(O)COP(=O)(O)O",      # 臂内含 COOH（需磷酸降级，非整分子磷酸）
]


@pytest.mark.parametrize("smiles", EXCLUDED)
def test_phosphate_detector_excludes(smiles: str) -> None:
    from namepredict.layer0.preprocessor import preprocess
    from namepredict.layer0.salt import dissociate_salt
    from namepredict.layer1.phosphate import detect_phosphate_whole

    mol = preprocess(smiles)
    organic, _salt = dissociate_salt(mol)
    assert detect_phosphate_whole(organic) is None, smiles
