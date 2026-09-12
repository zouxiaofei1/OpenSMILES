# IUPAC: P-67.1.3
# Layer: L1,L2,L3,L5
"""磷酸整名（kind=phosphate 正规 FgSpec 集成，P 中心无机功能母体 + O-侧递归命名）。

磷酸注册为 FG_SPECS 的 phosphate 类（p41=9/path=(1)，P-67.1.3.2 归入酯类），
走 L1 检测 → L2 P-44 选择 → L5 phosphate_names 组装：链=[P 单原子]、
owned={P+4O}，O–R 臂由 L3 iter_claims 作 o_side 取代基整体递归命名（任意复杂：
支链/芳环/含 OH/胺/手性/糖环），L5 按 n_oh 插 hydrogen/dihydrogen 组装。

臂内更高优先级 FG 存在时磷酸落选、降级为 phosphonooxy/膦酸氧基前缀
（P-67.1.5.1）；C–P 膦酸 / P–O–P 焦磷酸 / P(III) 不产出磷酸条目。
"""
from __future__ import annotations

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
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
    # 臂内/分子内更高优先级 FG（羧酸、羧酸酯）→ 磷酸落选并降级为 phosphonooxy 前缀（P-67.1.5.1）
    ("O=C(O)COP(=O)(O)O", "2-phosphonooxyacetic acid", "2-膦酸氧基乙酸"),
    ("CCCCCCCC[C@H](O)[C@H](CCCCCCCC(=O)O)OP(=O)(O)O",
     "(9S,10S)-10-hydroxy-9-phosphonooxyoctadecanoic acid", "(9S,10S)-10-羟基-9-膦酸氧基十八酸"),
    ("CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(=O)COP(=O)(O)O",
     "2-oxo-3-phosphonooxypropyl (9Z)-octadec-9-enoate", None),
    ("CCCCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCCCC)COP(=O)(O)O",
     "(2R)-1-decanoyloxy-3-phosphonooxypropan-2-yl octadecanoate", None),
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


# 结构近邻负例：不产出磷酸条目。P(III) / P–O–P 焦磷酸 / P–C 膦酸。
EXCLUDED = [
    "OP(O)O",                 # 亚磷酸 P(III)，无 P=O
    "O=P(O)(O)OP(=O)(O)O",   # 焦磷酸 P–O–P
    "CCP(=O)(O)O",            # 膦酸 P–C
]


@pytest.mark.parametrize("smiles", EXCLUDED)
def test_phosphate_detector_excludes(smiles: str) -> None:
    from namepredict.layer0.preprocessor import preprocess
    from namepredict.layer0.salt import dissociate_salt
    from namepredict.layer1.phosphate import phosphate_entries

    mol = preprocess(smiles)
    organic, _salt = dissociate_salt(mol)
    assert phosphate_entries(organic) == [], smiles
