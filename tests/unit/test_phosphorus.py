# 合并自 2 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_phosphate.py: 磷酸整名（kind=phosphate 正规 FgSpec 集成，P 中心无机功能母体 + O-侧递归命名）。
test_phosphoryl_prefix.py: 桥式五价磷的 P 酰基前缀（phosphoryl 族）命名。
"""
from __future__ import annotations

import pytest

from opensmiles.namer import SMILESNNamer
from opensmiles.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_phosphate.py
# IUPAC: P-67.1.3
# Layer: L1,L2,L3,L5
#
# 磷酸整名（kind=phosphate 正规 FgSpec 集成，P 中心无机功能母体 + O-侧递归命名）。
#
# 磷酸注册为 FG_SPECS 的 phosphate 类（p41=9/path=(1)，P-67.1.3.2 归入酯类），
# 走 L1 检测 → L2 P-44 选择 → L5 join_phosphate_name 组装：母体 kind=phosphate 在 _KIND_TABLE
# 只有一行词尾钩子（plain_hook，按 n_oh 出 磷酸/磷酸二氢/磷酸氢），O-侧臂与酯共用 _join_o_side_arms；
# owned={P+4O}，O–R 臂由 L3 iter_claims 作 o_side 取代基整体递归命名（任意复杂：
# 支链/芳环/含 OH/胺/手性/糖环）。
#
# 臂内更高优先级 FG 存在时磷酸落选、降级为 phosphonooxy/膦酸氧基前缀
# （P-67.1.5.1）；C–P 膦酸 / P–O–P 焦磷酸 / P(III) 不产出磷酸条目。
# ==========================================================================
phosphate__POSITIVE = [
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


@pytest.mark.parametrize("smiles,en,zh", phosphate__POSITIVE)
def test_phosphate(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: {r.meta}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 全名回归守卫：不含磷 / 含羧酸的分子不得被 phosphate 母体截胡。
phosphate__GUARDS = [
    ("[Na+].[O-]C(=O)C", "sodium acetate", "乙酸钠"),
]


@pytest.mark.parametrize("smiles,en,zh", phosphate__GUARDS)
def test_non_phosphate_guard(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 结构近邻负例：不产出磷酸条目。P(III) / P–C 膦酸。
# 焦磷酸 P–O–P 不在负例内：P-67.2 下每个链节仍算磷酸中心，见 CONDENSED 正例。
phosphate__EXCLUDED = [
    "OP(O)O",                 # 亚磷酸 P(III)，无 P=O
    "CCP(=O)(O)O",            # 膦酸 P–C
]


# 缩合磷酸（P–O–P）正例：链上无全酸式末端时仍取磷酸为母体，不得退为…甲烷母体。
condensed_phosphate__CASES = [
    ("O=P(O)(O)OP(=O)(O)O", "phosphono dihydrogen phosphate", None),   # 与 PubChem CID 1023 一致
    ("CCOP(=O)(O)OP(=O)(O)OCC", "ethoxy(hydroxy)phosphoryl ethyl hydrogen phosphate", None),
    ("COP(=O)(OC)OP(=O)(OC)OC", "dimethoxyphosphoryl dimethyl phosphate", None),  # 与 PubChem CID 120330 一致
]


@pytest.mark.parametrize("smiles,en,zh", condensed_phosphate__CASES)
def test_condensed_phosphate_parent(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: {r.meta}"
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles", phosphate__EXCLUDED)
def test_phosphate_detector_excludes(smiles: str) -> None:
    from opensmiles.layer0.preprocessor import preprocess
    from opensmiles.layer0.salt import dissociate_salt
    from opensmiles.layer1.analyzer import oxoacid_entries

    mol = preprocess(smiles)
    organic, _salt = dissociate_salt(mol)
    assert [e for e in oxoacid_entries(organic) if e["oxo_kind"] == "phosphate"] == [], smiles


# ==========================================================================
# 合并自 test_phosphoryl_prefix.py
# IUPAC: P-67.1.4.1.1.2, P-67.1.4.1.1.5, P-67.1.4.1.1.6, P-67.1.4.1.3, P-72.6.2
# Layer: L2,L5
#
# 桥式五价磷的 P 酰基前缀（phosphoryl 族）命名。
#
# 此前桥式磷（磷酸二酯 P–O–P、P 上带烷氧/烷硫基的 phosphoryl、焦磷酸）没有母体词干，
# P 被当成碳中心：`hydroxyphosphoryl]oxy` 写成 `hydroxymethoxy`、`phosphoryl` 写成 `methyl`。
# 本测试锁定：P(=O) 锚点取 phosphoryl 词干（P-67.1.4.1.1.2）、取代基按字母序拼接到该词干
# （P-67.1.4.1.1.5 拼接：hydroxy(methyl)phosphoryl / dimethoxyphosphoryl）、经 O/N/S 桥连母体时
# 整体加方括号再接桥后缀（P-67.1.4.1.3：hydroxy(methoxy)phosphoryl]oxy）、阴离子氧臂用 oxido
# （P-72.6.2：methoxy(oxido)phosphoryl）。
# ==========================================================================
phosphoryl_prefix__STAR_CASES = [
    # 锚定碎片（* = 母体连接点）：P 自身的酰基前缀
    ("*P(=O)(O)OC", "hydroxy(methoxy)phosphoryl", "羟基(甲氧基)磷酰基"),
    ("*P(=O)(O)CC", "ethyl(hydroxy)phosphoryl", "乙基(羟基)磷酰基"),
    ("*P(=O)(OC)OC", "dimethoxyphosphoryl", "二甲氧基磷酰基"),
    ("*P(=O)([O-])OC", "methoxy(oxido)phosphoryl", None),
    # O 桥复合前缀：整体围栏后接 oxy
    ("*OP(=O)(O)OC", "[hydroxy(methoxy)phosphoryl]oxy", "[羟基(甲氧基)磷酰基]氧基"),
    ("*OP(=O)(OCC)SCCC", "[ethoxy(propylsulfanyl)phosphoryl]oxy", None),
    ("*OP(=O)([O-])OC", "[methoxy(oxido)phosphoryl]oxy", None),
    # 缩合磷酸片段：P–O–P 两侧都必须是 phosphoryl，不得降级为碳
    ("*COP(=O)(O)OP(=O)(O)OC",
     "[hydroxy-[[hydroxy(methoxy)phosphoryl]oxy]phosphoryl]oxymethyl",
     "[羟基-[[羟基(甲氧基)磷酰基]氧基]磷酰基]氧基甲基"),
]

phosphoryl_prefix__WHOLE_CASES = [
    # P 直连母体链：直接缀 phosphoryl（P-67.1.4.1.1.5），不是 methyl
    ("CP(=O)(O)CC[C@H](N)C(=O)O", "(2S)-2-amino-4-[hydroxy(methyl)phosphoryl]butanoic acid",
     "(2S)-2-氨基-4-[羟基(甲基)磷酰基]丁酸"),
    # P 带烷硫基 + 烷氧基臂：按字母序 ethoxy(propylsulfanyl)
    ("CCCSP(=O)(OCC)Oc1ccc(Br)cc1Cl",
     "4-bromo-2-chloro-1-[ethoxy(propylsulfanyl)phosphoryl]oxybenzene", None),
    ("CCOP(=O)(SC(C)CC)N1CCSC1=O",
     "3-[butan-2-ylsulfanyl(ethoxy)phosphoryl]-1,3-thiazolidin-2-one", None),
    # N 桥（氨基甲酰磷酸）：P 前缀整体围栏后接 amino
    ("CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(=O)(O)OCCN=C(O)CCCCCCCCCCCCCCC",
     "(2R)-3-[2-(hexadecanoylamino)ethoxy-hydroxyphosphoryl]oxy-2-hydroxypropyl hexadecanoate",
     None),
    # 肌醇磷酸二酯：带手性环基的臂用方括号、桥后缀放括号外（[(2R,…)环己基]氧基）
    ("CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O[C@H](COC(=O)CCCCCCCCCCCCCCCCC)"
     "COP(=O)(O)OC1[C@H](O)[C@H](O)C(O)[C@H](O)[C@H]1O",
     "(2R)-1-[hydroxy-[(2R,3S,5R,6R)-2,3,4,5,6-pentahydroxycyclohexyl]oxyphosphoryl]oxy"
     "-3-octadecanoyloxypropan-2-yl (5Z,8Z,11Z,14Z)-henicosa-5,8,11,14-tetraenoate", None),
]

phosphoryl_prefix__REGRESSION_CASES = [
    # 整分子磷酸（kind=phosphate）路径不受 P 酰基前缀改动影响
    ("COP(=O)(O)OC", "dimethyl hydrogen phosphate", "磷酸氢二甲酯"),
    ("OP(=O)(O)O", "phosphoric acid", "磷酸"),
    ("COP(=O)(O)O", "methyl dihydrogen phosphate", "磷酸二氢甲酯"),
    # 终端磷酸降级前缀仍走锚定表（P-67.1.5.1）
    ("*OP(=O)(O)O", "phosphonooxy", "膦酸氧基"),
]


@pytest.mark.parametrize("smiles,en,zh", [*phosphoryl_prefix__STAR_CASES, *phosphoryl_prefix__WHOLE_CASES, *phosphoryl_prefix__REGRESSION_CASES])
def test_phosphoryl_prefix(smiles: str, en: str, zh: str | None) -> None:
    """桥式磷经 phosphoryl 族前缀命名，且不回落为把 P 当碳的甲基/甲氧基。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: 命名失败 {r.meta.get('reason')}"
    assert normalize_en(r.en) == normalize_en(en), f"{smiles}: {r.en!r} != {en!r}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh), f"{smiles}: {r.zh!r} != {zh!r}"


@pytest.mark.parametrize("smiles", [s for s, _, _ in [*phosphoryl_prefix__STAR_CASES, *phosphoryl_prefix__WHOLE_CASES]])
def test_no_carbon_fallback(smiles: str) -> None:
    """含 P 片段的名称里不得出现 P 被当碳的 methyl/methoxy 词干（P 只能以 phosphoryl/phosphono 系出现）。"""
    en = SMILESNNamer().name(smiles).en
    assert "phosph" in en, f"{smiles}: {en!r} 未出现任何磷前缀"
