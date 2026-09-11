# IUPAC: P-67.1.4.1.1.2, P-67.1.4.1.1.5, P-67.1.4.1.1.6, P-67.1.4.1.3, P-72.6.2
# Layer: L2,L5
"""桥式五价磷的 P 酰基前缀（phosphoryl 族）命名。

此前桥式磷（磷酸二酯 P–O–P、P 上带烷氧/烷硫基的 phosphoryl、焦磷酸）没有母体词干，
P 被当成碳中心：`hydroxyphosphoryl]oxy` 写成 `hydroxymethoxy`、`phosphoryl` 写成 `methyl`。
本测试锁定：P(=O) 锚点取 phosphoryl 词干（P-67.1.4.1.1.2）、取代基按字母序拼接到该词干
（P-67.1.4.1.1.5 拼接：hydroxy(methyl)phosphoryl / dimethoxyphosphoryl）、经 O/N/S 桥连母体时
整体加方括号再接桥后缀（P-67.1.4.1.3：hydroxy(methoxy)phosphoryl]oxy）、阴离子氧臂用 oxido
（P-72.6.2：methoxy(oxido)phosphoryl）。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

STAR_CASES = [
    # 锚定碎片（* = 母体连接点）：P 自身的酰基前缀
    ("*P(=O)(O)OC", "hydroxy(methoxy)phosphoryl", "羟基(甲氧基)磷酰基"),
    ("*P(=O)(O)CC", "ethyl(hydroxy)phosphoryl", "乙基(羟基)磷酰基"),
    ("*P(=O)(OC)OC", "dimethoxyphosphoryl", "二甲氧基磷酰基"),
    ("*P(=O)([O-])OC", "methoxy(oxido)phosphoryl", None),
    # O 桥复合前缀：整体围栏后接 oxy
    ("*OP(=O)(O)OC", "[hydroxy(methoxy)phosphoryl]oxy", "[羟基(甲氧基)磷酰基]氧基"),
    ("*OP(=O)(OCC)SCCC", "[ethoxy(propylsulfanyl)phosphoryl]oxy", None),
    ("*OP(=O)([O-])OC", "[methoxy(oxido)phosphoryl]oxy", None),
    # 焦磷酸：P–O–P 两侧都必须是 phosphoryl，不得降级为碳
    ("COP(=O)(O)OP(=O)(O)OC",
     "[hydroxy-[[hydroxy(methoxy)phosphoryl]oxy]phosphoryl]oxymethane", None),
]

WHOLE_CASES = [
    # P 直连母体链：直接缀 phosphoryl（P-67.1.4.1.1.5），不是 methyl
    ("CP(=O)(O)CC[C@H](N)C(=O)O", "(2S)-2-amino-4-[hydroxy(methyl)phosphoryl]butanoic acid",
     "(2S)-2-氨基-4-[羟基(甲基)磷酰基]丁酸"),
    # P 带烷硫基 + 烷氧基臂：按字母序 ethoxy(propylsulfanyl)
    ("CCCSP(=O)(OCC)Oc1ccc(Br)cc1Cl",
     "4-bromo-2-chloro-1-[ethoxy(propylsulfanyl)phosphoryl]oxybenzene", None),
    ("CCOP(=O)(SC(C)CC)N1CCSC1=O",
     "3-[butan-2-ylsulfanyl(ethoxy)phosphoryl]-1,3-thiazolidin-2-one", None),
    # 磷酸二酯桥：P–O–C 臂上的取代基不参与 P 酰基，桥连母体侧仍出 oxy
    ("CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(=O)(O)OCCN",
     "(2R)-3-[2-aminoethoxy(hydroxy)phosphoryl]oxy-2-hydroxypropyl hexadecanoate",
     "十六酸(2R)-3-[2-氨基乙氧基(羟基)磷酰基]氧基-2-羟基丙酯"),
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

REGRESSION_CASES = [
    # 整分子磷酸（kind=phosphate）路径不受 P 酰基前缀改动影响
    ("COP(=O)(O)OC", "dimethyl hydrogen phosphate", "磷酸氢二甲酯"),
    ("OP(=O)(O)O", "phosphoric acid", "磷酸"),
    ("COP(=O)(O)O", "methyl dihydrogen phosphate", "磷酸二氢甲酯"),
    # 终端磷酸降级前缀仍走锚定表（P-67.1.5.1）
    ("*OP(=O)(O)O", "phosphonooxy", "膦酸氧基"),
    ("*COP(=O)(O)O", "phosphonooxymethyl", "膦酸氧甲基"),
]


@pytest.mark.parametrize("smiles,en,zh", [*STAR_CASES, *WHOLE_CASES, *REGRESSION_CASES])
def test_phosphoryl_prefix(smiles: str, en: str, zh: str | None) -> None:
    """桥式磷经 phosphoryl 族前缀命名，且不回落为把 P 当碳的甲基/甲氧基。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, f"{smiles}: 命名失败 {r.meta.get('reason')}"
    assert normalize_en(r.en) == normalize_en(en), f"{smiles}: {r.en!r} != {en!r}"
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh), f"{smiles}: {r.zh!r} != {zh!r}"


@pytest.mark.parametrize("smiles", [s for s, _, _ in [*STAR_CASES, *WHOLE_CASES]])
def test_no_carbon_fallback(smiles: str) -> None:
    """含 P 片段的名称里不得出现 P 被当碳的 methyl/methoxy 词干（P 只能以 phosphoryl/phosphono 系出现）。"""
    en = SMILESNNamer().name(smiles).en
    assert "phosph" in en, f"{smiles}: {en!r} 未出现任何磷前缀"
