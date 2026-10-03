# 同位素命名：氘/氚词前缀、母体核素描述符、卤素核素前置。
"""tests/unit/test_isotope_naming.py — 同位素修饰化合物命名回归（P-8）。"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# IUPAC: P-8（P-81 符号 / P-82 同位素取代）
# Layer: L0（preprocessor 采集同位素到重原子属性）/ L5（assembler 渲染）
#
# ²H/³H 走 deuterio/tritio 词前缀：母体骨架上的氘带母体位次、倍数外置；
# 单碳母体（甲烷）与复合臂（甲基/甲氧基）走 P-16.5.1.3.1「首基平铺、
# 余基括注」。非氢核素走 (位次-质量数元素计数) 描述符，紧贴母体词干前，
# 单原子不写位次与计数。卤素核素走 (核素)卤素前缀 前置式，外层圆括号
# 不因核素描述符升级为方括号（P-82.3 排序、P-82.2.1 括号紧贴）。
#
# zh=None：金标中文「氘」与「氘代」两套用词并存且未判分，故只断言英文；
# 中英文均已与金标逐字一致者才写 zh。
# ==========================================================================
isotope__CASES = [
    # 氘：骨架多取代，倍数前置、位次全列
    ("[2H]C([2H])(C(=O)O)C([2H])([2H])C(=O)O", "2,2,3,3-tetradeuteriobutanedioic acid", None),
    ("[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]", "1,1,1,3,3,3-hexadeuteriopropan-2-one", None),
    ("[2H]c1c([2H])c([2H])c(CC(N)C(=O)O)c([2H])c1[2H]",
     "2-amino-3-(2,3,4,5,6-pentadeuteriophenyl)propanoic acid",
     "2-氨基-3-(2,3,4,5,6-五氘代苯基)丙酸"),
    # 氘：单碳母体 → 首基平铺、余基括注、无位次
    ("[2H]C(Cl)(Cl)Cl", "trichloro(deuterio)methane", None),
    # 氘：复合臂内倍数随词干入括号
    ("[2H]C([2H])(C1=C(N=CS1)CC)Br", "5-[bromo(dideuterio)methyl]-4-ethyl-1,3-thiazole", None),
    ("[2H]COC([2H])[C@H]1C[C@H](N(C1)C(=O)OC(C)(C)C)C=O",
     "tert-butyl (2S,4S)-4-[deuterio(deuteriomethoxy)methyl]-2-formylpyrrolidine-1-carboxylate", None),
    # 13C：同元素多原子 → 位次列表 + 计数下标，紧贴母体词干
    ("OC[C@H]1O[13C@H](O)[13C@H](O)[C@@H](O)[C@@H]1O",
     "(2S,3R,4S,5S,6R)-6-(hydroxymethyl)(2,3-13C2)oxane-2,3,4,5-tetrol",
     "(2S,3R,4S,5S,6R)-6-(羟甲基)(2,3-13C2)氧杂环己烷-2,3,4,5-四醇"),
    # 14C：单原子且位次不可定位（腈碳）→ 省位次与计数
    ("C1=C(C=C(C(=C1N)Cl)N)[14C]#N", "3,5-diamino-4-chloro(14C)benzonitrile",
     "3,5-二氨基-4-氯(14C)苯甲腈"),
    # 卤素核素：核素前置在普通卤素前缀名之前
    ("CNc1ccc(-c2nc3ccc(O)cc3s2)cc1[18F]", "2-[3-(18F)fluoro-4-(methylamino)phenyl]-1,3-benzothiazol-6-ol",
     "2-[3-(18F)氟-4-(甲氨基)苯基]-1,3-苯并噻唑-6-醇"),
    ("COC(=O)[C@H]1[C@@H](c2ccc([123I])cc2)C[C@@H]2CC[C@H]1N2CCCF",
     "methyl (1R,2S,3S,5S)-8-(3-fluoropropyl)-3-(4-(123I)iodophenyl)-8-azabicyclo[3.2.1]octane-2-carboxylate",
     "(1R,2S,3S,5S)-8-(3-氟丙基)-3-(4-(123I)碘苯基)-8-氮杂双环[3.2.1]辛烷-2-羧酸甲酯"),
]

# 无同位素的卤素/烷烃不得因同位素分支改动而变名（负向守卫）。
isotope_guard__CASES = [
    ("CNc1ccc(-c2nc3ccc(O)cc3s2)cc1F", "2-[3-fluoro-4-(methylamino)phenyl]-1,3-benzothiazol-6-ol", None),
    ("COC(=O)[C@H]1[C@@H](c2ccc(I)cc2)C[C@@H]2CC[C@H]1N2CCCF",
     "methyl (1R,2S,3S,5S)-8-(3-fluoropropyl)-3-(4-iodophenyl)-8-azabicyclo[3.2.1]octane-2-carboxylate", None),
    ("C(Cl)(Cl)(Cl)Cl", "tetrachloromethane", "四氯甲烷"),
]


@pytest.mark.parametrize("smiles,en,zh", isotope__CASES)
def test_isotope_naming(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("smiles,en,zh", isotope_guard__CASES)
def test_non_isotope_halo_unchanged(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
