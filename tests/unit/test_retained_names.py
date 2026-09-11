# IUPAC: P-62.2.1.1,P-66.1.1.4,P-29.2
# Layer: L3,L5
"""保留名（retained names）层：取代基保留名 anilino 与官能团衍生物保留前缀 carbamoyl/carbamoylamino/
carbamoyloxy/carbamothioylamino/sulfamoyl 必须真正进入输出，不得退化为系统名或逐原子拼装。

- P-62.2.1.1：phenylamino = anilino*（Glossary 818）；N-苯环带取代基时取代基前置于 anilino
  （2-methylanilino* = (2-methylphenyl)amino，Glossary 447；4-[(4-hydroxyanilino)methyl]phenol 为 PIN）。
- P-66.1.1.4.1：氨基甲酸（carbamic acid）的酰基保留前缀 carbamoyl；P-66.1.1.6 明确 ureido 不再
  使用，优选 carbamoylamino（P_1 附录）。
- P-66.1.1.4.2：磺酰胺保留前缀 sulfamoyl；(phenylamino)sulfonyl = phenylsulfamoyl*
  （Glossary 820），N-取代基与 sulfamoyl 融合（丁基(甲基)sulfamoyl）。

说明：gold/ChEBI 全量对 anilino 取 41:0、carbamoyl 取 54 处、sulfamoyl 取 20 处，故这些保留名
是 benchmark 与 IUPAC 一致的口径；而 vinyl/isobutyl/tosyl 等 general/not_rec 级保留名 gold 一律
取系统名（见 anchored_table.resolve_name），不在此列。
"""

from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# anilino：裸苯基、带环取代基、N-取代三种形态均取保留名 anilino/苯胺基
ANILINO = [
    ("O=C(O)c1ccccc1Nc1ccccc1", "2-anilinobenzoic acid", "2-苯胺基苯甲酸"),
    ("O=C(O)c1ccccc1Nc1ccc(Cl)cc1", "2-(4-chloroanilino)benzoic acid", "2-(4-氯苯胺基)苯甲酸"),
    ("Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl",
     "2-(2,6-dichloro-3-methylanilino)benzoic acid", "2-(2,6-二氯-3-甲基苯胺基)苯甲酸"),
    ("O=C(O)c1ccccc1N(C(=O)C)c1ccccc1",
     "2-(N-acetylanilino)benzoic acid", "2-(N-乙酰基苯胺基)苯甲酸"),
]


@pytest.mark.parametrize("smiles,en,zh", ANILINO)
def test_anilino_retained(smiles: str, en: str, zh: str) -> None:
    """anilino 取代基（裸/环取代/N-取代）输出保留名而非 phenylamino。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 锚定表保留前缀：片段级精确断言（碳/杂原子锚点）
RETAINED_FRAG = [
    ("*C(=O)N", "carbamoyl", "氨基甲酰基"),
    ("*NC(=O)N", "carbamoylamino", "氨基甲酰氨基"),
    ("*OC(=O)N", "carbamoyloxy", "氨基甲酰氧基"),
    ("*NC(N)=S", "carbamothioylamino", "氨基硫代羰基氨基"),
    ("*S(=O)(=O)N", "sulfamoyl", "氨磺酰基"),
]


@pytest.mark.parametrize("smiles,en,zh", RETAINED_FRAG)
def test_retained_fragment_prefixes(smiles: str, en: str, zh: str) -> None:
    """官能团衍生物保留前缀片段名精确匹配（不退化为 amino(oxo)methyl 等拼装式）。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


# 整分子：保留前缀嵌入组装名
RETAINED_WHOLE = [
    ("O=C(N)c1ccccc1C(=O)O", "2-carbamoylbenzoic acid", "2-氨基甲酰基苯甲酸"),
    ("N=C(O)N[C@@H](CS)C(=O)O",
     "(2R)-2-(carbamoylamino)-3-sulfanylpropanoic acid", "(2R)-2-(氨基甲酰氨基)-3-巯基丙酸"),
]


@pytest.mark.parametrize("smiles,en,zh", RETAINED_WHOLE)
def test_retained_prefix_whole_molecule(smiles: str, en: str, zh: str) -> None:
    """整分子组装名使用 carbamoyl/carbamoylamino 保留前缀。"""
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


def test_sulfamoyl_n_substituted_fuses() -> None:
    """N-取代磺酰胺取 sulfamoyl 且 N-取代基与之融合（P-66.1.1.4.2）。"""
    r = SMILESNNamer().name("COC1=CC=C(C(=O)NCCS(NCC=2C=NC=CC2)(=O)=O)C=C1")
    assert r.success
    assert normalize_en(r.en) == normalize_en(
        "4-methoxy-N-[2-(pyridin-3-ylmethylsulfamoyl)ethyl]benzamide")
    assert normalize_zh(r.zh) == normalize_zh(
        "4-甲氧基-N-[2-(吡啶-3-基甲氨基磺酰基)乙基]苯甲酰胺")


def test_general_level_retained_stays_systematic() -> None:
    """general/not_rec 级保留名（vinyl/isobutyl/tosyl）在 general 模式下仍取系统名，避免回退。"""
    cases = [
        ("C=Cc1ccccc1", "ethenylbenzene"),          # vinyl → gold 全量取 ethenyl（非 vinyl）
        ("CC(C)Cc1ccccc1", "2-methylpropyl"),       # isobutyl → 2-methylpropyl（非 isobutyl）
    ]
    for smiles, token in cases:
        r = SMILESNNamer().name(smiles)
        assert r.success
        assert token in normalize_en(r.en)
