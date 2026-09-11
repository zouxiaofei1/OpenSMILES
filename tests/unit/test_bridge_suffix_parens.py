# IUPAC: P-14.4,P-16.5.1,P-63.2.2.1,P-65.6.3.2.3
# Layer: L3,L5
"""O/S/N 桥取代基前缀的括注边界与嵌套层级。

规则（gold/ChEBI 与 Blue Book 一致）：
- 前端 R 为自带围栏的复合取代基时，围栏只括前端，-oxy/-sulfanyl/-amino 留在括号外：
  `...oxan-2-yl]oxymethyl`、`...cyclohex-2-en-1-yl]amino]cyclohexyl`、`...chromen-7-yl]oxy`。
  依据 P-63.2.2.1.1 的 `(pyridin-2-yl)oxy`、P-65.6.3.2.3 的
  `3-[(pyridine-3-carbonyl)oxy]propanoic acid`。
- 前端为简单基（methyl/benzyl/直链 propan-2-yl）或 retained 单词酰基（acetyl/benzoyl/
  hexadecanoyl）时整段融合平铺：`benzylamino`、`propan-2-yloxy`、`hexadecanoyloxy`。
- 前端围栏已是方括号且桥为 -amino 时整体再括一层（P-16.5.2 嵌套标记）：
  `[[(1S)-1-carboxy-3-phenylpropyl]amino]propanoyl`。
- 中文侧与英文侧同形（括号闭在前端「基」后，氧基/硫基/氨基留括号外）。

每个期望值均取自 benchmarks/merged_benchmark.json 的同分子 gold（见行末 id）。
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer5.assembler_prefixes import (
    _front_needs_enclosure,
    _split_bridge_suffix,
)
from namepredict.namer import SMILESNNamer

# (smiles, expected_en, expected_zh_or_None, gold_id)
CASES = [
    # 糖苷 O 桥 + 甲基：前端括号闭在 -yl 后，oxy 与 methyl 融合在外
    (
        "C[C@@H]1O[C@H](OC[C@H]2O[C@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@@H](O)[C@H](O)[C@@H]1O",
        "(2S,3R,4S,5S,6R)-6-[[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]"
        "oxane-2,3,4,5-tetrol",
        None,
        "chebi-2801",
    ),
    # N 桥复合前端：前端方括号 + amino 整体再括一层（[[X]amino]cyclohexyl）
    (
        "OCC1=C[C@H](N[C@H]2C[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)"
        "[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
        "(2R,3R,4S,5S,6R)-2-[(1R,2R,3S,4S,6R)-2,3-dihydroxy-6-(hydroxymethyl)-4-"
        "[[(1S,4R,5S,6S)-4,5,6-trihydroxy-3-(hydroxymethyl)cyclohex-2-en-1-yl]amino]cyclohexyl]"
        "oxy-6-(hydroxymethyl)oxane-3,4,5-triol",
        None,
        "chebi-2376",
    ),
    # O 桥复合前端（无立体描述符）：括号闭在 -yl 后（P-63.2.2.1.1 的 (pyridin-2-yl)oxy）
    (
        "O1C(=CC=C1)CNC(COC=1C=CC=C2C=CC(=NC12)N1CCCC1)=O",
        "N-(furan-2-ylmethyl)-2-(2-pyrrolidin-1-ylquinolin-8-yl)oxyacetamide",
        "N-(呋喃-2-基甲基)-2-(2-吡咯烷-1-基喹啉-8-基)氧基乙酰胺",
        "tiers-102983",
    ),
    # 苄基型前端融合平铺（甲基直接接桥后缀）
    ("OCCNCc1ccccc1", "2-(benzylamino)ethanol", "2-(苄氨基)乙醇", None),
    # 直链 -yl 前端平铺
    ("BrC=1C=C(C=O)C=C(C1)OC(C)C", "3-bromo-5-propan-2-yloxybenzaldehyde", None, None),
    # 磺酰基桥 + 直链 -yl 前端：英文侧平铺（中文侧仍括注）
    (
        "CC(C)S(=O)(=O)c1ccc(C(N)=O)cc1",
        "4-propan-2-ylsulfonylbenzamide",
        None,
        None,
    ),
]


@pytest.mark.parametrize("smiles,en,zh,gold_id", CASES)
def test_bridge_suffix_parens(smiles: str, en: str, zh: str | None, gold_id: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


# 拆分判据单元测试：(stem, 期望拆分结果)
SPLIT_CASES = [
    ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yloxy",
     ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yl", "oxy")),
    ("(8R,9S,13S,14S)-13-methyl-17-oxo-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]"
     "phenanthren-3-yloxy",
     ("(8R,9S,13S,14S)-13-methyl-17-oxo-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]"
      "phenanthren-3-yl", "oxy")),
    ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yloxy", ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yl", "oxy")),
    ("(1S)-1-carboxy-3-phenylpropylamino", ("(1S)-1-carboxy-3-phenylpropyl", "amino")),
    # 不拆：简单保留基、直链 -yl、retained 单词酰基、苄基型前端
    ("benzyloxy", None),
    ("propan-2-yloxy", None),
    ("1,3-dihydroxypropan-2-yloxy", None),
    ("hexadecanoyloxy", None),
    ("(3,4-dichlorophenyl)methylamino", None),
    ("(1-methylimidazol-2-yl)sulfanylmethyl", None),
]


@pytest.mark.parametrize("stem,expected", SPLIT_CASES)
def test_split_bridge_suffix(stem: str, expected: tuple[str, str] | None) -> None:
    assert _split_bridge_suffix(stem) == expected


@pytest.mark.parametrize(
    "base,suf,expected",
    [
        ("(2S,3R,4S,5S,6R)-3,4,5-trihydroxy-6-methyloxan-2-yl", "oxy", True),
        ("5,6-dihydroxy-4-oxo-2-phenylchromen-7-yl", "oxy", True),
        ("(1S)-1-carboxy-3-phenylpropyl", "amino", True),
        ("propan-2-yl", "oxy", False),
        ("1,3-dihydroxypropan-2-yl", "oxy", False),
        ("hexadecanoyl", "oxy", False),
        ("(9Z)-octadec-9-enoyl", "oxy", True),  # 带位次的系统酰基名取 P-65.6.3.2.3 的 [(X)oxy] 式
        ("[4-methyl-2-(trifluoromethyl)phenyl]methyl", "sulfanyl", False),
    ],
)
def test_front_needs_enclosure(base: str, suf: str, expected: bool) -> None:
    assert _front_needs_enclosure(base, suf) is expected
