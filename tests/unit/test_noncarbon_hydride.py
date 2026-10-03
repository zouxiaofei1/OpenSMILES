"""非碳母体氢化物（P-21 / P-41 类 21–39）：杂原子自任母体、余者作取代基，非标准键数带 λ。

Layer: L1（heterane FG）/ L2（单原子骨架）/ L5（取名）
"""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

n = SMILESNNamer()


CASES = [
    # 四配位硫：四条臂各成一取代基，键数 4 非标准 → λ4
    ("CCS(C)(C(=C)C)SC",
     "ethyl-methyl-methylsulfanyl-prop-1-en-2-yl-lambda4-sulfane"),
    # 三配位碘：非标准 → λ3；片段态取 -ylidene
    ("C1=CSC=C1/C=C(/C=O)\\I=N",
     "(2Z)-2-(imino-lambda3-iodanyl)-3-thiophen-3-ylprop-2-enal"),
    # λ4 硫作取代基（片段态转 -yl）
    ("CS(C)(C)C1=CC=C(C(=C)C1)N",
     "6-methylidene-4-(trimethyl-lambda4-sulfanyl)cyclohexa-1,3-dien-1-amine"),
    # λ6 硫（全氟）作取代基
    ("CC(CSC1=CC=CC=C1C(F)(F)F)(C(=O)NC2=CC(=CC=C2)S(F)(F)(F)(F)F)O",
     "2-hydroxy-2-methyl-N-[3-(pentafluoro-lambda6-sulfanyl)phenyl]-3-[2-(trifluoromethyl)phenyl]sulfanylpropanamide"),
]


@pytest.mark.parametrize("smiles,gold", CASES)
def test_heterane_lambda(smiles, gold):
    """杂原子烃母体带 λ 命名（P-14.1.3 / P-21.1.2.1）。"""
    assert n.name(smiles).en == gold


# ── 边界：标准价与其它 FG 通路不受影响 ──────────
BOUNDARY = [
    ("CSC", "methylsulfanylmethane"),        # 硫醚 S(II) 标准价：P-41 类 41 排在碳之后，仍以碳为母体
    ("CCO", "ethanol"),                      # 醇不受影响
    ("CC(=O)O", "acetic acid"),              # 羧酸不受影响
    ("C[S+](C)(C)C", "tetramethylsulfanium"),  # 阳离子走 -ium 通路，不标 λ
]


@pytest.mark.parametrize("smiles,gold", BOUNDARY)
def test_heterane_does_not_claim_standard_valence(smiles, gold):
    """标准价杂原子与既有 FG 通路不被杂原子烃抢走。"""
    assert n.name(smiles).en == gold


def test_oxo_center_stays_with_oxoacid_path():
    """带 =O 的磷/硫仍由含氧酸/磷酰通路命名，不出 λ（P-41 类 7/9 高于类 21–39）。"""
    assert "lambda" not in n.name("O=P(O)(OC)OC").en
    assert "lambda" not in n.name("CS(C)(=O)=O").en


def test_thioxo_phosphorus_stays_phosphinothioyl():
    """带 =S 的磷仍走 phosphinothioyl 保留前缀（P-67.1.4.1.1.4），不标 λ。"""
    en = n.name("COP(=S)(OC)OC").en
    assert "phosphinothioyl" in en and "lambda" not in en


def test_normalize_folds_lambda_superscript():
    """λ⁵ 与 λ5 判分等价（英文的 lambda5 与 λ5 各自保留）。"""
    assert normalize_zh("三甲基-λ⁴-硫代") == normalize_zh("三甲基-λ4-硫代")
    assert normalize_en("trimethyl-lambda4-sulfanyl") == "trimethyl-lambda4-sulfanyl"


# ── 均一杂原子链（P-21.2.2 + P-31.1 不饱和）──────────
CHAIN_CASES = [
    ("NN", "diazane"),
    ("NNN", "triazane"),
    ("NNNNNNNNN", "nonaazane"),          # 倍数词末元音不省：nonaazane 非 nonazane（P-16.7.2）
    ("SS", "disulfane"),
    ("SSSS", "tetrasulfane"),            # 末端 -SH 的潜在官能性被忽略（P-21.2.2）
    ("[SiH3][SiH2][SiH2][SiH2][SiH3]", "pentasilane"),
    ("PP", "diphosphane"),
    ("PPP", "triphosphane"),
    ("CSC", "methylsulfanylmethane"),    # 硫醚：标准价且端原子带碳 → 不是杂原子链
]


@pytest.mark.parametrize("smiles,gold", CHAIN_CASES)
def test_heterane_chain_saturated(smiles, gold):
    """均一无环母体氢化物（P-21.2.2）。"""
    assert n.name(smiles).en == gold


CHAIN_UNSAT_CASES = [
    ("N=N", "diazene"),                  # P-14.3.4.2(d)：二核单不饱和省位次
    ("NN=N", "triazene"),                # 三核同理
    ("NN=NNN", "pentaaz-2-ene"),         # 规则原文例句：位次依杂链编号取最低
    ("N#N", "diazyne"),
    ("S=S", "disulfene"),
    ("[SiH2]=[SiH2]", "disilene"),       # 规则原文例句
    ("P=P", "diphosphene"),
    ("CC=CC", "but-2-ene"),              # 碳链不饱和不受影响
]


@pytest.mark.parametrize("smiles,gold", CHAIN_UNSAT_CASES)
def test_heterane_chain_unsaturation(smiles, gold):
    """杂原子链的不饱和态与碳同样处理（P-31.1.2.2.1）。"""
    assert n.name(smiles).en == gold


# ── 取代链与支链：与碳链同一套「枚举链候选 → 支链递归成前缀」路径 ──
CHAIN_SUBST_CASES = [
    ("CSSSC", "1,3-dimethyltrisulfane"),      # 取代链：位次按烃的方式从一端编号
    ("COOOC", "1,3-dimethyltrioxidane"),
    ("CNNN", "1-methyltriazane"),             # P-68.3.1.4.1 原文例：氮链的取代基用数字位次，非 N-
    ("C[SiH2][SiH2][SiH3]", "1-methyltrisilane"),
    ("[SiH3][SiH]([SiH3])[SiH2][SiH3]", "2-silyltetrasilane"),   # 支链：最长链作母体，支链为 silyl 前缀
    ("[SiH3][SiH2][SiH2][SiH]([SiH2][SiH3])[SiH2][SiH3]", "3-disilanylhexasilane"),  # 支链本身是链：disilane → disilanyl
]


@pytest.mark.parametrize("smiles,gold", CHAIN_SUBST_CASES)
def test_heterane_chain_substituted_and_branched(smiles, gold):
    """取代链与支链（P-21.2.2 预选名可取代；有支链则最长链作母体、支链作前缀）。"""
    assert n.name(smiles).en == gold


def test_heterane_chain_short_chalcogen_not_parent():
    """连续 1–2 个硫族原子不作母体氢化物（P-68.4.0），仍走既有二硫化物/过氧化物通路。"""
    assert n.name("CSSC").en == "(methyldisulfanyl)methane"
    assert n.name("COOC").en == "methoxyoxymethane"


def test_heterane_chain_excludes_ring_and_charged():
    """环内杂原子与带电杂原子链不算无环母体氢化物（P-21.2.2；叠氮 N=[N+]=[N-]）。"""
    assert "1,2,3-triazol-1-yl" in n.name("C1=CN(N=N1)CCCNC2=C(C=C(C=N2)Br)F").en  # 环内 N-N-N 不得被当链
    assert n.name("NNCC").en == "hydrazinylethane"


# ── 环系 λ（P-22.2.7 / P-23.6 / P-25.6）──────────
RING_LAMBDA_CASES = [
    ("O=P1(NCCCl)OCCCN1", "1,3,2lambda5-oxazaphosphinane"),          # P-22.2.7.1：λ 紧跟杂原子位次
    ("O=S1(=O)N2CN3CCN(C2)CN1C3", "9lambda6-thia-"),                 # P-23.6.1：λ 置于 'a' 前缀之前
    ("N(=[N+]=[N-])I1OC(C2=C1C=CC=C2)=O", "1-azido-1lambda3-"),      # P-25.6：λ 在稠环系统名首
]


@pytest.mark.parametrize("smiles,fragment", RING_LAMBDA_CASES)
def test_ring_lambda_present(smiles, fragment):
    """环系非标准键数标 λ，且各自落在规则要求的位置。"""
    assert fragment in n.name(smiles).en


def test_ring_lambda_absent_for_standard_valence():
    """标准价杂环不带 λ（噻吩/吡啶）。"""
    assert "lambda" not in n.name("c1ccsc1").en
    assert "lambda" not in n.name("c1ccncc1").en


def test_heterane_chain_lambda():
    """链上非标准键数原子带 λ 位次（P-21.2.4：λn 置于各位次之后）。"""
    assert n.name("[SH2](S)S").en == "2lambda4-trisulfane"   # P-21.2.4 原文例
    assert n.name("SSSS").en == "tetrasulfane"               # 标准价链不带 λ


def test_heterane_chain_excludes_non_parent_hydrides():
    """不是母体氢化物的同元素相邻原子不算链：过氧 -O-OH、二硫化物 R-S-S-R。"""
    assert "oxidane" not in n.name("CC(=O)OO").en        # 过氧乙酸走含氧酸通路
    assert "sulfane" not in n.name("CSSC").en            # 二甲基二硫醚：端原子带碳
