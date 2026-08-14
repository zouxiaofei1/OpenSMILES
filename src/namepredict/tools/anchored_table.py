"""简单取代基的锚定 canonical-SMILES 查表。
键映射到 (registry_key, en, zh, paren, kind) 元组；命中即证明锚定键唯一标识简单取代基，未命中则回退到完整命名路径。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.layer3.submol_build import build_anchor_submol


# 锚定 SMILES → (registry_key | None, en, zh, requires_parentheses, kind)
ANCHOR_TABLE: dict[str, tuple[str | None, str, str, bool, str]] = {
    # 直链正烷基（C1–C12；对 name_mode 不敏感）
    "*C": (None, "methyl", "甲基", False, "alkyl"),
    "*CC": (None, "ethyl", "乙基", False, "alkyl"),
    "*CCC": (None, "propyl", "丙基", False, "alkyl"),
    "*CCCC": (None, "butyl", "丁基", False, "alkyl"),
    "*CCCCC": (None, "pentyl", "戊基", False, "alkyl"),
    "*CCCCCC": (None, "hexyl", "己基", False, "alkyl"),
    "*CCCCCCC": (None, "heptyl", "庚基", False, "alkyl"),
    "*CCCCCCCC": (None, "octyl", "辛基", False, "alkyl"),
    "*CCCCCCCCC": (None, "nonyl", "壬基", False, "alkyl"),
    "*CCCCCCCCCC": (None, "decyl", "癸基", False, "alkyl"),
    "*CCCCCCCCCCC": (None, "undecyl", "十一烷基", False, "alkyl"),
    # 支链保留基（registry 键，对 name_mode 敏感）
    "*C(C)C": ("isopropyl", None, None, False, "alkyl"),
    "*C(C)(C)C": ("tert-butyl", None, None, False, "alkyl"),
    "*CC(C)C": ("isobutyl", None, None, False, "alkyl"),
    "*C(C)CC": ("sec-butyl", None, None, False, "alkyl"),
    "*CC(C)(C)C": ("neopentyl", None, None, False, "alkyl"),
    "*CCC(C)C": ("isopentyl", None, None, False, "alkyl"),
    "*C(C)(C)CC": ("2-methylbutan-2-yl", None, None, False, "alkyl"),
    # 烯基保留基
    "*C=C": ("vinyl", None, None, False, "alkyl"),
    "*C=C-C": ("allyl", None, None, False, "alkyl"),
    "*C=CC": ("allyl", None, None, False, "alkyl"),
    "*CC=C": ("allyl", None, None, False, "alkyl"),
    "*C(=C)C": ("isopropenyl", None, None, False, "alkyl"),
    # 卤代烷基（复合前缀需括号）
    "*CCl": (None, "chloromethyl", "氯甲基", True, "alkyl"),
    "*CBr": (None, "bromomethyl", "溴甲基", True, "alkyl"),
    "*CCCl": (None, "2-chloroethyl", "2-氯乙基", True, "alkyl"),
    "*CCCCl": (None, "3-chloropropyl", "3-氯丙基", True, "alkyl"),
    "*CCCCCl": (None, "4-chlorobutyl", "4-氯丁基", True, "alkyl"),
    "*CCCBr": (None, "3-bromopropyl", "3-溴丙基", True, "alkyl"),
    "*CCCCBr": (None, "4-bromobutyl", "4-溴丁基", True, "alkyl"),
    # 1-环烷基乙基（C1(ring)-C(C)H-）
    "*C(C)C1CCCCC1": (None, "1-cyclohexylethyl", "1-环己基乙基", True, "alkyl"),
    # CF3（三氟甲基）
    "*C(F)(F)F": (None, "trifluoromethyl", "三氟甲基", False, "alkyl"),
    # 杂原子前缀 — 烷氧基 / 硫基 / 亚磺酰基 / 磺酰基（registry 键）
    "*OC": ("methoxy", None, None, False, "leaf"),
    "*SC": ("methylsulfanyl", None, None, False, "leaf"),
    "*S(C)=O": ("methylsulfinyl", None, None, False, "leaf"),
    "*S(C)(=O)=O": ("methylsulfonyl", None, None, False, "leaf"),
    "*S(=O)(=O)O": ("sulfo", None, None, False, "leaf"),
    "*S(=O)(=O)c1ccc(C)cc1": ("tosyl", None, None, False, "leaf"),
    "*S(=O)(=O)C(F)(F)F": ("triflyl", None, None, False, "leaf"),
    # 含氮端基官能团（registry 键）
    "*N=O": ("nitroso", None, None, False, "leaf"),
    "*N=[N+]=[N-]": ("azido", None, None, False, "leaf"),
    "*[N+]#[C-]": ("isocyano", None, None, False, "leaf"),
    "*C#N": ("cyano", None, None, False, "leaf"),
    # 哌啶基
    "*C1CCNCC1": (None, "piperidin-4-yl", "哌啶-4-基", True, "alkyl"),
    "*C1CCCNC1": (None, "piperidin-3-yl", "哌啶-3-基", True, "alkyl"),
    "*[C@@H]1CCCCN1": (None, "piperidin-2-yl", "哌啶-2-基", True, "alkyl"),
    # 烯基支链保留基
    "*C=C(C)C": ("isobutyl", None, None, False, "alkyl"),
    "*C=CC(C)C": ("isopentyl", None, None, False, "alkyl"),
    "*CC=C(C)C": ("3-methylbut-2-enyl", None, None, False, "alkyl"),
    # 环烷基（对 name_mode 不敏感）
    "*C1CC1": (None, "cyclopropyl", "环丙基", False, "alkyl"),
    "*C1CCC1": (None, "cyclobutyl", "环丁基", False, "alkyl"),
    "*C1CCCC1": (None, "cyclopentyl", "环戊基", False, "alkyl"),
    "*C1CCCCC1": (None, "cyclohexyl", "环己基", False, "alkyl"),
    "*C1CCCCCC1": (None, "cycloheptyl", "环庚基", False, "alkyl"),
    "*C1CCCCCCC1": (None, "cyclooctyl", "环辛基", False, "alkyl"),
    # 芳基
    "*c1ccccc1": (None, "phenyl", "苯基", False, "aryl"),
    "*c1ccc(Cl)cc1": (None, "4-chlorophenyl", "4-氯苯基", True, "aryl"),
    "*c1cccc(Cl)c1": (None, "3-chlorophenyl", "3-氯苯基", True, "aryl"),
    "*c1ccccc1Cl": (None, "2-chlorophenyl", "2-氯苯基", True, "aryl"),
    # 单原子卤素（始终单键；无键型歧义）
    "*F": (None, "fluoro", "氟", False, "halo"),
    "*Cl": (None, "chloro", "氯", False, "halo"),
    "*Br": (None, "bromo", "溴", False, "halo"),
    "*I": (None, "iodo", "碘", False, "halo"),
    # 键拓扑唯一的多原子端基官能团
    "*[N+](=O)[O-]": (None, "nitro", "硝基", False, "leaf"),
    "*N=C=O": (None, "isocyanato", "异氰酸根合", False, "leaf"),
    "*N=C=S": (None, "isothiocyanato", "异硫氰酸根合", False, "leaf"),
    # ── 保留的芳烷基 / 不饱和基（registry 键）──
    "*Cc1ccccc1": ("benzyl", None, None, False, "aryl"),
    "*CCc1ccccc1": ("phenethyl", None, None, True, "aryl"),
    "*C(c1ccccc1)c1ccccc1": ("benzhydryl", None, None, False, "aryl"),
    "*C(c1ccccc1)(c1ccccc1)c1ccccc1": ("trityl", None, None, False, "aryl"),
    "*CC#C": ("propargyl", None, None, False, "alkyl"),
    "*C=Cc1ccccc1": ("cinnamyl", None, None, False, "aryl"),
    "*Cc1ccoc1": ("furfuryl", None, None, True, "aryl"),
    "*Cc1ccsc1": ("thenyl", None, None, True, "aryl"),
    # ── 杂芳基（在注册位置以 registry 键索引）──
    "*c1ccoc1": ("furyl", None, None, True, "aryl"),
    "*c1ccsc1": ("thienyl", None, None, True, "aryl"),
    "*c1ccccn1": ("pyridyl", None, None, True, "aryl"),
    "*c1ccc2ncccc2c1": ("quinolyl", None, None, True, "aryl"),
    "*c1nccc2ccccc12": ("isoquinolyl", None, None, True, "aryl"),
    "*c1ccc2cc3ccccc3cc2c1": ("anthryl", None, None, True, "aryl"),
    "*c1cc2ccccc2c2ccccc12": ("phenanthryl", None, None, True, "aryl"),
    "*C1C2CC3CC(C2)CC1C3": ("adamantyl", None, None, True, "aryl"),
    # ── 杂芳基其他位置（非注册，描述性名称）──
    "*c1cccnc1": (None, "pyridin-3-yl", "吡啶-3-基", True, "aryl"),
    "*c1ccncc1": (None, "pyridin-4-yl", "吡啶-4-基", True, "aryl"),
    "*c1ccco1": (None, "furan-3-yl", "呋喃-3-基", True, "aryl"),
    "*c1cccs1": (None, "thiophen-3-yl", "噻吩-3-基", True, "aryl"),
    "*c1ccc2cccnc2c1": (None, "quinolin-3-yl", "喹啉-3-基", True, "aryl"),
    "*c1cccc2cccnc12": (None, "quinolin-4-yl", "喹啉-4-基", True, "aryl"),
    "*c1cc2ccccc2cn1": (None, "isoquinolin-3-yl", "异喹啉-3-基", True, "aryl"),
    "*c1cncc2ccccc12": (None, "isoquinolin-4-yl", "异喹啉-4-基", True, "aryl"),
    "*c1c2ccccc2cc2ccccc12": (None, "anthracen-9-yl", "蒽-9-基", True, "aryl"),
    "*c1ccc2ccc3ccccc3c2c1": (None, "phenanthren-1-yl", "菲-1-基", True, "aryl"),
    "*c1ccc2c(ccc3ccccc32)c1": (None, "phenanthren-2-yl", "菲-2-基", True, "aryl"),
    "*C12CC3CC(CC(C3)C1)C2": (None, "adamantan-1-yl", "金刚烷-1-基", True, "aryl"),
    # ── 甲苯基（o/m/p 拆分；NOT_RECOMMENDED，描述性名称）──
    "*c1ccccc1C": (None, "2-methylphenyl", "2-甲基苯基", True, "aryl"),
    "*c1cccc(C)c1": (None, "3-methylphenyl", "3-甲基苯基", True, "aryl"),
    "*c1ccc(C)cc1": (None, "4-methylphenyl", "4-甲基苯基", True, "aryl"),
    # ── 杂原子前缀 — O / S（registry 键）──
    "*OO": ("hydroperoxy", None, None, False, "leaf"),
    "*SCC": ("ethylsulfanyl", None, None, False, "leaf"),
    "*S": ("sulfanyl", None, None, False, "leaf"),
    # ── 杂原子前缀 — N（registry 键）──
    "*NN": ("hydrazinyl", None, None, False, "leaf"),
    "*Nc1ccccc1": ("anilino", None, None, False, "leaf"),
    "*NC": ("methylamino", None, None, False, "leaf"),
    "*N(C)C": ("dimethylamino", None, None, False, "leaf"),
    "*N=N": ("diazenyl", None, None, False, "leaf"),
    "*[N+]=[N-]": ("diazo", None, None, False, "leaf"),
    "*C(N)=O": ("carbamoyl", None, None, False, "leaf"),
    "*S(N)(=O)=O": ("sulfamoyl", None, None, False, "leaf"),
    "*C(=N)N": ("amidino", None, None, False, "leaf"),
    "*NC(=N)N": ("guanidino", None, None, False, "leaf"),
    "*NC(N)=O": ("ureido", None, None, False, "leaf"),
    # ── P / B / Se 前缀（registry 键）──
    "*P": ("phosphanyl", None, None, False, "leaf"),
    "*B": ("boranyl", None, None, False, "leaf"),
    "*[SeH]": ("selanyl", None, None, False, "leaf"),
    # ── 多卤代烷基（registry 键；CF3 型）──
    "*C(Cl)(Cl)Cl": ("trichloromethyl", None, None, False, "alkyl"),
    "*C(Br)(Br)Br": ("tribromomethyl", None, None, False, "alkyl"),
    "*C(F)F": ("difluoromethyl", None, None, False, "alkyl"),
    "*C(F)(F)C(F)(F)F": ("pentafluoroethyl", None, None, False, "alkyl"),
    # ── 酰基前缀（registry 键）──
    "*C=O": ("formyl", None, None, False, "leaf"),
    "*C(C)=O": ("acetyl", None, None, False, "leaf"),
    "*C(=O)c1ccccc1": ("benzoyl", None, None, False, "leaf"),
    "*C(=O)CC": ("propionyl", None, None, False, "leaf"),
    "*C(=O)CCC": ("butyryl", None, None, False, "leaf"),
    "*C(=O)C(C)C": ("isobutyryl", None, None, False, "leaf"),
    "*C(=O)CCCC": ("valeryl", None, None, False, "leaf"),
    "*C(=O)C(N)=O": ("oxamoyl", None, None, False, "leaf"),
    "*C(=O)OC": ("methoxycarbonyl", None, None, False, "leaf"),
    "*C(=O)OCC": ("ethoxycarbonyl", None, None, False, "leaf"),
}


# ── 保留取代基注册表（IUPAC 2013 P-29/P-57；general/pin 双语名称）──

class IupacLevel(Enum):
    PIN = "pin"
    GENERAL = "general"
    NOT_RECOMMENDED = "not_rec"


@dataclass(frozen=True)
class RetainedSubstituent:
    en: str
    zh: str
    systematic_en: str
    systematic_zh: str
    level: IupacLevel


def _build_registry() -> dict[str, RetainedSubstituent]:
    P, G, N = IupacLevel.PIN, IupacLevel.GENERAL, IupacLevel.NOT_RECOMMENDED
    return {
        "tert-butyl": RetainedSubstituent( "tert-butyl", "叔丁基", "1,1-dimethylethyl", "1,1-二甲基乙基", P, ),
        "isopropyl": RetainedSubstituent( "isopropyl", "异丙基", "propan-2-yl", "丙-2-基", G,  ),
        "isobutyl": RetainedSubstituent("isobutyl", "异丁基", "2-methylpropyl", "2-甲基丙基", N, ),
        "sec-butyl": RetainedSubstituent("sec-butyl", "仲丁基", "butan-2-yl", "丁-2-基", N, ),
        "neopentyl": RetainedSubstituent( "neopentyl", "新戊基", "2,2-dimethylpropyl", "2,2-二甲基丙基", N, ),
        "isopentyl": RetainedSubstituent( "isopentyl", "异戊基", "3-methylbutyl", "3-甲基丁基", N,  ),
        "2-methylbutan-2-yl": RetainedSubstituent("2-methylbutan-2-yl", "2-甲基丁-2-基", "2-methylbutan-2-yl", "2-甲基丁-2-基", P,  ),
        "3-methylbut-2-enyl": RetainedSubstituent(   "3-methylbut-2-enyl", "3-甲基丁-2-烯基",  "3-methylbut-2-enyl", "3-甲基丁-2-烯基", P,  ),
        "phenethyl": RetainedSubstituent(  "phenethyl", "2-苯乙基", "2-phenylethyl", "2-苯乙基", G,  ),
        "benzhydryl": RetainedSubstituent(    "benzhydryl", "二苯甲基", "diphenylmethyl", "二苯甲基", G,    ),
        "trityl": RetainedSubstituent(   "trityl", "三苯甲基", "triphenylmethyl", "三苯甲基", G,   ),
        "vinyl": RetainedSubstituent(  "vinyl", "乙烯基", "ethenyl", "乙烯基", G,   ),
        "allyl": RetainedSubstituent(   "allyl", "烯丙基", "prop-2-en-1-yl", "丙-2-烯-1-基", G,   ),
        "isopropenyl": RetainedSubstituent(  "isopropenyl", "异丙烯基", "prop-1-en-2-yl", "丙-1-烯-2-基", G,   ),
        "propargyl": RetainedSubstituent(   "propargyl", "炔丙基", "prop-2-yn-1-yl", "丙-2-炔-1-基", G,  ),
        "crotyl": RetainedSubstituent(   "crotyl", "巴豆基", "(E)-but-2-enyl", "(E)-丁-2-烯基", G,  ),
        "cinnamyl": RetainedSubstituent(   "cinnamyl", "肉桂基", "3-phenylprop-2-en-1-yl", "3-苯基丙-2-烯-1-基", G,  ),
        "phenyl": RetainedSubstituent(  "phenyl", "苯基", "phenyl", "苯基", P,  ),
        "benzyl": RetainedSubstituent(     "benzyl", "苄基", "phenylmethyl", "苯甲基", P, ),
        "tolyl": RetainedSubstituent(  "tolyl", "甲苯基", "methylphenyl", "甲基苯基", N,   ),
        "furfuryl": RetainedSubstituent(   "furfuryl", "糠基", "furan-2-ylmethyl", "呋喃-2-基甲基", N, ),
        "thenyl": RetainedSubstituent(    "thenyl", "噻吩甲基", "thiophen-2-ylmethyl", "噻吩-2-基甲基", N, ),
        "furyl": RetainedSubstituent(   "furyl", "呋喃基", "furan-2-yl", "呋喃-2-基", G, ),
        "thienyl": RetainedSubstituent( "thienyl", "噻吩基", "thiophen-2-yl", "噻吩-2-基", G,  ),
        "pyridyl": RetainedSubstituent(  "pyridyl", "吡啶基", "pyridin-2-yl", "吡啶-2-基", G,  ),
        "quinolyl": RetainedSubstituent(  "quinolyl", "喹啉基", "quinolin-2-yl", "喹啉-2-基", G,  ),
        "isoquinolyl": RetainedSubstituent( "isoquinolyl", "异喹啉基", "isoquinolin-1-yl", "异喹啉-1-基", G,),
        "anthryl": RetainedSubstituent( "anthryl", "蒽基", "anthracen-2-yl", "蒽-2-基", G,  ),
        "phenanthryl": RetainedSubstituent(   "phenanthryl", "菲基", "phenanthren-9-yl", "菲-9-基", G,  ),
        "adamantyl": RetainedSubstituent(  "adamantyl", "金刚烷基", "adamantan-2-yl", "金刚烷-2-基", G,  ),
        "piperidyl": RetainedSubstituent(   "piperidyl", "哌啶基", "piperidin-2-yl", "哌啶-2-基", G,   ),
        "methoxy": RetainedSubstituent(   "methoxy", "甲氧基", "methoxy", "甲氧基", P, ),
        "hydroperoxy": RetainedSubstituent( "hydroperoxy", "氢过氧基", "hydroperoxy", "氢过氧基", P,  ),
        "methylsulfanyl": RetainedSubstituent( "methylsulfanyl", "甲硫基", "methylsulfanyl", "甲硫基", P, ),
        "ethylsulfanyl": RetainedSubstituent(      "ethylsulfanyl", "乙硫基", "ethylsulfanyl", "乙硫基", P,  ),
        "sulfanyl": RetainedSubstituent(     "sulfanyl", "硫烷基", "sulfanyl", "硫烷基", P,   ),
        "selanyl": RetainedSubstituent(     "selanyl", "硒烷基", "selanyl", "硒烷基", P,   ),
        "methylsulfinyl": RetainedSubstituent(     "methylsulfinyl", "甲亚磺酰基", "methanesulfinyl", "甲亚磺酰基", P, ),
        "methylsulfonyl": RetainedSubstituent(   "methylsulfonyl", "甲磺酰基", "methanesulfonyl", "甲磺酰基", P, ),
        "sulfo": RetainedSubstituent(     "sulfo", "磺基", "sulfo", "磺基", P, ),
        "mesyl": RetainedSubstituent(    "mesyl", "甲磺酰基", "methanesulfonyl", "甲磺酰基", N,   ),
        "tosyl": RetainedSubstituent(   "tosyl", "对甲苯磺酰基", "4-methylbenzenesulfonyl", "4-甲基苯磺酰基", N, ),
        "triflyl": RetainedSubstituent(  "triflyl", "三氟甲磺酰基", "trifluoromethanesulfonyl", "三氟甲磺酰基", N, ),
        "nitroso": RetainedSubstituent(  "nitroso", "亚硝基", "nitroso", "亚硝基", P, ),
        "azido": RetainedSubstituent(   "azido", "叠氮基", "azido", "叠氮基", P, ),
        "hydrazinyl": RetainedSubstituent(   "hydrazinyl", "肼基", "hydrazinyl", "肼基", P, ),
        "anilino": RetainedSubstituent(  "anilino", "苯胺基", "phenylamino", "苯氨基", P,  ),
        "methylamino": RetainedSubstituent(     "methylamino", "甲氨基", "methylamino", "甲氨基", P,),
        "dimethylamino": RetainedSubstituent(     "dimethylamino", "二甲氨基", "dimethylamino", "二甲氨基", P, ),
        "diazenyl": RetainedSubstituent(    "diazenyl", "二氮烯基", "diazenyl", "二氮烯基", P,),
        "diazo": RetainedSubstituent(   "diazo", "重氮基", "diazo", "重氮基", P,  ),
        "cyano": RetainedSubstituent( "cyano", "氰基", "cyano", "氰基", P,  ),
        "isocyano": RetainedSubstituent(     "isocyano", "异氰基", "isocyano", "异氰基", P, ),
        "carbamoyl": RetainedSubstituent(    "carbamoyl", "氨基甲酰基", "aminocarbonyl", "氨基羰基", P,  ),
        "sulfamoyl": RetainedSubstituent(     "sulfamoyl", "氨基磺酰基", "aminosulfonyl", "氨基磺酰基", P, ),
        "amidino": RetainedSubstituent(    "amidino", "脒基", "carbaminidoyl", "甲脒基", N,),
        "guanidino": RetainedSubstituent( "guanidino", "胍基", "carbamimidamido", "胍基", N, ),
        "ureido": RetainedSubstituent(     "ureido", "脲基", "carbamoylamino", "脲基", N, ),
        "phosphanyl": RetainedSubstituent(   "phosphanyl", "膦基", "phosphanyl", "膦基", P,  ),
        "trichloromethyl": RetainedSubstituent( "trichloromethyl", "三氯甲基", "trichloromethyl", "三氯甲基", P,),
        "tribromomethyl": RetainedSubstituent(   "tribromomethyl", "三溴甲基", "tribromomethyl", "三溴甲基", P, ),
        "difluoromethyl": RetainedSubstituent(     "difluoromethyl", "二氟甲基", "difluoromethyl", "二氟甲基", P, ),
        "pentafluoroethyl": RetainedSubstituent(  "pentafluoroethyl", "五氟乙基", "pentafluoroethyl", "五氟乙基", P,  ),
        "formyl": RetainedSubstituent(  "formyl", "甲酰基", "formyl", "甲酰基", P, ),
        "acetyl": RetainedSubstituent("acetyl", "乙酰基", "acetyl", "乙酰基", P, ),
        "benzoyl": RetainedSubstituent(  "benzoyl", "苯甲酰基", "benzoyl", "苯甲酰基", P,),
        "propionyl": RetainedSubstituent("propionyl", "丙酰基", "propanoyl", "丙酰基", G, ),
        "butyryl": RetainedSubstituent("butyryl", "丁酰基", "butanoyl", "丁酰基", G, ),
        "isobutyryl": RetainedSubstituent( "isobutyryl", "异丁酰基", "2-methylpropanoyl", "2-甲基丙酰基", G,),
        "valeryl": RetainedSubstituent(  "valeryl", "戊酰基", "pentanoyl", "戊酰基", G, ),
        "oxamoyl": RetainedSubstituent("oxamoyl", "草氨酰基", "oxamoyl", "草氨酰基", P, ),
        "methoxycarbonyl": RetainedSubstituent(  "methoxycarbonyl", "甲氧羰基", "methoxycarbonyl", "甲氧羰基", P, ),
        "ethoxycarbonyl": RetainedSubstituent( "ethoxycarbonyl", "乙氧羰基", "ethoxycarbonyl", "乙氧羰基", P, ),
    }


_REGISTRY: dict[str, RetainedSubstituent] = _build_registry()


def get_retained(key: str) -> RetainedSubstituent | None:
    """按 registry 键查找保留取代基。"""
    return _REGISTRY.get(key)


def resolve_name(key: str, *, name_mode: str = "general") -> tuple[str, str]:
    """返回给定命名模式下 registry 键对应的 (en, zh)。

    - "general"：始终为保留名/常用名
    - "pin"：除非条目本身为 PIN 级，否则用系统名
    """
    entry = _REGISTRY[key]
    if name_mode == "pin" and entry.level != IupacLevel.PIN:
        return entry.systematic_en, entry.systematic_zh
    return entry.en, entry.zh


def pick_root(mol: Mol, atoms: frozenset[int]) -> int:
    """取代基侧键合原子：`atoms` 中与外部重原子（母体）成键的原子。
    回退到最小索引。"""
    for a in atoms:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in atoms:
                return a
    return min(atoms)


def anchored_key(mol: Mol, atoms: frozenset[int], attach_old: int | None = None) -> str | None:
    from rdkit.Chem import MolToSmiles

    if attach_old is None:
        attach_old = pick_root(mol, atoms)
    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> tuple[str, tuple[str | None, str, str, bool, str]] | None:
    key = anchored_key(mol, atoms, attach_old)
    if key is None:
        return None
    hit = ANCHOR_TABLE.get(key)
    return None if hit is None else (key, hit)


def _resolve(
    hit: tuple[str | None, str, str, bool, str], *, name_mode: str,
) -> tuple[str, str, bool, str]:
    reg_key, en, zh, paren, kind = hit
    if reg_key is not None:
        en, zh = resolve_name(reg_key, name_mode=name_mode)
    return en, zh, paren, kind


def anchored_entry(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None, *, name_mode: str = "general",
) -> tuple[str, str, bool, str] | None:
    """在 name_mode 下将原子集解析为 (en, zh, paren, kind)，无命中则返回 None。"""
    got = _table_hit(mol, atoms, attach_old)
    return None if got is None else _resolve(got[1], name_mode=name_mode)


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
    *, name_mode: str = "general",
) -> tuple[str, str, bool] | None:
    """查找取代基原子集；在 name_mode 下返回 (en, zh, paren)。

    registry 键条目通过文件内 resolve_name 解析，因此 pin 模式对 isopropyl 会生成 propan-2-yl 等。无命中时返回 None。
    """
    entry = anchored_entry(mol, atoms, attach_old, name_mode=name_mode)
    return None if entry is None else (entry[0], entry[1], entry[2])
