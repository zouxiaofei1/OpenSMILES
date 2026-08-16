"""简单取代基的锚定 canonical-SMILES 查表：锚定键 → (en, zh, paren, kind)。
未命中回退完整命名路径；保留取代基经 anchored 字段反查索引，非保留取代基
内联存储；构建期校验 canonical 对拍、键唯一性与跨表冲突。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.layer3.submol_build import build_anchor_submol


# 非保留取代基：锚定 SMILES → (en, zh, requires_parentheses, kind)
# 保留取代基（registry 条目）的锚定键在 RetainedSubstituent.anchored 中声明。
ANCHOR_TABLE: dict[str, tuple[str, str, bool, str]] = {
    #基础取代基
    "*F": ("fluoro", "氟", False, "halo"),
    "*Cl": ("chloro", "氯", False, "halo"),
    "*Br": ("bromo", "溴", False, "halo"),
    "*I": ("iodo", "碘", False, "halo"),
    "*[N+](=O)[O-]": ("nitro", "硝基", False, "leaf"),
    "*N=C=O": ("isocyanato", "异氰酸根合", False, "leaf"),
    "*N=C=S": ("isothiocyanato", "异硫氰酸根合", False, "leaf"),
    "*C": ("methyl", "甲基", False, "alkyl"),
    "*CC": ("ethyl", "乙基", False, "alkyl"),
    "*CCC": ("propyl", "丙基", False, "alkyl"),
    "*CCCC": ("butyl", "丁基", False, "alkyl"),
    "*CCCCC": ("pentyl", "戊基", False, "alkyl"),
    "*CCCCCC": ("hexyl", "己基", False, "alkyl"),
    "*CCCCCCC": ("heptyl", "庚基", False, "alkyl"),
    "*CCCCCCCC": ("octyl", "辛基", False, "alkyl"),
    "*CCCCCCCCC": ("nonyl", "壬基", False, "alkyl"),
    "*CCCCCCCCCC": ("decyl", "癸基", False, "alkyl"),
    "*CCCCCCCCCCC": ("undecyl", "十一烷基", False, "alkyl"),
    "*C(F)(F)F": ("trifluoromethyl", "三氟甲基", False, "alkyl"),
    "*C1CC1": ("cyclopropyl", "环丙基", False, "alkyl"),
    "*C1CCC1": ("cyclobutyl", "环丁基", False, "alkyl"),
    "*C1CCCC1": ("cyclopentyl", "环戊基", False, "alkyl"),
    "*C1CCCCC1": ("cyclohexyl", "环己基", False, "alkyl"),
    "*C1CCCCCC1": ("cycloheptyl", "环庚基", False, "alkyl"),
    "*C1CCCCCCC1": ("cyclooctyl", "环辛基", False, "alkyl"),
    "*c1ccccc1": ("phenyl", "苯基", False, "aryl"),
    #因未知原因删了会regress的取代基
    "*CCCl": ("2-chloroethyl", "2-氯乙基", True, "alkyl"),
    "*CCCCl": ("3-chloropropyl", "3-氯丙基", True, "alkyl"),
    "*CCCCCl": ("4-chlorobutyl", "4-氯丁基", True, "alkyl"),
    "*CCCBr": ("3-bromopropyl", "3-溴丙基", True, "alkyl"),
    "*CCCCBr": ("4-bromobutyl", "4-溴丁基", True, "alkyl"),
    "*c1cccnc1": ("pyridin-3-yl", "吡啶-3-基", True, "aryl"),
    "*c1ccncc1": ("pyridin-4-yl", "吡啶-4-基", True, "aryl"),
    "*c1ccco1": ("furan-3-yl", "呋喃-3-基", True, "aryl"),
    "*c1cccs1": ("thiophen-3-yl", "噻吩-3-基", True, "aryl"),
    "*c1ccc2cccnc2c1": ("quinolin-3-yl", "喹啉-3-基", True, "aryl"),
    "*c1cccc2cccnc12": ("quinolin-4-yl", "喹啉-4-基", True, "aryl"),
    "*c1cc2ccccc2cn1": ("isoquinolin-3-yl", "异喹啉-3-基", True, "aryl"),
    "*c1cncc2ccccc12": ("isoquinolin-4-yl", "异喹啉-4-基", True, "aryl"),
    "*c1c2ccccc2cc2ccccc12": ("anthracen-9-yl", "蒽-9-基", True, "aryl"),
    "*c1ccc2ccc3ccccc3c2c1": ("phenanthren-1-yl", "菲-1-基", True, "aryl"),
    "*c1ccc2c(ccc3ccccc32)c1": ("phenanthren-2-yl", "菲-2-基", True, "aryl"),
    "*C1CCNCC1": ("piperidin-4-yl", "哌啶-4-基", True, "alkyl"),
    "*C1CCCNC1": ("piperidin-3-yl", "哌啶-3-基", True, "alkyl"),
    "*[C@@H]1CCCCN1": ("piperidin-2-yl", "哌啶-2-基", True, "alkyl"),
    "*CCl": ("chloromethyl", "氯甲基", True, "alkyl"),
    "*CBr": ("bromomethyl", "溴甲基", True, "alkyl"),
    #因iupac规则保留的取代基
}


# ── 保留取代基注册表（IUPAC 2013 P-29/P-57；general/pin 双语名称）──

class IupacLevel(Enum):
    """IUPAC 命名层级：PIN / general / 不推荐。"""

    PIN = "pin"
    GENERAL = "general"
    NOT_RECOMMENDED = "not_rec"


@dataclass(frozen=True)
class RetainedSubstituent:
    """一个保留取代基条目的双语名称、层级与锚定键。"""

    en: str
    zh: str
    systematic_en: str
    systematic_zh: str
    level: IupacLevel
    anchored: tuple[str, ...] = ()  # 本保留基对应的锚定 canonical-SMILES 键
    paren: bool = False  # 作前缀时是否需要括号
    kind: str = "leaf"  # alkyl/aryl/halo/leaf


def _build_registry() -> dict[str, RetainedSubstituent]:
    """构建保留取代基注册表（键 → RetainedSubstituent 条目）。"""
    P, G, N = IupacLevel.PIN, IupacLevel.GENERAL, IupacLevel.NOT_RECOMMENDED
    return {
        # 支链 / 不饱和烷基
        "tert-butyl": RetainedSubstituent( "tert-butyl", "叔丁基", "1,1-dimethylethyl", "1,1-二甲基乙基", P, anchored=("*C(C)(C)C", ), paren=False, kind="alkyl", ),
        "isopropyl": RetainedSubstituent( "isopropyl", "异丙基", "propan-2-yl", "丙-2-基", G, anchored=("*C(C)C", ), paren=False, kind="alkyl", ),
        "isobutyl": RetainedSubstituent("isobutyl", "异丁基", "2-methylpropyl", "2-甲基丙基", N, anchored=("*CC(C)C", ), paren=False, kind="alkyl", ),
        "sec-butyl": RetainedSubstituent("sec-butyl", "仲丁基", "butan-2-yl", "丁-2-基", N, anchored=("*C(C)CC", ), paren=False, kind="alkyl", ),
        "neopentyl": RetainedSubstituent( "neopentyl", "新戊基", "2,2-dimethylpropyl", "2,2-二甲基丙基", N, anchored=("*CC(C)(C)C", ), paren=False, kind="alkyl", ),
        "isopentyl": RetainedSubstituent( "isopentyl", "异戊基", "3-methylbutyl", "3-甲基丁基", N, anchored=("*CCC(C)C", ), paren=False, kind="alkyl", ),
        "2-methylbutan-2-yl": RetainedSubstituent("2-methylbutan-2-yl", "2-甲基丁-2-基", "2-methylbutan-2-yl", "2-甲基丁-2-基", P, anchored=("*C(C)(C)CC", ), paren=False, kind="alkyl", ),
        "3-methylbut-2-enyl": RetainedSubstituent("3-methylbut-2-enyl", "3-甲基丁-2-烯基", "3-methylbut-2-enyl", "3-甲基丁-2-烯基", P, anchored=("*CC=C(C)C", ), paren=False, kind="alkyl", ),
        "vinyl": RetainedSubstituent( "vinyl", "乙烯基", "ethenyl", "乙烯基", G, anchored=("*C=C", ), paren=False, kind="alkyl", ),
        "allyl": RetainedSubstituent( "allyl", "烯丙基", "prop-2-en-1-yl", "丙-2-烯-1-基", G, anchored=("*CC=C", ), paren=False, kind="alkyl", ),
        "isopropenyl": RetainedSubstituent( "isopropenyl", "异丙烯基", "prop-1-en-2-yl", "丙-1-烯-2-基", G, anchored=("*C(=C)C", ), paren=False, kind="alkyl", ),
        "propargyl": RetainedSubstituent( "propargyl", "炔丙基", "prop-2-yn-1-yl", "丙-2-炔-1-基", G, anchored=("*CC#C", ), paren=False, kind="alkyl", ),
        "difluoromethyl": RetainedSubstituent( "difluoromethyl", "二氟甲基", "difluoromethyl", "二氟甲基", P, anchored=("*C(F)F", ), paren=False, kind="alkyl", ),
        "pentafluoroethyl": RetainedSubstituent( "pentafluoroethyl", "五氟乙基", "pentafluoroethyl", "五氟乙基", P, anchored=("*C(F)(F)C(F)(F)F", ), paren=False, kind="alkyl", ),
        "trichloromethyl": RetainedSubstituent( "trichloromethyl", "三氯甲基", "trichloromethyl", "三氯甲基", P, anchored=("*C(Cl)(Cl)Cl", ), paren=False, kind="alkyl", ),
        "tribromomethyl": RetainedSubstituent( "tribromomethyl", "三溴甲基", "tribromomethyl", "三溴甲基", P, anchored=("*C(Br)(Br)Br", ), paren=False, kind="alkyl", ),
        # 芳烷基 / 芳基
        "phenethyl": RetainedSubstituent( "phenethyl", "2-苯乙基", "2-phenylethyl", "2-苯乙基", G, anchored=("*CCc1ccccc1", ), paren=True, kind="aryl", ),
        "benzhydryl": RetainedSubstituent( "benzhydryl", "二苯甲基", "diphenylmethyl", "二苯甲基", G, anchored=("*C(c1ccccc1)c1ccccc1", ), paren=False, kind="aryl", ),
        "trityl": RetainedSubstituent( "trityl", "三苯甲基", "triphenylmethyl", "三苯甲基", G, anchored=("*C(c1ccccc1)(c1ccccc1)c1ccccc1", ), paren=False, kind="aryl", ),
        "cinnamyl": RetainedSubstituent( "cinnamyl", "肉桂基", "3-phenylprop-2-en-1-yl", "3-苯基丙-2-烯-1-基", G, anchored=("*C=Cc1ccccc1", ), paren=False, kind="aryl", ),
        "benzyl": RetainedSubstituent( "benzyl", "苄基", "phenylmethyl", "苯甲基", P, anchored=("*Cc1ccccc1", ), paren=False, kind="aryl", ),
        "furfuryl": RetainedSubstituent( "furfuryl", "糠基", "furan-2-ylmethyl", "呋喃-2-基甲基", N, anchored=("*Cc1ccoc1", ), paren=True, kind="aryl", ),
        "thenyl": RetainedSubstituent( "thenyl", "噻吩甲基", "thiophen-2-ylmethyl", "噻吩-2-基甲基", N, anchored=("*Cc1ccsc1", ), paren=True, kind="aryl", ),
        "furyl": RetainedSubstituent( "furyl", "呋喃基", "furan-2-yl", "呋喃-2-基", G, anchored=("*c1ccoc1", ), paren=True, kind="aryl", ),
        "thienyl": RetainedSubstituent( "thienyl", "噻吩基", "thiophen-2-yl", "噻吩-2-基", G, anchored=("*c1ccsc1", ), paren=True, kind="aryl", ),
        "pyridyl": RetainedSubstituent( "pyridyl", "吡啶基", "pyridin-2-yl", "吡啶-2-基", G, anchored=("*c1ccccn1", ), paren=True, kind="aryl", ),
        "quinolyl": RetainedSubstituent( "quinolyl", "喹啉基", "quinolin-2-yl", "喹啉-2-基", G, anchored=("*c1ccc2ncccc2c1", ), paren=True, kind="aryl", ),
        "isoquinolyl": RetainedSubstituent( "isoquinolyl", "异喹啉基", "isoquinolin-1-yl", "异喹啉-1-基", G, anchored=("*c1nccc2ccccc12", ), paren=True, kind="aryl", ),
        "anthryl": RetainedSubstituent( "anthryl", "蒽基", "anthracen-2-yl", "蒽-2-基", G, anchored=("*c1ccc2cc3ccccc3cc2c1", ), paren=True, kind="aryl", ),
        "phenanthryl": RetainedSubstituent( "phenanthryl", "菲基", "phenanthren-9-yl", "菲-9-基", G, anchored=("*c1cc2ccccc2c2ccccc12", ), paren=True, kind="aryl", ),
        "adamantyl": RetainedSubstituent( "adamantyl", "金刚烷基", "adamantan-2-yl", "金刚烷-2-基", G, anchored=("*C1C2CC3CC(C2)CC1C3", ), paren=True, kind="aryl", ),
        # 杂原子前缀 — O / S
        "methoxy": RetainedSubstituent( "methoxy", "甲氧基", "methoxy", "甲氧基", P, anchored=("*OC", ), paren=False, kind="leaf", ),
        "hydroperoxy": RetainedSubstituent( "hydroperoxy", "氢过氧基", "hydroperoxy", "氢过氧基", P, anchored=("*OO", ), paren=False, kind="leaf", ),
        "methylsulfanyl": RetainedSubstituent( "methylsulfanyl", "甲硫基", "methylsulfanyl", "甲硫基", P, anchored=("*SC", ), paren=False, kind="leaf", ),
        "ethylsulfanyl": RetainedSubstituent( "ethylsulfanyl", "乙硫基", "ethylsulfanyl", "乙硫基", P, anchored=("*SCC", ), paren=False, kind="leaf", ),
        "sulfanyl": RetainedSubstituent( "sulfanyl", "硫烷基", "sulfanyl", "硫烷基", P, anchored=("*S", ), paren=False, kind="leaf", ),
        "selanyl": RetainedSubstituent( "selanyl", "硒烷基", "selanyl", "硒烷基", P, anchored=("*[SeH]", ), paren=False, kind="leaf", ),
        "methylsulfinyl": RetainedSubstituent( "methylsulfinyl", "甲亚磺酰基", "methanesulfinyl", "甲亚磺酰基", P, anchored=("*S(C)=O", ), paren=False, kind="leaf", ),
        "methylsulfonyl": RetainedSubstituent( "methylsulfonyl", "甲磺酰基", "methanesulfonyl", "甲磺酰基", P, anchored=("*S(C)(=O)=O", ), paren=False, kind="leaf", ),
        "sulfo": RetainedSubstituent( "sulfo", "磺基", "sulfo", "磺基", P, anchored=("*S(=O)(=O)O", ), paren=False, kind="leaf", ),
        "tosyl": RetainedSubstituent( "tosyl", "对甲苯磺酰基", "4-methylbenzenesulfonyl", "4-甲基苯磺酰基", N, anchored=("*S(=O)(=O)c1ccc(C)cc1", ), paren=False, kind="leaf", ),
        "triflyl": RetainedSubstituent( "triflyl", "三氟甲磺酰基", "trifluoromethanesulfonyl", "三氟甲磺酰基", N, anchored=("*S(=O)(=O)C(F)(F)F", ), paren=False, kind="leaf", ),
        # 杂原子前缀 — N
        "nitroso": RetainedSubstituent( "nitroso", "亚硝基", "nitroso", "亚硝基", P, anchored=("*N=O", ), paren=False, kind="leaf", ),
        "azido": RetainedSubstituent( "azido", "叠氮基", "azido", "叠氮基", P, anchored=("*N=[N+]=[N-]", ), paren=False, kind="leaf", ),
        "hydrazinyl": RetainedSubstituent( "hydrazinyl", "肼基", "hydrazinyl", "肼基", P, anchored=("*NN", ), paren=False, kind="leaf", ),
        "anilino": RetainedSubstituent( "anilino", "苯胺基", "phenylamino", "苯氨基", P, anchored=("*Nc1ccccc1", ), paren=False, kind="leaf", ),
        "methylamino": RetainedSubstituent( "methylamino", "甲氨基", "methylamino", "甲氨基", P, anchored=("*NC", ), paren=False, kind="leaf", ),
        "dimethylamino": RetainedSubstituent( "dimethylamino", "二甲氨基", "dimethylamino", "二甲氨基", P, anchored=("*N(C)C", ), paren=False, kind="leaf", ),
        "diazenyl": RetainedSubstituent( "diazenyl", "二氮烯基", "diazenyl", "二氮烯基", P, anchored=("*N=N", ), paren=False, kind="leaf", ),
        "diazo": RetainedSubstituent( "diazo", "重氮基", "diazo", "重氮基", P, anchored=("*[N+]=[N-]", ), paren=False, kind="leaf", ),
        "cyano": RetainedSubstituent( "cyano", "氰基", "cyano", "氰基", P, anchored=("*C#N", ), paren=False, kind="leaf", ),
        "isocyano": RetainedSubstituent( "isocyano", "异氰基", "isocyano", "异氰基", P, anchored=("*[N+]#[C-]", ), paren=False, kind="leaf", ),
        "carbamoyl": RetainedSubstituent( "carbamoyl", "氨基甲酰基", "aminocarbonyl", "氨基羰基", P, anchored=("*C(N)=O", ), paren=False, kind="leaf", ),
        "sulfamoyl": RetainedSubstituent( "sulfamoyl", "氨基磺酰基", "aminosulfonyl", "氨基磺酰基", P, anchored=("*S(N)(=O)=O", ), paren=False, kind="leaf", ),
        "amidino": RetainedSubstituent( "amidino", "脒基", "carbaminidoyl", "甲脒基", N, anchored=("*C(=N)N", ), paren=False, kind="leaf", ),
        "guanidino": RetainedSubstituent( "guanidino", "胍基", "carbamimidamido", "胍基", N, anchored=("*NC(=N)N", ), paren=False, kind="leaf", ),
        "ureido": RetainedSubstituent( "ureido", "脲基", "carbamoylamino", "脲基", N, anchored=("*NC(N)=O", ), paren=False, kind="leaf", ),
        # P / B / Se 前缀
        "phosphanyl": RetainedSubstituent( "phosphanyl", "膦基", "phosphanyl", "膦基", P, anchored=("*P", ), paren=False, kind="leaf", ),
        "boranyl": RetainedSubstituent( "boranyl", "硼烷基", "boranyl", "硼烷基", P, anchored=("*B", ), paren=False, kind="leaf", ),
        # 酰基前缀
        "formyl": RetainedSubstituent( "formyl", "甲酰基", "formyl", "甲酰基", P, anchored=("*C=O", ), paren=False, kind="leaf", ),
        "acetyl": RetainedSubstituent( "acetyl", "乙酰基", "acetyl", "乙酰基", P, anchored=("*C(C)=O", ), paren=False, kind="leaf", ),
        "benzoyl": RetainedSubstituent( "benzoyl", "苯甲酰基", "benzoyl", "苯甲酰基", P, anchored=("*C(=O)c1ccccc1", ), paren=False, kind="leaf", ),
        "propionyl": RetainedSubstituent( "propionyl", "丙酰基", "propanoyl", "丙酰基", G, anchored=("*C(=O)CC", ), paren=False, kind="leaf", ),
        "butyryl": RetainedSubstituent( "butyryl", "丁酰基", "butanoyl", "丁酰基", G, anchored=("*C(=O)CCC", ), paren=False, kind="leaf", ),
        "isobutyryl": RetainedSubstituent( "isobutyryl", "异丁酰基", "2-methylpropanoyl", "2-甲基丙酰基", G, anchored=("*C(=O)C(C)C", ), paren=False, kind="leaf", ),
        "valeryl": RetainedSubstituent( "valeryl", "戊酰基", "pentanoyl", "戊酰基", G, anchored=("*C(=O)CCCC", ), paren=False, kind="leaf", ),
        "oxamoyl": RetainedSubstituent( "oxamoyl", "草氨酰基", "oxamoyl", "草氨酰基", P, anchored=("*C(=O)C(N)=O", ), paren=False, kind="leaf", ),
        "methoxycarbonyl": RetainedSubstituent( "methoxycarbonyl", "甲氧羰基", "methoxycarbonyl", "甲氧羰基", P, anchored=("*C(=O)OC", ), paren=False, kind="leaf", ),
        "ethoxycarbonyl": RetainedSubstituent( "ethoxycarbonyl", "乙氧羰基", "ethoxycarbonyl", "乙氧羰基", P, anchored=("*C(=O)OCC", ), paren=False, kind="leaf", ),
    }


_REGISTRY: dict[str, RetainedSubstituent] = _build_registry()


def _canon(smi: str) -> str:
    """锚定 SMILES 的 canonical 形式；解析失败抛错。"""
    from rdkit.Chem import MolFromSmiles, MolToSmiles

    m = MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"[anchored_table] 无法解析锚定键: {smi!r}")
    return MolToSmiles(m)


def _build_anchor_index() -> dict[str, str]:
    """从 registry 反查锚定 canonical-SMILES → registry key。

    校验：anchored 键必须已是 canonical 形式（防 *C=C-C 类死条目）、
    跨 key 唯一（一个锚定键只映射一个保留基）。"""
    index: dict[str, str] = {}
    for key, entry in _REGISTRY.items():
        for smi in entry.anchored:
            canon = _canon(smi)
            if canon != smi:
                raise ValueError(f"[anchored_table] 非 canonical 锚定键 {key!r}: {smi!r} -> {canon!r}")
            prev = index.get(smi)
            if prev is not None:
                raise ValueError(f"[anchored_table] 锚定键冲突 {smi!r}: {prev!r} vs {key!r}")
            index[smi] = key
    return index


_ANCHOR_INDEX: dict[str, str] = _build_anchor_index()


def _validate_inline_table() -> None:
    """内联表键必须 canonical 且不与 registry 锚定冲突。"""
    for smi in ANCHOR_TABLE:
        canon = _canon(smi)
        if canon != smi:
            raise ValueError(f"[anchored_table] 内联键非 canonical: {smi!r} -> {canon!r}")
        if smi in _ANCHOR_INDEX:
            raise ValueError(f"[anchored_table] 内联键与 registry 锚定冲突: {smi!r}")


_validate_inline_table()


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
    """生成原子集的锚定 canonical-SMILES 键，失败时返回 None。"""
    from rdkit.Chem import MolToSmiles

    if attach_old is None:
        attach_old = pick_root(mol, atoms)
    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> tuple[str, str | tuple[str, str, bool, str]] | None:
    """查 registry 锚定索引，再查内联表；返回 (key, payload)。"""
    key = anchored_key(mol, atoms, attach_old)
    if key is None:
        return None
    reg_key = _ANCHOR_INDEX.get(key)
    if reg_key is not None:
        return (key, reg_key)
    hit = ANCHOR_TABLE.get(key)
    return None if hit is None else (key, hit)


def _resolve(
    hit: str | tuple[str, str, bool, str], *, name_mode: str,
) -> tuple[str, str, bool, str]:
    """将查表命中解析为 (en, zh, paren, kind)（区分 registry 键与内联条目）。"""
    if isinstance(hit, str):  # registry key → resolve_name + 条目自带 paren/kind
        en, zh = resolve_name(hit, name_mode=name_mode)
        e = _REGISTRY[hit]
        return en, zh, e.paren, e.kind
    en, zh, paren, kind = hit  # 内联条目
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


def anchored_whole_mol(mol: Mol, *, name_mode: str = "general") -> tuple[str, str, bool, str] | None:
    """整分子 canonical SMILES 命中锚定表时返回 (en, zh, paren, kind)。

    供顶层命名器入口使用：当输入分子本身就是带 * 锚点的一个锚定键
    （如 *C(C)C、*OC、*Cc1ccccc1）时直接返回保留名，避免自由基母体
    管线对表内基团的编号/骨架误认。普通分子（锚定表键均带 * 前缀）
    永不命中，安全走完整命名路径。
    """
    from rdkit.Chem import MolToSmiles

    key = MolToSmiles(mol)
    hit: str | tuple[str, str, bool, str] | None
    reg_key = _ANCHOR_INDEX.get(key)
    if reg_key is not None:
        hit = reg_key
    else:
        hit = ANCHOR_TABLE.get(key)
    return None if hit is None else _resolve(hit, name_mode=name_mode)
