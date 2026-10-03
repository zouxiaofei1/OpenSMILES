"""取代基锚定 canonical-SMILES 查表，未命中回退完整命名路径；
简单/保留取代基统一存 registry，构建期校验锚定键唯一。"""
from __future__ import annotations

import re
from dataclasses import dataclass

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, N, O, zh_bridge_root
from namepredict.tools import memo
from namepredict.tools.re import SUB_LOCANT_RE
from namepredict.layer3.submol_build import build_anchor_submol


# ── 取代基注册表（IUPAC 2013 P-29/P-57）──

@dataclass(frozen=True)
class RetainedSubstituent:
    """一个保留取代基条目的双语名称与锚定键。"""
    en: str
    zh: str
    anchored: tuple[str, ...] = ()  # 本保留基对应的锚定 canonical-SMILES 键
    paren: bool = False  # 作前缀时是否需要括号

# 叶子节点，均为不可再取代的取代基
_REGISTRY: dict[str, RetainedSubstituent] = {  
    "fluoro": RetainedSubstituent("fluoro", "氟", anchored=("*F",), paren=False),
    "chloro": RetainedSubstituent("chloro", "氯", anchored=("*Cl",), paren=False),
    "bromo": RetainedSubstituent("bromo", "溴", anchored=("*Br",), paren=False),
    "iodo": RetainedSubstituent("iodo", "碘", anchored=("*I",), paren=False),
    "nitro": RetainedSubstituent("nitro", "硝基", anchored=("*[N+](=O)[O-]",), paren=False),
    "oxo": RetainedSubstituent("oxo", "氧代", anchored=("*=O",), paren=False),
    "isocyanato": RetainedSubstituent("isocyanato", "异氰酸基", anchored=("*N=C=O",), paren=False),
    "isothiocyanato": RetainedSubstituent("isothiocyanato", "异硫氰酸基", anchored=("*N=C=S",), paren=False),
    "cyanato": RetainedSubstituent("cyanato", "氰氧基", anchored=("*OC#N",), paren=False),  # P-67.1.4.2：-O-C#N 由 cyanic acid 衍生的取代基前缀
    "thiocyanato": RetainedSubstituent("thiocyanato", "硫氰基", anchored=("*SC#N",), paren=False),  # P-67.1.4.2：-S-C#N 由 thiocyanic acid 衍生的取代基前缀
    "methyl": RetainedSubstituent("methyl", "甲基", anchored=("*C",), paren=False),
    "methylidene": RetainedSubstituent("methylidene", "亚甲基", anchored=("*=C",), paren=False),
    "ethylidene": RetainedSubstituent("ethylidene", "亚乙基", anchored=("*=CC",), paren=False),
    "propylidene": RetainedSubstituent("propylidene", "亚丙基", anchored=("*=CCC",), paren=False),
    "cyclopropylidene": RetainedSubstituent("cyclopropylidene", "环丙亚基", anchored=("*=C1CC1",), paren=False),
    "cyclohexylidene": RetainedSubstituent("cyclohexylidene", "环己亚基", anchored=("*=C1CCCCC1",), paren=False),
    "sulfanylidene": RetainedSubstituent("sulfanylidene", "硫代", anchored=("*=S",), paren=False),  # P-66.1.2：=S 作前缀统一「硫代」（金标 0 处「硫烷亚基」）
    "diaminomethylidene": RetainedSubstituent("diaminomethylidene", "二氨基亚甲基", anchored=("*C(=N)N",), paren=True),  # 脒/胍残基 C(=N)N 按 P-66.1.1 取亚基式
    "tert-butyl": RetainedSubstituent("tert-butyl", "叔丁基", anchored=("*C(C)(C)C",), paren=False),  # 支链 / 不饱和烷基
    "isopropyl": RetainedSubstituent("propan-2-yl", "丙-2-基", anchored=("*C(C)C",), paren=False),
    "isobutyl": RetainedSubstituent("2-methylpropyl", "2-甲基丙基", anchored=("*CC(C)C",), paren=False),
    "sec-butyl": RetainedSubstituent("butan-2-yl", "丁-2-基", anchored=("*C(C)CC",), paren=False),
    "neopentyl": RetainedSubstituent("2,2-dimethylpropyl", "2,2-二甲基丙基", anchored=("*CC(C)(C)C",), paren=False),
    "isopentyl": RetainedSubstituent("3-methylbutyl", "3-甲基丁基", anchored=("*CCC(C)C",), paren=False),
    "vinyl": RetainedSubstituent("ethenyl", "乙烯基", anchored=("*C=C",), paren=False),
    "allyl": RetainedSubstituent("prop-2-enyl", "丙-2-烯基", anchored=("*CC=C",), paren=False),
    "isopropenyl": RetainedSubstituent("prop-1-en-2-yl", "丙-1-烯-2-基", anchored=("*C(=C)C",), paren=False),
    "propargyl": RetainedSubstituent("prop-2-ynyl", "丙-2-炔基", anchored=("*CC#C",), paren=False),
    "benzyl": RetainedSubstituent("benzyl", "苄基", anchored=("*Cc1ccccc1",), paren=False),
    "methoxy": RetainedSubstituent("methoxy", "甲氧基", anchored=("*OC",), paren=False),
    "hydroperoxy": RetainedSubstituent("hydroperoxy", "氢过氧基", anchored=("*OO",), paren=False),
    "hydroxy": RetainedSubstituent("hydroxy", "羟基", anchored=("*O",), paren=False),
    "oxidanyl": RetainedSubstituent("oxidanyl", "氧基", anchored=("*[O]",), paren=False),
    "oxido": RetainedSubstituent("oxido", "氧基", anchored=("*[O-]",), paren=False),  # 去质子酚氧负离子（P-66.1.1.4 阴离子前缀）
    "methylsulfanyl": RetainedSubstituent("methylsulfanyl", "甲硫基", anchored=("*SC",), paren=False),
    "ethylsulfanyl": RetainedSubstituent("ethylsulfanyl", "乙硫基", anchored=("*SCC",), paren=False),
    "sulfanyl": RetainedSubstituent("sulfanyl", "硫基", anchored=("*S",), paren=False),  # -SH 端点与 -S- 桥共用「硫基」，不按取代基语境分化
    "sulfido": RetainedSubstituent("sulfido", "硫代", anchored=("*[S-]",), paren=False),  # 去质子硫醇负离子（P-66.1.1.4 阴离子前缀）
    "selanyl": RetainedSubstituent("selanyl", "氢硒基", anchored=("*[SeH]",), paren=False),
    "methylselanyl": RetainedSubstituent("methylselanyl", "甲硒基", anchored=("*[Se]C",), paren=False),  # P-63.2 同类硫的 methylsulfanyl 对应的硒类似物
    "methylsulfinyl": RetainedSubstituent("methylsulfinyl", "甲基亚磺酰", anchored=("*S(C)=O",), paren=False),
    "methylsulfonyl": RetainedSubstituent("methylsulfonyl", "甲磺酰基", anchored=("*S(C)(=O)=O",), paren=False),
    "methanesulfonamido": RetainedSubstituent("methanesulfonamido", "甲磺酰胺基", anchored=("*NS(C)(=O)=O",), paren=True),  # P-66.1.1.4.3：甲磺酰基的 N-酰基残基作前缀
    "sulfo": RetainedSubstituent("sulfo", "磺基", anchored=("*S(=O)(=O)O",), paren=False),
    "sulfonato": RetainedSubstituent("sulfonato", "磺酸根", anchored=("*S(=O)(=O)[O-]",), paren=False),
    "sulfinato": RetainedSubstituent("sulfinato", "亚磺酸根", anchored=("*S(=O)[O-]",), paren=False),
    "sulfonatooxy": RetainedSubstituent("sulfonatooxy", "磺酸氧基", anchored=("*OS(=O)(=O)[O-]",), paren=False),
    "carboxy": RetainedSubstituent("carboxy", "羧基", anchored=("*C(=O)O",), paren=False),
    "carboxylato": RetainedSubstituent("carboxylato", "羧酸根", anchored=("*C(=O)[O-]",), paren=False),  # 去质子羧基（P-66.1.1.4）
    "carbamoyl": RetainedSubstituent("carbamoyl", "氨基甲酰基", anchored=("*C(N)=O",), paren=False),  # P-66.1.1.4.1 氨基甲酸酰基的保留前缀 carbamoyl
     "carbamoylamino": RetainedSubstituent("carbamoylamino", "氨基甲酰氨基", anchored=("*NC(N)=O",), paren=True),  # P-66.1.1.6 ureido 已不推荐，优选 carbamoylamino
    "azaniumyl": RetainedSubstituent("azaniumyl", "铵基", anchored=("*[NH3+]",), paren=False),
    "methylazaniumyl": RetainedSubstituent("methylazaniumyl", "甲基铵基", anchored=("*[NH2+]C",), paren=False),
    "dimethylazaniumyl": RetainedSubstituent("dimethylazaniumyl", "二甲基铵基", anchored=("*[NH+](C)C",), paren=True),
    "trimethylazaniumyl": RetainedSubstituent("trimethylazaniumyl", "三甲基铵基", anchored=("*[N+](C)(C)C",), paren=True),
    "carbamoyloxy": RetainedSubstituent("carbamoyloxy", "氨基甲酰氧基", anchored=("*OC(N)=O",), paren=False),  # 氨基甲酸 O-酯残基（P-66.1.1.4.1）
    "carbamothioylamino": RetainedSubstituent("carbamothioylamino", "氨基硫代羰基氨基", anchored=("*NC(N)=S",), paren=True),  # 硫代氨基甲酸残基（P-66.1.1.4：carbamothioyl）
    "sulfamoyl": RetainedSubstituent("sulfamoyl", "氨磺酰基", anchored=("*S(N)(=O)=O",), paren=False),  # P-66.1.1.4.2 磺酰胺基，N-取代时基名前移
    "phosphono": RetainedSubstituent("phosphono", "膦酸", anchored=("*P(=O)(O)O",), paren=False),  # P-102 phosphono 为 -PO(OH)2，P 直连母体
    "phosphonato": RetainedSubstituent("phosphonato", "膦酸根", anchored=("*P(=O)([O-])O", "*P(=O)([O-])[O-]"), paren=False),  # P-102：phosphonato 表示 -PO(O-)2（单/双阴离子）
    "phosphonooxy": RetainedSubstituent("phosphonooxy", "膦酸氧基", anchored=("*OP(=O)(O)O",), paren=False),  # P-67.1.5.1 磷酸降级前缀 phosphonooxy
    "phosphonatooxy": RetainedSubstituent("phosphonatooxy", "膦酸氧基", anchored=("*OP(=O)([O-])O", "*OP(=O)([O-])[O-]"), paren=False),
    "phosphonooxymethyl": RetainedSubstituent("phosphonooxymethyl", "膦酸氧甲基", anchored=("*COP(=O)(O)O",), paren=True),
    "phosphonatooxymethyl": RetainedSubstituent("phosphonatooxymethyl", "膦酸氧甲基", anchored=("*COP(=O)([O-])O",), paren=True),
    "borono": RetainedSubstituent("borono", "硼酸", anchored=("*B(O)O",), paren=False),  # P-67.1.4.2：-B(OH)2 预选前缀 borono
    "trimethylsilyl": RetainedSubstituent("trimethylsilyl", "三甲基甲硅烷基", anchored=("*[Si](C)(C)C",), paren=False),  # P-67.1.4.2：-Si(CH3)3 由 silane 衍生
    "difluoroboranyl": RetainedSubstituent("difluoroboranyl", "二氟硼烷基", anchored=("*B(F)F",), paren=False),  # P-68：-BF2 由 borane 衍生的 boranyl 前缀
    "diphenylboranyl": RetainedSubstituent("diphenylboranyl", "二苯基硼烷基", anchored=("*B(c1ccccc1)c1ccccc1",), paren=False),  # P-68 例：bis(diphenylboranyl)
    "formyl": RetainedSubstituent("formyl", "甲酰基", anchored=("*C=O",), paren=False),
    "carboxymethyl": RetainedSubstituent("carboxymethyl", "羧甲基", anchored=("*CC(=O)O",), paren=True),
    "hydroxymethyl": RetainedSubstituent("hydroxymethyl", "羟甲基", anchored=("*CO",), paren=True),
    "nitroso": RetainedSubstituent("nitroso", "亚硝基", anchored=("*N=O",), paren=False),
    "azido": RetainedSubstituent("azido", "叠氮基", anchored=("*N=[N+]=[N-]",), paren=False),
    "amino": RetainedSubstituent("amino", "氨基", anchored=("*N",), paren=False),
    "hydrazinyl": RetainedSubstituent("hydrazinyl", "肼基", anchored=("*NN",), paren=False),
    "anilino": RetainedSubstituent("anilino", "苯胺基", anchored=("*Nc1ccccc1",), paren=False),
    "diazenyl": RetainedSubstituent("diazenyl", "二氮烯基", anchored=("*N=N",), paren=False),
    "diazo": RetainedSubstituent("diazo", "重氮基", anchored=("*=[N+]=[N-]",), paren=False),  # P-66.3：重氮基 C=N+=N-，取代体碳以双键连 N+（原单键锚不可达）
    "diazonio": RetainedSubstituent("diazonio", "重氮鎓基", anchored=("*[N+]#N",), paren=False),  # P-65.3：重氮鎓 Ar-N2+，N+ 直连母体作 diazonio 前缀
    "cyano": RetainedSubstituent("cyano", "氰基", anchored=("*C#N",), paren=False),
    "isocyano": RetainedSubstituent("isocyano", "异氰基", anchored=("*[N+]#[C-]",), paren=False),
}


def _canon(smi: str) -> str:
    """锚定 SMILES 的 canonical 形式；解析失败抛错。"""
    from rdkit.Chem import MolFromSmiles, MolToSmiles

    m = MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"[anchored_table] 无法解析锚定键: {smi!r}")
    return MolToSmiles(m)


def _build_anchor_index() -> dict[str, str]:
    """反查锚定键→registry key，校验 canonical 且唯一。"""
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


def _resolve_name(key: str) -> tuple[str, str]:
    """返回 registry 键对应的 (en, zh)。"""
    entry = _REGISTRY[key]
    return entry.en, entry.zh

def anchored_key(mol: Mol, atoms: frozenset[int], attach_old: int | None = None) -> str | None:
    """生成原子集的锚定 canonical-SMILES 键，失败时返回 None。"""
    return memo.by_key("anchored_key", (id(mol), atoms, attach_old),
                       lambda: _anchored_key_uncached(mol, atoms, attach_old), mol)


def _anchored_key_uncached(mol: Mol, atoms: frozenset[int], attach_old: int) -> str | None:
    """构建锚定子分子并取 canonical SMILES（无记忆版）。"""
    from rdkit.Chem import MolToSmiles

    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None

_WHOLE_ONLY_KEYS = frozenset({"*O", "*[O]", "*N"})  # 整分子顶层才命中的锚定键：游离 O/N 自由基（表 2.1 去氢），L3 跳过


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> str | None:
    """查 registry 锚定索引返回 key；无命中返回 None。"""
    key = anchored_key(mol, atoms, attach_old)
    if key is None or key in _WHOLE_ONLY_KEYS:
        return None
    return _ANCHOR_INDEX.get(key)


def _attached_via_nitrogen(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> bool:
    """块连接原子的块外邻居是否为氮（经氮连母体，即胍的亚氨基桥）。"""
    if attach_old is None:
        return False
    return any(nb.GetAtomicNum() == N and nb.GetIdx() not in atoms
               for nb in mol.GetAtomWithIdx(attach_old).GetNeighbors())


def _ester_o_side(mol: Mol, atoms: frozenset[int], attach_old: int) -> tuple[frozenset[int], int] | None:
    """块恰为 -C(=O)-O-R（连接点=羰基碳）时返回 (O 侧原子集, 单键氧)，否则 None。"""
    at = mol.GetAtomWithIdx(attach_old)
    if at.GetAtomicNum() != C:
        return None
    dbl: list[int] = []
    sgl: list[int] = []
    for nb in at.GetNeighbors():
        if nb.GetIdx() not in atoms:  # 块外邻居（母体侧）不参与判定
            continue
        if nb.GetAtomicNum() != O:
            return None
        bond = mol.GetBondBetweenAtoms(attach_old, nb.GetIdx())
        (dbl if bond.GetBondType() == BondType.DOUBLE else sgl).append(nb.GetIdx())
    if len(dbl) != 1 or len(sgl) != 1:
        return None
    front = frozenset(a for a in atoms if a not in (attach_old, dbl[0]))
    return (front, sgl[0]) if len(front) > 1 else None


def _alkoxycarbonyl(mol: Mol, atoms: frozenset[int], attach_old: int) -> tuple[str, str, bool] | None:
    """-C(=O)-O-R 收成 <R>oxycarbonyl / <R>氧羰基（P-65.1.7.2 酯作前缀）。"""
    from namepredict.layer3.as_substituent import name_as_substituent

    side = _ester_o_side(mol, atoms, attach_old)
    if side is None:
        return None
    front, o_idx = side
    hit = name_as_substituent(mol, o_idx, front)
    if hit is None or not hit[0].endswith("oxy"):
        return None
    return hit[0] + "carbonyl", zh_bridge_root(hit[1]) + "羰基", False


_RING_YL_EN_RE = re.compile(r"^(.+?)-(\d+[a-z]?)-yl$")  # 环胺 N-侧基名（pyrrolidin-1-yl）
_FREE_VALENCE_YL_RE = re.compile(r"-\d+[a-z]?-yl$")  # 带位次的自由价词尾（…pentan-2-yl / …triazin-2-yl）
_AMINO_EN, _AMINO_ZH = "amino", "氨基"
_ANILINO_EN, _ANILINO_ZH = "anilino", "苯胺基"
_PHENYL_EN, _PHENYL_ZH = "phenyl", "苯基"



def carbamoyl_prefix_name(en: str, zh: str, *, in_ring: bool) -> tuple[str, str] | None:
    """N-侧胺名 → carbamoyl 一族前缀名（P-65.2.1.5）；不合式返回 None。"""
    if not en or not zh:
        return None
    if in_ring:  # N 为环员：按酰基取环母体名 <母体>-<位次>-carbonyl
        m = _RING_YL_EN_RE.match(en)
        if m is None or not zh.endswith("基"):
            return None
        stem, loc = m.group(1), m.group(2)
        parent = stem if stem.endswith("e") else stem + "e"  # -yl 词干补回母体氢化物词尾
        return f"{parent}-{loc}-carbonyl", zh[:-1] + "羰基"
    if en.endswith(_AMINO_EN):  # N-取代胺名去「氨基」接 carbamoyl
        stem_en, stem_zh = en[: -len(_AMINO_EN)], zh[: -len(_AMINO_ZH)]
        fenced = stem_en.startswith("(") and stem_en.endswith(")")
        if "(" in stem_en and not fenced:
            if stem_en.endswith(")"):  # 取代胺名（methyl(propyl)）：围栏须连 carbamoyl 一并括起（P-16.5.1.1）
                return f"[{stem_en}carbamoyl]", f"[{stem_zh}氨基甲酰基]"
            if _FREE_VALENCE_YL_RE.search(stem_en):  # 自由价词尾（…pentan-2-yl）：围栏止于词尾（P-16.5.1.2）
                return f"[{stem_en}]carbamoyl", f"[{stem_zh}]氨基甲酰基"
            return stem_en + "carbamoyl", stem_zh + "氨基甲酰基"  # 前导位次隔开的取代基括号无需围栏
        if len(SUB_LOCANT_RE.findall(stem_en)) >= 2 and not fenced:  # 复合环名前导多位次段：括起消歧（P-16.5.1.3.1）
            return f"({stem_en})carbamoyl", f"({stem_zh})氨基甲酰基"
        return stem_en + "carbamoyl", stem_zh + "氨基甲酰基"
    if en.endswith(_ANILINO_EN):  # N-芳基用 anilino 保留名，须换回「芳基」再接 carbamoyl
        aryl_en = en[: -len(_ANILINO_EN)] + _PHENYL_EN
        aryl_zh = zh[: -len(_ANILINO_ZH)] + _PHENYL_ZH
        if aryl_en[:1].isdigit():  # 前导位次须括起消歧（P-16.5.1.3.1）
            aryl_en, aryl_zh = f"({aryl_en})", f"({aryl_zh})"
        return aryl_en + "carbamoyl", aryl_zh + "氨基甲酰基"
    return None


def _fence_arm(name: str) -> str:
    """复合臂名围栏：名内已含任何括号改方括号，否则圆括号（P-16.5.1.3.1）。"""
    return f"[{name}]" if any(ch in name for ch in "()[]") else f"({name})"


def _amidine_arm(mol: Mol, root: int, atoms: frozenset[int], hub: frozenset[int]) -> tuple[str, str, bool] | None:
    """把母体氮外侧键起的臂命名（返回 en, zh, 是否须围栏）；失败返回 None。"""
    from namepredict.layer3.as_substituent import name_as_substituent
    from namepredict.tools.block_cut import cut_block

    arm = cut_block(mol, root, hub)
    arm = frozenset(arm) & atoms if arm else arm
    if not arm:
        return None
    return name_as_substituent(mol, root, arm)


def _amidinium_cation_arm(mol: Mol, atoms: frozenset[int], attach_old: int) -> tuple[str, str] | None:
    """root 为脒碳、块外双键连阳离子氮（azanium 母体）→ diaminomethylidene(<N 取代基>) 名。"""
    at = mol.GetAtomWithIdx(attach_old)
    if at.GetAtomicNum() != C:
        return None
    dbl_cation = None
    sgl_ns: list[int] = []
    for n in at.GetNeighbors():
        bt = mol.GetBondBetweenAtoms(attach_old, n.GetIdx()).GetBondType()
        if n.GetIdx() not in atoms:
            if n.GetAtomicNum() == N and n.GetFormalCharge() == 1 and bt == BondType.DOUBLE:
                dbl_cation = n.GetIdx()  # 阳离子亚胺氮即 azanium 母体，位于块外
            continue
        if n.GetAtomicNum() == N and bt == BondType.SINGLE:
            sgl_ns.append(n.GetIdx())
    if dbl_cation is None or not sgl_ns:
        return None
    if any(n.GetAtomicNum() != 1 and n.GetIdx() not in (dbl_cation, *sgl_ns)
           for n in at.GetNeighbors()):
        return None  # 脒碳还挂别的重原子：非本式
    hub = frozenset((attach_old, dbl_cation, *sgl_ns))
    roots = [n.GetIdx() for a in sgl_ns for n in mol.GetAtomWithIdx(a).GetNeighbors()
             if n.GetAtomicNum() != 1 and n.GetIdx() in atoms and n.GetIdx() not in hub]
    if not roots:
        return "diaminomethylidene", "二氨基亚甲基"  # P-34：亚胺氮留在阳离子母体，块内为二氨基亚甲基
    arms = [_amidine_arm(mol, r, atoms, hub) for r in roots]
    if any(a is None for a in arms):
        return None
    en = ", ".join(sorted(a[0] for a in arms))
    zh = "、".join(sorted(a[1] for a in arms))
    return f"diaminomethylidene({en})", f"二氨基亚甲基({zh})"


def _amidine_group(na: list[tuple[str, str, bool]], ni: list[tuple[str, str, bool]],
                   suffix: str, conj: str, idx: int) -> str:
    """氨基/亚氨基侧取代基名 → N/N' 位次前缀段 + 保留名（P-66.4.1.3.1）；idx 选 en/zh。"""
    if not na and not ni:
        return suffix
    na_s = [a[idx] for a in na]
    ni_s = [a[idx] for a in ni]
    if len(na_s) == 1 and len(ni_s) == 1 and na_s[0] == ni_s[0] and na_s[0].isalpha():
        return f"N,N'-{conj}{na_s[0]}{suffix}"  # 同名简单基：合并计数（N,N'-二甲基）
    segs: list[str] = []
    for tokens, arms in ((("N",) * len(na), na), (("N'",) * len(ni), ni)):
        if arms:
            segs.append(f"{','.join(tokens)}-{'-'.join(_fence_arm(a[idx]) if a[2] else a[idx] for a in arms)}")
    return "-".join(segs) + suffix


def _n_substituted_amidine(mol: Mol, atoms: frozenset[int], attach_old: int) -> tuple[str, str] | None:
    """块为 N-取代脒（root 氮连 C(=N)N 且任一脒氮带取代基）时的 carbamimidoyl 式前缀名。"""
    if attach_old is None or attach_old not in atoms:
        return None
    nb_atom = mol.GetAtomWithIdx(attach_old)
    if nb_atom.GetAtomicNum() != N or nb_atom.GetFormalCharge() != 0 \
            or nb_atom.GetIsAromatic() or nb_atom.IsInRing():
        return None
    c_idx = next((n.GetIdx() for n in nb_atom.GetNeighbors()
                  if n.GetAtomicNum() == C and n.GetIdx() in atoms), None)
    if c_idx is None:
        return None
    dbl_n = sgl_n = None
    for n in mol.GetAtomWithIdx(c_idx).GetNeighbors():
        if n.GetAtomicNum() != N or n.GetIdx() not in atoms or n.GetIdx() == attach_old:
            continue  # attach_old 即桥氮，另两个氮才是脒的亚氨基/氨基
        bt = mol.GetBondBetweenAtoms(c_idx, n.GetIdx()).GetBondType()
        if bt == BondType.DOUBLE:
            dbl_n = n.GetIdx()
        elif bt == BondType.SINGLE:
            sgl_n = n.GetIdx()
    if dbl_n is None or sgl_n is None:
        return None
    if any(n.GetAtomicNum() != 1 and n.GetIdx() not in (attach_old, dbl_n, sgl_n)
           for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()):
        return None  # 脒碳还挂别的重原子：非本式
    hub = frozenset((attach_old, c_idx, dbl_n, sgl_n))
    nb_arms = [n.GetIdx() for n in nb_atom.GetNeighbors()
               if n.GetAtomicNum() != 1 and n.GetIdx() in atoms and n.GetIdx() not in hub]
    na_arms = [n.GetIdx() for n in mol.GetAtomWithIdx(sgl_n).GetNeighbors()
               if n.GetAtomicNum() != 1 and n.GetIdx() in atoms and n.GetIdx() not in hub]
    ni_arms = [n.GetIdx() for n in mol.GetAtomWithIdx(dbl_n).GetNeighbors()
               if n.GetAtomicNum() != 1 and n.GetIdx() in atoms and n.GetIdx() not in hub]
    if not (nb_arms or na_arms or ni_arms):
        return None  # 无取代：留给表内 diaminomethylidene / 递归路径
    en_arms = [_amidine_arm(mol, a, atoms, hub) for a in (*nb_arms, *na_arms, *ni_arms)]
    if any(n is None for n in en_arms):
        return None
    nb_en, na_en, ni_en = en_arms[:len(nb_arms)], \
        en_arms[len(nb_arms):len(nb_arms) + len(na_arms)], en_arms[len(nb_arms) + len(na_arms):]
    group_en = _amidine_group(na_en, ni_en, "carbamimidoyl", "di", 0)
    group_zh = _amidine_group(na_en, ni_en, "氨基甲亚氨酰基", "二", 1)
    nb_names = [_fence_arm(a[0]) if a[2] else a[0] for a in nb_en]
    nb_zhens = [_fence_arm(a[1]) if a[2] else a[1] for a in nb_en]
    en_parts = sorted([group_en, *nb_names], key=lambda s: s.lstrip("N,%-'"))
    zh_parts = sorted([group_zh, *nb_zhens], key=lambda s: s.lstrip("N,%-'"))
    if len(en_parts) == 1:
        return _fence_arm(en_parts[0]) + "amino", _fence_arm(zh_parts[0]) + "氨基"
    return (en_parts[0] + _fence_arm(en_parts[1]) + "amino",
            zh_parts[0] + _fence_arm(zh_parts[1]) + "氨基")


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
) -> tuple[str, str, bool] | None:
    """查找取代基原子集，返回 (en, zh, paren)；无命中返回 None。"""
    reg_key = _table_hit(mol, atoms, attach_old)
    if reg_key is not None:
        # P-34/P-66.4.1.3.1：-C(=NH)NH2 不经 N 连母体（C/S/O/P…）一律取 carbamimidoyl；经 N 连（胍亚氨基桥）仍用 diaminomethylidene
        if reg_key == "diaminomethylidene" and not _attached_via_nitrogen(mol, atoms, attach_old):
            return "carbamimidoyl", "氨基甲亚氨酰基", False
        en, zh = _resolve_name(reg_key)
        return en, zh, _REGISTRY[reg_key].paren
    if attach_old is None:
        return None
    amidine = _n_substituted_amidine(mol, atoms, attach_old)
    if amidine is not None:
        return amidine[0], amidine[1], True
    amidinium = _amidinium_cation_arm(mol, atoms, attach_old)
    if amidinium is not None:
        return amidinium[0], amidinium[1], True
    return _alkoxycarbonyl(mol, atoms, attach_old)
