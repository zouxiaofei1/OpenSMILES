"""取代基锚定 canonical-SMILES 查表，未命中回退完整命名路径；
简单/保留取代基统一存 registry，构建期校验锚定键唯一。"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem import BondType, Mol

from namepredict.constants import C, O, zh_bridge_root
from namepredict.tools import memo
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
    "sulfanyl": RetainedSubstituent("sulfanyl", "巯基", anchored=("*S",), paren=False),
    "sulfido": RetainedSubstituent("sulfido", "硫代", anchored=("*[S-]",), paren=False),  # 去质子硫醇负离子（P-66.1.1.4 阴离子前缀）
    "selanyl": RetainedSubstituent("selanyl", "氢硒基", anchored=("*[SeH]",), paren=False),
    "methylsulfinyl": RetainedSubstituent("methylsulfinyl", "甲基亚磺酰", anchored=("*S(C)=O",), paren=False),
    "methylsulfonyl": RetainedSubstituent("methylsulfonyl", "甲磺酰基", anchored=("*S(C)(=O)=O",), paren=False),
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
    "formyl": RetainedSubstituent("formyl", "甲酰", anchored=("*C=O",), paren=False),
    "carboxymethyl": RetainedSubstituent("carboxymethyl", "羧甲基", anchored=("*CC(=O)O",), paren=True),
    "hydroxymethyl": RetainedSubstituent("hydroxymethyl", "羟甲基", anchored=("*CO",), paren=True),
    "nitroso": RetainedSubstituent("nitroso", "亚硝基", anchored=("*N=O",), paren=False),
    "azido": RetainedSubstituent("azido", "叠氮基", anchored=("*N=[N+]=[N-]",), paren=False),
    "amino": RetainedSubstituent("amino", "氨基", anchored=("*N",), paren=False),
    "hydrazinyl": RetainedSubstituent("hydrazinyl", "肼基", anchored=("*NN",), paren=False),
    "anilino": RetainedSubstituent("anilino", "苯胺基", anchored=("*Nc1ccccc1",), paren=False),
    "diazenyl": RetainedSubstituent("diazenyl", "二氮烯基", anchored=("*N=N",), paren=False),
    "diazo": RetainedSubstituent("diazo", "重氮基", anchored=("*[N+]=[N-]",), paren=False),
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


def resolve_name(key: str) -> tuple[str, str]:
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


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
) -> tuple[str, str, bool] | None:
    """查找取代基原子集，返回 (en, zh, paren)；无命中返回 None。"""
    reg_key = _table_hit(mol, atoms, attach_old)
    if reg_key is not None:
        en, zh = resolve_name(reg_key)
        return en, zh, _REGISTRY[reg_key].paren
    if attach_old is None:
        return None
    return _alkoxycarbonyl(mol, atoms, attach_old)
