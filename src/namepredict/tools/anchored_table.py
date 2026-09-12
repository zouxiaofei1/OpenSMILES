"""取代基锚定 canonical-SMILES 查表（registry 经 anchored 反查索引），未命中回退完整命名路径；
简单/保留取代基统一存 registry；构建期校验 canonical 对拍与锚定键唯一性。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

from namepredict.tools import memo
from namepredict.layer3.submol_build import build_anchor_submol


# ── 取代基注册表（IUPAC 2013 P-29/P-57；general/pin 双语名称）──

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
    """构建取代基注册表（键 → RetainedSubstituent 条目）；简单取代基一并入表。"""
    P, G, N = IupacLevel.PIN, IupacLevel.GENERAL, IupacLevel.NOT_RECOMMENDED
    out: dict[str, RetainedSubstituent] = {  # 基础取代基（原 ANCHOR_TABLE 并入 registry）
        "fluoro": RetainedSubstituent( "fluoro", "氟", "fluoro", "氟", P, anchored=("*F", ), paren=False, kind="halo", ),
        "chloro": RetainedSubstituent( "chloro", "氯", "chloro", "氯", P, anchored=("*Cl", ), paren=False, kind="halo", ),
        "bromo": RetainedSubstituent( "bromo", "溴", "bromo", "溴", P, anchored=("*Br", ), paren=False, kind="halo", ),
        "iodo": RetainedSubstituent( "iodo", "碘", "iodo", "碘", P, anchored=("*I", ), paren=False, kind="halo", ),
        "nitro": RetainedSubstituent( "nitro", "硝基", "nitro", "硝基", P, anchored=("*[N+](=O)[O-]", ), paren=False, kind="leaf", ),
        "oxo": RetainedSubstituent( "oxo", "氧代", "oxo", "氧代", P, anchored=("*=O", ), paren=False, kind="leaf", ),
        "nitro": RetainedSubstituent( "nitro", "硝基", "nitro", "硝基", P, anchored=("*[N+](=O)[O-]", ), paren=False, kind="leaf", ),
        "isocyanato": RetainedSubstituent( "isocyanato", "异氰酸基", "isocyanato", "异氰酸基", P, anchored=("*N=C=O", ), paren=False, kind="leaf", ),
        "isothiocyanato": RetainedSubstituent( "isothiocyanato", "异硫氰酸基", "isothiocyanato", "异硫氰酸基", P, anchored=("*N=C=S", ), paren=False, kind="leaf", ),
        "methyl": RetainedSubstituent( "methyl", "甲基", "methyl", "甲基", P, anchored=("*C", ), paren=False, kind="alkyl", ),
        "methylidene": RetainedSubstituent( "methylidene", "亚甲基", "methylidene", "亚甲基", P, anchored=("*=C", ), paren=False, kind="alkyl", ),
        "ethylidene": RetainedSubstituent( "ethylidene", "亚乙基", "ethylidene", "亚乙基", P, anchored=("*=CC", ), paren=False, kind="alkyl", ),
        "propylidene": RetainedSubstituent( "propylidene", "亚丙基", "propylidene", "亚丙基", P, anchored=("*=CCC", ), paren=False, kind="alkyl", ),
        "cyclopropylidene": RetainedSubstituent( "cyclopropylidene", "环丙亚基", "cyclopropylidene", "环丙亚基", P, anchored=("*=C1CC1", ), paren=False, kind="alkyl", ),
        "cyclohexylidene": RetainedSubstituent( "cyclohexylidene", "环己亚基", "cyclohexylidene", "环己亚基", P, anchored=("*=C1CCCCC1", ), paren=False, kind="alkyl", ),
        "sulfanylidene": RetainedSubstituent( "sulfanylidene", "硫烷亚基", "sulfanylidene", "硫烷亚基", P, anchored=("*=S", ), paren=False, kind="leaf", ),
        "diaminomethylidene": RetainedSubstituent( "diaminomethylidene", "二氨基亚甲基", "diaminomethylidene", "二氨基亚甲基", P, anchored=("*C(=N)N", ), paren=True, kind="leaf", ),  # 脒/胍残基 C(=N)N：gold 按 P-66.1.1 取亚基式（diaminomethylideneamino），不取等价的 amino(imino)methylamino（两者互变、分子式相同，取测试集口径）。
        "ethyl": RetainedSubstituent( "ethyl", "乙基", "ethyl", "乙基", P, anchored=("*CC", ), paren=False, kind="alkyl", ),
        "propyl": RetainedSubstituent( "propyl", "丙基", "propyl", "丙基", P, anchored=("*CCC", ), paren=False, kind="alkyl", ),
        "butyl": RetainedSubstituent( "butyl", "丁基", "butyl", "丁基", P, anchored=("*CCCC", ), paren=False, kind="alkyl", ),
        "tert-butyl": RetainedSubstituent( "tert-butyl", "叔丁基", "tert-butyl", "叔丁基", P, anchored=("*C(C)(C)C", ), paren=False, kind="alkyl", ),  # 支链 / 不饱和烷基
        "isopropyl": RetainedSubstituent( "propan-2-yl", "丙-2-基", "propan-2-yl", "丙-2-基", G, anchored=("*C(C)C", ), paren=False, kind="alkyl", ),
        "isobutyl": RetainedSubstituent("isobutyl", "异丁基", "2-methylpropyl", "2-甲基丙基", N, anchored=("*CC(C)C", ), paren=False, kind="alkyl", ),
        "sec-butyl": RetainedSubstituent("sec-butyl", "仲丁基", "butan-2-yl", "丁-2-基", N, anchored=("*C(C)CC", ), paren=False, kind="alkyl", ),
        "neopentyl": RetainedSubstituent( "neopentyl", "新戊基", "2,2-dimethylpropyl", "2,2-二甲基丙基", N, anchored=("*CC(C)(C)C", ), paren=False, kind="alkyl", ),
        "isopentyl": RetainedSubstituent( "isopentyl", "异戊基", "3-methylbutyl", "3-甲基丁基", N, anchored=("*CCC(C)C", ), paren=False, kind="alkyl", ),
        "vinyl": RetainedSubstituent( "vinyl", "乙烯基", "ethenyl", "乙烯基", G, anchored=("*C=C", ), paren=False, kind="alkyl", ),
         "allyl": RetainedSubstituent( "allyl", "烯丙基", "prop-2-enyl", "丙-2-烯基", G, anchored=("*CC=C", ), paren=False, kind="alkyl", ),
        "isopropenyl": RetainedSubstituent( "prop-1-en-2-yl", "异丙烯基", "prop-1-en-2-yl", "丙-1-烯-2-基", G, anchored=("*C(=C)C", ), paren=False, kind="alkyl", ),
               "propargyl": RetainedSubstituent( "propargyl", "炔丙基", "prop-2-ynyl", "丙-2-炔基", G, anchored=("*CC#C", ), paren=False, kind="alkyl", ),
        "benzyl": RetainedSubstituent( "benzyl", "苄基", "benzyl", "苄基", P, anchored=("*Cc1ccccc1", ), paren=False, kind="aryl", ),
        "methoxy": RetainedSubstituent( "methoxy", "甲氧基", "methoxy", "甲氧基", P, anchored=("*OC", ), paren=False, kind="leaf", ),
        "hydroperoxy": RetainedSubstituent( "hydroperoxy", "氢过氧基", "hydroperoxy", "氢过氧基", P, anchored=("*OO", ), paren=False, kind="leaf", ),
        "hydroxy": RetainedSubstituent( "hydroxy", "羟基", "hydroxy", "羟基", P, anchored=("*O", ), paren=False, kind="leaf", ),
        "oxidanyl": RetainedSubstituent( "oxidanyl", "氧基", "oxidanyl", "氧基", P, anchored=("*[O]", ), paren=False, kind="leaf", ),
        "methylsulfanyl": RetainedSubstituent( "methylsulfanyl", "甲硫基", "methylsulfanyl", "甲硫基", P, anchored=("*SC", ), paren=False, kind="leaf", ),
        "ethylsulfanyl": RetainedSubstituent( "ethylsulfanyl", "乙硫基", "ethylsulfanyl", "乙硫基", P, anchored=("*SCC", ), paren=False, kind="leaf", ),
        "sulfanyl": RetainedSubstituent( "sulfanyl", "巯基", "sulfanyl", "巯基", P, anchored=("*S", ), paren=False, kind="leaf", ),
        "selanyl": RetainedSubstituent( "selanyl", "氢硒基", "selanyl", "氢硒基", P, anchored=("*[SeH]", ), paren=False, kind="leaf", ),
        "methylsulfinyl": RetainedSubstituent( "methylsulfinyl", "甲基亚磺酰", "methylsulfinyl", "甲基亚磺酰", P, anchored=("*S(C)=O", ), paren=False, kind="leaf", ),
        "methylsulfonyl": RetainedSubstituent( "methylsulfonyl", "甲磺酰基", "methylsulfonyl", "甲磺酰基", P, anchored=("*S(C)(=O)=O", ), paren=False, kind="leaf", ),
        "sulfo": RetainedSubstituent( "sulfo", "磺基", "sulfo", "磺基", P, anchored=("*S(=O)(=O)O", ), paren=False, kind="leaf", ),
        "tosyl": RetainedSubstituent( "tosyl", "对甲苯磺酰基", "4-methylbenzenesulfonyl", "4-甲基苯磺酰基", N, anchored=("*S(=O)(=O)c1ccc(C)cc1", ), paren=False, kind="leaf", ),
        "carboxy": RetainedSubstituent( "carboxy", "羧基", "carboxy", "羧基", P, anchored=("*C(=O)O", ), paren=False, kind="leaf", ),
        "carbamoyl": RetainedSubstituent( "carbamoyl", "氨基甲酰基", "carbamoyl", "氨基甲酰基", P, anchored=("*C(N)=O", ), paren=False, kind="leaf", ),  # P-66.1.1.4.1 氨基甲酸（carbamic acid）的酰基保留前缀；gold/ChEBI 全量 54 处取 carbamoyl，不取 aminocarbonyl/amino(oxo)methyl
        "carbamoylamino": RetainedSubstituent( "carbamoylamino", "氨基甲酰氨基", "carbamoylamino", "氨基甲酰氨基", P, anchored=("*NC(N)=O", ), paren=True, kind="leaf", ),  # P-66.1.1.6 ureido 在 IUPAC 已不推荐（P_1 附录：ureido/ureylene 不用），优选 carbamoylamino
        # 铵/𬭩型阳离子取代基（P-62.4.1：铵 azanium 去氢得 azaniumyl 型前缀）。gold 全量 174 例含 azanium*，        # 现行管线把这些带电 N 片段整体丢弃（no_coverage_gate / coverage_complete 误判），故按锚定叶子入表。
        "azaniumyl": RetainedSubstituent( "azaniumyl", "铵基", "azaniumyl", "铵基", P, anchored=("*[NH3+]", ), paren=False, kind="leaf", ),
        "methylazaniumyl": RetainedSubstituent( "methylazaniumyl", "甲基铵基", "methylazaniumyl", "甲基铵基", P, anchored=("*[NH2+]C", ), paren=False, kind="leaf", ),
        "dimethylazaniumyl": RetainedSubstituent( "dimethylazaniumyl", "二甲基铵基", "dimethylazaniumyl", "二甲基铵基", P, anchored=("*[NH+](C)C", ), paren=True, kind="leaf", ),
        "trimethylazaniumyl": RetainedSubstituent( "trimethylazaniumyl", "三甲基铵基", "trimethylazaniumyl", "三甲基铵基", P, anchored=("*[N+](C)(C)C", ), paren=True, kind="leaf", ),
        "carbamoyloxy": RetainedSubstituent( "carbamoyloxy", "氨基甲酰氧基", "carbamoyloxy", "氨基甲酰氧基", P, anchored=("*OC(N)=O", ), paren=False, kind="leaf", ),  # 氨基甲酸 O-酯残基（P-66.1.1.4.1）
        "carbamothioylamino": RetainedSubstituent( "carbamothioylamino", "氨基硫代羰基氨基", "carbamothioylamino", "氨基硫代羰基氨基", P, anchored=("*NC(N)=S", ), paren=True, kind="leaf", ),  # 硫代氨基甲酸残基（P-66.1.1.4：carbamothioyl）
        "sulfamoyl": RetainedSubstituent( "sulfamoyl", "氨磺酰基", "sulfamoyl", "氨磺酰基", P, anchored=("*S(N)(=O)=O", ), paren=False, kind="leaf", ),  # P-66.1.1.4.2 磺酰胺（sulfamoyl = H2N-SO2-）；N-取代时基名随取代基前移
        "phosphono": RetainedSubstituent( "phosphono", "膦酸", "phosphono", "膦酸", P, anchored=("*P(=O)(O)O", ), paren=False, kind="leaf", ),  # P-102：phosphono 表示 -PO(OH)2，P 直连母体（对比 O 桥的 phosphonooxy）
        "phosphonato": RetainedSubstituent( "phosphonato", "膦酸根", "phosphonato", "膦酸根", P, anchored=("*P(=O)([O-])O", "*P(=O)([O-])[O-]"), paren=False, kind="leaf", ),  # P-102：phosphonato 表示 -PO(O-)2（单/双阴离子）
        "phosphonooxy": RetainedSubstituent( "phosphonooxy", "膦酸氧基", "phosphonooxy", "膦酸氧基", P, anchored=("*OP(=O)(O)O", ), paren=False, kind="leaf", ),  # 磷酸降级前缀（P-67.1.5.1：羧酸等更高优先级 FG 存在时磷酸以 phosphonooxy 前缀表达）
        "phosphonatooxy": RetainedSubstituent( "phosphonatooxy", "膦酸氧基", "phosphonatooxy", "膦酸氧基", P, anchored=("*OP(=O)([O-])O", "*OP(=O)([O-])[O-]"), paren=False, kind="leaf", ),
        "phosphonooxymethyl": RetainedSubstituent( "phosphonooxymethyl", "膦酸氧甲基", "phosphonooxymethyl", "膦酸氧甲基", P, anchored=("*COP(=O)(O)O", ), paren=True, kind="leaf", ),
        "phosphonatooxymethyl": RetainedSubstituent( "phosphonatooxymethyl", "膦酸氧甲基", "phosphonatooxymethyl", "膦酸氧甲基", P, anchored=("*COP(=O)([O-])O", ), paren=True, kind="leaf", ),
        "formyl": RetainedSubstituent( "formyl", "甲酰", "formyl", "甲酰", P, anchored=("*C=O", ), paren=False, kind="leaf", ),
        "carboxymethyl": RetainedSubstituent( "carboxymethyl", "羧甲基", "carboxymethyl", "羧甲基", P, anchored=("*CC(=O)O", ), paren=True, kind="leaf", ),
        "hydroxymethyl": RetainedSubstituent( "hydroxymethyl", "羟甲基", "hydroxymethyl", "羟甲基", P, anchored=("*CO", ), paren=True, kind="leaf", ),
        "nitroso": RetainedSubstituent( "nitroso", "亚硝基", "nitroso", "亚硝基", P, anchored=("*N=O", ), paren=False, kind="leaf", ),
        "azido": RetainedSubstituent( "azido", "叠氮基", "azido", "叠氮基", P, anchored=("*N=[N+]=[N-]", ), paren=False, kind="leaf", ),
        "amino": RetainedSubstituent( "amino", "氨基", "amino", "氨基", P, anchored=("*N", ), paren=False, kind="leaf", ),
        "hydrazinyl": RetainedSubstituent( "hydrazinyl", "肼基", "hydrazinyl", "肼基", P, anchored=("*NN", ), paren=False, kind="leaf", ),
        "anilino": RetainedSubstituent( "anilino", "苯胺基", "phenylamino", "苯氨基", P, anchored=("*Nc1ccccc1", ), paren=False, kind="leaf", ),
        "diazenyl": RetainedSubstituent( "diazenyl", "二氮烯基", "diazenyl", "二氮烯基", P, anchored=("*N=N", ), paren=False, kind="leaf", ),
        "diazo": RetainedSubstituent( "diazo", "重氮基", "diazo", "重氮基", P, anchored=("*[N+]=[N-]", ), paren=False, kind="leaf", ),
        "cyano": RetainedSubstituent( "cyano", "氰基", "cyano", "氰基", P, anchored=("*C#N", ), paren=False, kind="leaf", ),
        "isocyano": RetainedSubstituent( "isocyano", "异氰基", "isocyano", "异氰基", P, anchored=("*[N+]#[C-]", ), paren=False, kind="leaf", ),
    }
    return out


_REGISTRY: dict[str, RetainedSubstituent] = _build_registry()


def _canon(smi: str) -> str:
    """锚定 SMILES 的 canonical 形式；解析失败抛错。"""
    from rdkit.Chem import MolFromSmiles, MolToSmiles

    m = MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"[anchored_table] 无法解析锚定键: {smi!r}")
    return MolToSmiles(m)


def _build_anchor_index() -> dict[str, str]:
    """从 registry 反查锚定 canonical-SMILES → registry key，并校验 anchored 键已是 canonical 形式且跨 key 唯一。"""
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


def resolve_name(key: str, *, name_mode: str = "general") -> tuple[str, str]:
    """返回 registry 键对应的 (en, zh)：仅 PIN 级条目取保留名，其余（含 general 级）取系统名。"""
    entry = _REGISTRY[key]
    if name_mode == "pin" and entry.level != IupacLevel.PIN:
        return entry.systematic_en, entry.systematic_zh  # PIN 模式：非 PIN 级保留名回落系统名
    if entry.level == IupacLevel.PIN:
        return entry.en, entry.zh  # 唯一有独立保留名的 PIN 级条目为 anilino（苯胺基）；gold 全量 41:0 取保留名
    return entry.systematic_en, entry.systematic_zh  # general/not_rec 级：gold 对 vinyl/isobutyl/tosyl 等一律取系统名


def pick_root(mol: Mol, atoms: frozenset[int]) -> int:
    """取代基侧键合原子：atoms 中与外部重原子（母体）成键者，无则回退最小索引。"""
    for a in atoms:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in atoms:
                return a
    return min(atoms)


def anchored_key(mol: Mol, atoms: frozenset[int], attach_old: int | None = None) -> str | None:
    """生成原子集的锚定 canonical-SMILES 键，失败时返回 None。

    同一个 claim 会在查表与递归拆分之间被反复求键（实测约一半重复），而键只依赖
    (mol, 原子集, 连接原子)，故按此记忆；仍限定在单次命名内，不跨分子共享。
    """
    root = pick_root(mol, atoms) if attach_old is None else attach_old
    return memo.by_key("anchored_key", (id(mol), atoms, root),
                       lambda: _anchored_key_uncached(mol, atoms, root), mol)


def _anchored_key_uncached(mol: Mol, atoms: frozenset[int], attach_old: int) -> str | None:
    """实际构建锚定子分子并取 canonical SMILES（无记忆版本，见 anchored_key）。"""
    from rdkit.Chem import MolToSmiles

    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None

_WHOLE_ONLY_KEYS = frozenset({"*O", "*[O]", "*N"})  # 仅整分子顶层命中的锚定键：单原子杂原子自由基（表 2.1 去氢）。L3 取代基查表（anchored_entry）跳过它们——游离 O/N 原子会被 _one_alkyl 误作侧链提取（酮羰基氧、酯氧、胺氮），命中 *O/*N 会错名成羟基/氨基。


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> str | None:
    """查 registry 锚定索引，返回 registry key；无命中返回 None。"""
    key = anchored_key(mol, atoms, attach_old)
    if key is None or key in _WHOLE_ONLY_KEYS:
        return None
    return _ANCHOR_INDEX.get(key)


def anchored_entry(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None, *, name_mode: str = "general",
) -> tuple[str, str, bool, str] | None:
    """在 name_mode 下将原子集解析为 (en, zh, paren, kind)，无命中则返回 None。"""
    reg_key = _table_hit(mol, atoms, attach_old)
    if reg_key is None:
        return None
    en, zh = resolve_name(reg_key, name_mode=name_mode)
    e = _REGISTRY[reg_key]
    return en, zh, e.paren, e.kind


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
    *, name_mode: str = "general",
) -> tuple[str, str, bool] | None:
    """查找取代基原子集，按 name_mode 经 resolve_name 返回 (en, zh, paren)；无命中返回 None。"""
    entry = anchored_entry(mol, atoms, attach_old, name_mode=name_mode)
    return None if entry is None else (entry[0], entry[1], entry[2])


def anchored_whole_mol(mol: Mol, *, name_mode: str = "general") -> tuple[str, str, bool, str] | None:
    """整分子 canonical SMILES 命中锚定表时返回 (en, zh, paren, kind)。"""
    from rdkit.Chem import MolToSmiles

    key = MolToSmiles(mol)
    reg_key = _ANCHOR_INDEX.get(key)
    if reg_key is None:
        return None
    en, zh = resolve_name(reg_key, name_mode=name_mode)
    e = _REGISTRY[reg_key]
    return en, zh, e.paren, e.kind
