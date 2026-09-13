"""取代基锚定 canonical-SMILES 查表（registry 经 anchored 反查索引），未命中回退完整命名路径；
简单/保留取代基统一存 registry；构建期校验 canonical 对拍与锚定键唯一性。"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem import Mol

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


_REGISTRY: dict[str, RetainedSubstituent] = {  # 基础取代基（原 ANCHOR_TABLE 并入 registry）
    "fluoro": RetainedSubstituent("fluoro", "氟", anchored=("*F",), paren=False),
    "chloro": RetainedSubstituent("chloro", "氯", anchored=("*Cl",), paren=False),
    "bromo": RetainedSubstituent("bromo", "溴", anchored=("*Br",), paren=False),
    "iodo": RetainedSubstituent("iodo", "碘", anchored=("*I",), paren=False),
    "nitro": RetainedSubstituent("nitro", "硝基", anchored=("*[N+](=O)[O-]",), paren=False),
    "oxo": RetainedSubstituent("oxo", "氧代", anchored=("*=O",), paren=False),
    "nitro": RetainedSubstituent("nitro", "硝基", anchored=("*[N+](=O)[O-]",), paren=False),
    "isocyanato": RetainedSubstituent("isocyanato", "异氰酸基", anchored=("*N=C=O",), paren=False),
    "isothiocyanato": RetainedSubstituent("isothiocyanato", "异硫氰酸基", anchored=("*N=C=S",), paren=False),
    "methyl": RetainedSubstituent("methyl", "甲基", anchored=("*C",), paren=False),
    "methylidene": RetainedSubstituent("methylidene", "亚甲基", anchored=("*=C",), paren=False),
    "ethylidene": RetainedSubstituent("ethylidene", "亚乙基", anchored=("*=CC",), paren=False),
    "propylidene": RetainedSubstituent("propylidene", "亚丙基", anchored=("*=CCC",), paren=False),
    "cyclopropylidene": RetainedSubstituent("cyclopropylidene", "环丙亚基", anchored=("*=C1CC1",), paren=False),
    "cyclohexylidene": RetainedSubstituent("cyclohexylidene", "环己亚基", anchored=("*=C1CCCCC1",), paren=False),
    "sulfanylidene": RetainedSubstituent("sulfanylidene", "硫烷亚基", anchored=("*=S",), paren=False),
    "diaminomethylidene": RetainedSubstituent("diaminomethylidene", "二氨基亚甲基", anchored=("*C(=N)N",), paren=True),  # 脒/胍残基 C(=N)N：gold 按 P-66.1.1 取亚基式（diaminomethylideneamino），不取等价的 amino(imino)methylamino（两者互变、分子式相同，取测试集口径）。
    "ethyl": RetainedSubstituent("ethyl", "乙基", anchored=("*CC",), paren=False),
    "propyl": RetainedSubstituent("propyl", "丙基", anchored=("*CCC",), paren=False),
    "butyl": RetainedSubstituent("butyl", "丁基", anchored=("*CCCC",), paren=False),
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
    "methylsulfanyl": RetainedSubstituent("methylsulfanyl", "甲硫基", anchored=("*SC",), paren=False),
    "ethylsulfanyl": RetainedSubstituent("ethylsulfanyl", "乙硫基", anchored=("*SCC",), paren=False),
    "sulfanyl": RetainedSubstituent("sulfanyl", "巯基", anchored=("*S",), paren=False),
    "selanyl": RetainedSubstituent("selanyl", "氢硒基", anchored=("*[SeH]",), paren=False),
    "methylsulfinyl": RetainedSubstituent("methylsulfinyl", "甲基亚磺酰", anchored=("*S(C)=O",), paren=False),
    "methylsulfonyl": RetainedSubstituent("methylsulfonyl", "甲磺酰基", anchored=("*S(C)(=O)=O",), paren=False),
    "sulfo": RetainedSubstituent("sulfo", "磺基", anchored=("*S(=O)(=O)O",), paren=False),
    "tosyl": RetainedSubstituent("4-methylbenzenesulfonyl", "4-甲基苯磺酰基", anchored=("*S(=O)(=O)c1ccc(C)cc1",), paren=False),
    "carboxy": RetainedSubstituent("carboxy", "羧基", anchored=("*C(=O)O",), paren=False),
    "carbamoyl": RetainedSubstituent("carbamoyl", "氨基甲酰基", anchored=("*C(N)=O",), paren=False),  # P-66.1.1.4.1 氨基甲酸（carbamic acid）的酰基保留前缀；gold/ChEBI 全量 54 处取 carbamoyl，不取 aminocarbonyl/amino(oxo)methyl
    "carbamoylamino": RetainedSubstituent("carbamoylamino", "氨基甲酰氨基", anchored=("*NC(N)=O",), paren=True),  # P-66.1.1.6 ureido 在 IUPAC 已不推荐（P_1 附录：ureido/ureylene 不用），优选 carbamoylamino
    # 铵/𬭩型阳离子取代基（P-62.4.1：铵 azanium 去氢得 azaniumyl 型前缀）。gold 全量 174 例含 azanium*，        # 现行管线把这些带电 N 片段整体丢弃（no_coverage_gate / coverage_complete 误判），故按锚定叶子入表。
    "azaniumyl": RetainedSubstituent("azaniumyl", "铵基", anchored=("*[NH3+]",), paren=False),
    "methylazaniumyl": RetainedSubstituent("methylazaniumyl", "甲基铵基", anchored=("*[NH2+]C",), paren=False),
    "dimethylazaniumyl": RetainedSubstituent("dimethylazaniumyl", "二甲基铵基", anchored=("*[NH+](C)C",), paren=True),
    "trimethylazaniumyl": RetainedSubstituent("trimethylazaniumyl", "三甲基铵基", anchored=("*[N+](C)(C)C",), paren=True),
    "carbamoyloxy": RetainedSubstituent("carbamoyloxy", "氨基甲酰氧基", anchored=("*OC(N)=O",), paren=False),  # 氨基甲酸 O-酯残基（P-66.1.1.4.1）
    "carbamothioylamino": RetainedSubstituent("carbamothioylamino", "氨基硫代羰基氨基", anchored=("*NC(N)=S",), paren=True),  # 硫代氨基甲酸残基（P-66.1.1.4：carbamothioyl）
    "sulfamoyl": RetainedSubstituent("sulfamoyl", "氨磺酰基", anchored=("*S(N)(=O)=O",), paren=False),  # P-66.1.1.4.2 磺酰胺（sulfamoyl = H2N-SO2-）；N-取代时基名随取代基前移
    "phosphono": RetainedSubstituent("phosphono", "膦酸", anchored=("*P(=O)(O)O",), paren=False),  # P-102：phosphono 表示 -PO(OH)2，P 直连母体（对比 O 桥的 phosphonooxy）
    "phosphonato": RetainedSubstituent("phosphonato", "膦酸根", anchored=("*P(=O)([O-])O", "*P(=O)([O-])[O-]"), paren=False),  # P-102：phosphonato 表示 -PO(O-)2（单/双阴离子）
    "phosphonooxy": RetainedSubstituent("phosphonooxy", "膦酸氧基", anchored=("*OP(=O)(O)O",), paren=False),  # 磷酸降级前缀（P-67.1.5.1：羧酸等更高优先级 FG 存在时磷酸以 phosphonooxy 前缀表达）
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


def resolve_name(key: str) -> tuple[str, str]:
    """返回 registry 键对应的 (en, zh)。"""
    entry = _REGISTRY[key]
    return entry.en, entry.zh

def anchored_key(mol: Mol, atoms: frozenset[int], attach_old: int | None = None) -> str | None:
    """生成原子集的锚定 canonical-SMILES 键，失败时返回 None。"""
    return memo.by_key("anchored_key", (id(mol), atoms, attach_old),
                       lambda: _anchored_key_uncached(mol, atoms, attach_old), mol)


def _anchored_key_uncached(mol: Mol, atoms: frozenset[int], attach_old: int) -> str | None:
    """实际构建锚定子分子并取 canonical SMILES（无记忆版本，见 anchored_key）。"""
    from rdkit.Chem import MolToSmiles

    anchor = build_anchor_submol(mol, atoms, attach_old)
    return MolToSmiles(anchor) if anchor is not None else None

_WHOLE_ONLY_KEYS = frozenset({"*O", "*[O]", "*N"})  # 仅整分子顶层命中的锚定键：单原子杂原子自由基（表 2.1 去氢）。L3 取代基查表（anchored_lookup）跳过它们——游离 O/N 原子会被 _one_alkyl 误作侧链提取（酮羰基氧、酯氧、胺氮），命中 *O/*N 会错名成羟基/氨基。


def _table_hit(mol: Mol, atoms: frozenset[int], attach_old: int | None) -> str | None:
    """查 registry 锚定索引，返回 registry key；无命中返回 None。"""
    key = anchored_key(mol, atoms, attach_old)
    if key is None or key in _WHOLE_ONLY_KEYS:
        return None
    return _ANCHOR_INDEX.get(key)


def anchored_lookup(
    mol: Mol, atoms: frozenset[int], attach_old: int | None = None,
) -> tuple[str, str, bool] | None:
    """查找取代基原子集，返回 (en, zh, paren)；无命中返回 None。"""
    reg_key = _table_hit(mol, atoms, attach_old)
    if reg_key is None:
        return None
    en, zh = resolve_name(reg_key)
    return en, zh, _REGISTRY[reg_key].paren
