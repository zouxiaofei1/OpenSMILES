"""取代基锚定 canonical-SMILES 查表：registry 条目经 anchored 字段反查索引。

未命中回退完整命名路径；简单取代基（氟/氯/溴/碘/硝基/正构烷基等）与保留
取代基统一存于 registry；构建期校验 canonical 对拍与锚定键唯一性。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import Mol

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
    out: dict[str, RetainedSubstituent] = {
        # 基础取代基（原 ANCHOR_TABLE 并入 registry）
        "fluoro": RetainedSubstituent( "fluoro", "氟", "fluoro", "氟", P, anchored=("*F", ), paren=False, kind="halo", ),
        "chloro": RetainedSubstituent( "chloro", "氯", "chloro", "氯", P, anchored=("*Cl", ), paren=False, kind="halo", ),
        "bromo": RetainedSubstituent( "bromo", "溴", "bromo", "溴", P, anchored=("*Br", ), paren=False, kind="halo", ),
        "iodo": RetainedSubstituent( "iodo", "碘", "iodo", "碘", P, anchored=("*I", ), paren=False, kind="halo", ),
        "nitro": RetainedSubstituent( "nitro", "硝基", "nitro", "硝基", P, anchored=("*[N+](=O)[O-]", ), paren=False, kind="leaf", ),
        "isocyanato": RetainedSubstituent( "isocyanato", "异氰酸基", "isocyanato", "异氰酸基", P, anchored=("*N=C=O", ), paren=False, kind="leaf", ),
        "isothiocyanato": RetainedSubstituent( "isothiocyanato", "异硫氰酸基", "isothiocyanato", "异硫氰酸基", P, anchored=("*N=C=S", ), paren=False, kind="leaf", ),
        "methyl": RetainedSubstituent( "methyl", "甲基", "methyl", "甲基", P, anchored=("*C", ), paren=False, kind="alkyl", ),
        "ethyl": RetainedSubstituent( "ethyl", "乙基", "ethyl", "乙基", P, anchored=("*CC", ), paren=False, kind="alkyl", ),
        "propyl": RetainedSubstituent( "propyl", "丙基", "propyl", "丙基", P, anchored=("*CCC", ), paren=False, kind="alkyl", ),
        "butyl": RetainedSubstituent( "butyl", "丁基", "butyl", "丁基", P, anchored=("*CCCC", ), paren=False, kind="alkyl", ),
        # 支链 / 不饱和烷基
        "tert-butyl": RetainedSubstituent( "tert-butyl", "叔丁基", "tert-butyl", "叔丁基", P, anchored=("*C(C)(C)C", ), paren=False, kind="alkyl", ),
        "isopropyl": RetainedSubstituent( "propan-2-yl", "丙-2-基", "propan-2-yl", "丙-2-基", G, anchored=("*C(C)C", ), paren=False, kind="alkyl", ),
        "isobutyl": RetainedSubstituent("isobutyl", "异丁基", "2-methylpropyl", "2-甲基丙基", N, anchored=("*CC(C)C", ), paren=False, kind="alkyl", ),
        "sec-butyl": RetainedSubstituent("sec-butyl", "仲丁基", "butan-2-yl", "丁-2-基", N, anchored=("*C(C)CC", ), paren=False, kind="alkyl", ),
        "neopentyl": RetainedSubstituent( "neopentyl", "新戊基", "2,2-dimethylpropyl", "2,2-二甲基丙基", N, anchored=("*CC(C)(C)C", ), paren=False, kind="alkyl", ),
        "isopentyl": RetainedSubstituent( "isopentyl", "异戊基", "3-methylbutyl", "3-甲基丁基", N, anchored=("*CCC(C)C", ), paren=False, kind="alkyl", ),
        "vinyl": RetainedSubstituent( "vinyl", "乙烯基", "ethenyl", "乙烯基", G, anchored=("*C=C", ), paren=False, kind="alkyl", ),
        "allyl": RetainedSubstituent( "allyl", "烯丙基", "prop-2-en-1-yl", "丙-2-烯-1-基", G, anchored=("*CC=C", ), paren=False, kind="alkyl", ),
        "isopropenyl": RetainedSubstituent( "prop-1-en-2-yl", "异丙烯基", "prop-1-en-2-yl", "丙-1-烯-2-基", G, anchored=("*C(=C)C", ), paren=False, kind="alkyl", ),
        "propargyl": RetainedSubstituent( "propargyl", "炔丙基", "prop-2-yn-1-yl", "丙-2-炔-1-基", G, anchored=("*CC#C", ), paren=False, kind="alkyl", ),
        "benzyl": RetainedSubstituent( "benzyl", "苄基", "benzyl", "苄基", P, anchored=("*Cc1ccccc1", ), paren=False, kind="aryl", ),
        "methoxy": RetainedSubstituent( "methoxy", "甲氧基", "methoxy", "甲氧基", P, anchored=("*OC", ), paren=False, kind="leaf", ),
        "hydroperoxy": RetainedSubstituent( "hydroperoxy", "氢过氧基", "hydroperoxy", "氢过氧基", P, anchored=("*OO", ), paren=False, kind="leaf", ),
        "hydroxy": RetainedSubstituent( "hydroxy", "羟基", "hydroxy", "羟基", P, anchored=("*O", ), paren=False, kind="leaf", ),
        "oxidanyl": RetainedSubstituent( "oxidanyl", "氧基", "oxidanyl", "氧基", P, anchored=("*[O]", ), paren=False, kind="leaf", ),
        "methylsulfanyl": RetainedSubstituent( "methylsulfanyl", "甲硫基", "methylsulfanyl", "甲硫基", P, anchored=("*SC", ), paren=False, kind="leaf", ),
        "ethylsulfanyl": RetainedSubstituent( "ethylsulfanyl", "乙硫基", "ethylsulfanyl", "乙硫基", P, anchored=("*SCC", ), paren=False, kind="leaf", ),
        "sulfanyl": RetainedSubstituent( "sulfanyl", "硫烷基", "sulfanyl", "硫烷基", P, anchored=("*S", ), paren=False, kind="leaf", ),
        "selanyl": RetainedSubstituent( "selanyl", "硒烷基", "selanyl", "硒烷基", P, anchored=("*[SeH]", ), paren=False, kind="leaf", ),
        "methylsulfinyl": RetainedSubstituent( "methylsulfinyl", "甲亚磺酰基", "methanesulfinyl", "甲亚磺酰基", P, anchored=("*S(C)=O", ), paren=False, kind="leaf", ),
        "methylsulfonyl": RetainedSubstituent( "methylsulfonyl", "甲磺酰基", "methanesulfonyl", "甲磺酰基", P, anchored=("*S(C)(=O)=O", ), paren=False, kind="leaf", ),
        "sulfo": RetainedSubstituent( "sulfo", "磺基", "sulfo", "磺基", P, anchored=("*S(=O)(=O)O", ), paren=False, kind="leaf", ),
        "tosyl": RetainedSubstituent( "tosyl", "对甲苯磺酰基", "4-methylbenzenesulfonyl", "4-甲基苯磺酰基", N, anchored=("*S(=O)(=O)c1ccc(C)cc1", ), paren=False, kind="leaf", ),
        "carboxy": RetainedSubstituent( "carboxy", "羧基", "carboxy", "羧基", P, anchored=("*C(=O)O", ), paren=False, kind="leaf", ),
        "carboxymethyl": RetainedSubstituent( "carboxymethyl", "羧甲基", "carboxymethyl", "羧甲基", P, anchored=("*CC(=O)O", ), paren=False, kind="leaf", ),
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


def get_retained(key: str) -> RetainedSubstituent | None:
    """按 registry 键查找保留取代基。"""
    return _REGISTRY.get(key)


def resolve_name(key: str, *, name_mode: str = "general") -> tuple[str, str]:
    """返回给定命名模式下 registry 键对应的 (en, zh)。

    - "general"：始终为保留名/常用名
    - "pin"：除非条目本身为 PIN 级，否则用系统名
    """
    entry = _REGISTRY[key]
    
    # if name_mode == "pin" and entry.level != IupacLevel.PIN:
    #     return entry.systematic_en, entry.systematic_zh
    return entry.systematic_en, entry.systematic_zh


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


# 仅整分子顶层命中的锚定键：单原子杂原子自由基（表 2.1 去氢）。
# L3 取代基查表（anchored_entry）跳过它们——游离 O/N 原子会被 _one_alkyl
# 误作侧链提取（酮羰基氧、酯氧、胺氮），命中 *O/*N 会错名成羟基/氨基。
_WHOLE_ONLY_KEYS = frozenset({"*O", "*[O]", "*N"})


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
    """查找取代基原子集；在 name_mode 下返回 (en, zh, paren)。

    registry 键条目通过文件内 resolve_name 解析，因此 pin 模式对 isopropyl 会生成 propan-2-yl 等。无命中时返回 None。
    """
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
