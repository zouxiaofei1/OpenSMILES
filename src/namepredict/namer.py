"""顶层命名管线：SMILES → L1–L5 双语 IUPAC（含缓存）。"""

from __future__ import annotations

import copy
import time
from collections import Counter

from rdkit import Chem

from namepredict.tools import memo
from namepredict.tools.common_names import CommonNameCache
from namepredict.constants import (
    BIS_EN, BIS_ZH, N_PREFIX_KINDS, OXO_CENTER_KINDS, SIMPLE_MAX_HEAVY, SIMPLE_MOLECULES,
)
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import _frag_role, dissociate_salt, pair_unique_ions
from namepredict.layer1.analyzer import analyze
from namepredict.tools.re import component_order_key
from namepredict.layer2.parent_select import finalize_parent_ownership, select_parent
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.locant_calc import prefix_locant_set, suffix_locant_set
from namepredict.layer4.numbering import number
from namepredict.layer5.assembler import assemble
from namepredict.types import NameResult


def _fail(time_ms: float = 0.0, reason: str = "parse", **meta) -> NameResult:
    """构造一个失败命名结果并附带原因与元数据。"""
    return NameResult(
        en="", zh="", success=False, source="iupac",
        time_ms=time_ms, meta={"reason": reason, **meta},
    )


def _elapsed_ms(t0: float) -> float:
    """计算自 t0 起已消耗的毫秒数。"""
    return (time.perf_counter() - t0) * 1000.0


def _simple_species(mol) -> tuple[str, str] | None:
    """≤3 重原子的单质/简单分子查表（constants.SIMPLE_MOLECULES）；未命中返回 None。"""
    if mol.GetNumHeavyAtoms() > SIMPLE_MAX_HEAVY:
        return None
    return SIMPLE_MOLECULES.get(Chem.MolToSmiles(mol))


def _chain_meta(numbered: dict) -> dict:
    """从编号结果提取母体链、kind 与 S 桥自含围栏标记。"""
    parent = numbered.get("parent") or {}
    return {"parent_chain": list(parent.get("chain") or []), "parent_kind": parent.get("kind"),
            "parent_labels": _label_list(parent),
            "bridge_self_enclosed": bool(numbered.get("bridge_self_enclosed"))}

def _label_list(parent: dict) -> list:
    """母体整体编号标签（稠环桥头 3a/6a）；长度不符时返回空表。"""
    labels = (parent.get("numbering_scaffold") or {}).get("labels") or []
    chain = parent.get("chain") or []
    return [int(x) if str(x).isdigit() else str(x) for x in labels] if len(labels) == len(chain) else []


def _ok_result(numbered: dict, *, t0: float) -> NameResult | None:
    """组装编号结果，成功且非空才返回（meta 附链元数据）。"""
    result = assemble(numbered, time_ms=_elapsed_ms(t0))
    if not result.success or not result.en:
        return None
    result.meta = {**(result.meta or {}), **_chain_meta(numbered),
                   "parent_substituent_count": len(numbered.get("substituents") or [])}
    return result


def _chain_set(parent: dict) -> set[int]:
    """取出母体的链原子索引集合。"""
    return set(parent.get("chain") or [])


def _remap_attach(parent: dict, s: dict) -> dict:
    """确保 attach_idx 位于母体链上，供 L4 orient 使用。"""
    if s.get("o_side"):
        return s  # ester O 侧烷基：连接点保留在酯 O 上，不做链重映射
    chain = _chain_set(parent)
    attach = s.get("attach_idx")
    if attach in chain:
        return s
    for alt in _remap_candidates(parent):
        if alt in chain:
            return {**s, "attach_idx": alt}
    return s


def _remap_candidates(parent: dict) -> list[int]:
    """可重挂的母体锚点：环附着原子，其次单锚点官能团锚点。"""
    ring = parent.get("ring_attach_idx")
    facts = parent.get("principal_expression_facts")
    out = [ring] if ring is not None else []
    if facts is not None and len(facts.anchor_atoms) == 1:
        out.extend(sorted(facts.anchor_atoms))
    return out




def _subs_for_numbering(parent: dict, subst: list[dict]) -> list[dict]:
    """筛选并重映射参与编号的取代基（O 侧/链上/N 端）。"""
    chain = _chain_set(parent)
    out: list[dict] = []
    for s in subst:
        s2 = _remap_attach(parent, s)
        if s2.get("o_side") or s2.get("attach_idx") in chain or s2.get("kind") in N_PREFIX_KINDS:
            out.append(s2)
    return out


def _assemble_candidate(parent, subst, *, t0: float) -> NameResult | None:
    """对单个候选执行编号+组装，失败返回 None（meta 附位次键）。"""
    try:
        numbered = number(parent, _subs_for_numbering(parent, subst))
    except (ValueError, KeyError, TypeError):
        return None
    hit = _ok_result(numbered, t0=t0)
    if hit is not None:
        hit.meta = {**(hit.meta or {}), "p44_1_1_key": suffix_locant_set(numbered),
                    "p45_2_2_key": prefix_locant_set(numbered)}
    return hit


def _prepare_candidate(
    info: dict, parent: dict, *, cache: CommonNameCache | None = None,
) -> tuple[dict, list[dict]]:
    """完成母体归属与取代基提取，返回 (母体, 取代基表)。"""
    mol = info["mol"]
    parent = finalize_parent_ownership(parent, mol)
    if not parent.get("owned_atoms"):
        return parent, []
    if not parent.get("chain") and not info.get("has_ring"):
        return parent, []
    return parent, extract_substituents(info, parent, cache=cache)

def _candidate_key(hit: NameResult) -> tuple:
    """候选裁决键：(P-44.1.1 后缀位次, P-45.2.2 前缀位次)。"""
    meta = hit.meta or {}
    return meta.get("p44_1_1_key") or (), meta.get("p45_2_2_key") or ()


def _best_hit(hits: list[tuple]) -> NameResult | None:
    """候选裁决：P-44.1.1 未决时按 P-45.2.2 前缀位次取最小。"""
    if not hits:
        return None
    if len({suffix for suffix, _, _, _ in hits}) > 1:
        return hits[0][3]  # P-44.1.1 已决：保持候选顺序（当前即首选）
    return min(hits, key=lambda t: (t[1], t[2]))[3]


def _try_phase(prepared, *, t0):
    """L4+L5入口：逐候选装配，取裁决最优者。"""
    hits = []
    for order, (parent, subst) in enumerate(prepared):
        hit = _assemble_candidate(parent, subst, t0=t0)
        if hit is not None:
            hits.append((*_candidate_key(hit), order, hit))
    return _best_hit(hits)


def _run_candidates(
    info: dict, *, t0: float, cache: CommonNameCache | None = None,
) -> NameResult:
    """_run_candidates"""
    phase = select_parent(info) or []#Layer2入口
    prepared = [_prepare_candidate(info, cand, cache=cache) for cand in phase]#L3
    hit = _try_phase(prepared, t0=t0)#L4入口
    return hit or _fail(_elapsed_ms(t0), "no_assemblable_candidate")


def _apply_salt_suffix(result: NameResult, salt: dict) -> NameResult:
    """把盐元数据组装为名称后缀（碱金属盐/HCl 盐），仅成功结果生效。"""
    if not result.success or not salt:
        return result
    if (result.meta or {}).get("parent_kind") in OXO_CENTER_KINDS:  # 含氧酸中心母体自带盐组装
        return result
    from namepredict.layer5.stems import join_metal_salt_names

    en, zh = join_metal_salt_names({"salt": salt}, result.en or "", result.zh or "")
    if en == (result.en or "") and zh == (result.zh or ""):
        return result
    out = copy.copy(result)
    out.en, out.zh = en, zh
    return out


COMPONENT_SEP = "; "  # 无法配对的剩余组分之间的连接符


def _name_component(frag, *, cache) -> tuple[str, str]:
    """单片段组分名：盐词表 → 独立命名 → 规范 SMILES 兜底（绝不丢弃）。"""
    role = _frag_role(frag)
    if role is not None:
        return role[1], role[2]
    sub = _name_mol(frag, cache=cache)
    if sub.success and (sub.en or "").strip() and (sub.zh or "").strip():
        return sub.en, sub.zh
    smi = Chem.MolToSmiles(frag)
    return smi, smi


def _merge_repeats(names: list[tuple[str, str]]) -> tuple[list[str], list[str]]:
    """重复组分加倍数围栏（P-16.3.2），返回 (en 词条, zh 词条)。"""
    counts = Counter(n_en for n_en, _ in names)
    seen: set[str] = set()
    out_en, out_zh = [], []
    for n_en, n_zh in names:
        if n_en in seen:
            continue
        seen.add(n_en)
        k = counts[n_en]
        b_en, b_zh = BIS_EN.get(k), BIS_ZH.get(k)
        if k > 1 and b_en and b_zh:
            out_en.append(f"{b_en}({n_en})")
            out_zh.append(f"{b_zh}({n_zh})")
            continue
        out_en.extend([n_en] * k)
        out_zh.extend([n_zh] * k)
    return out_en, out_zh


def _name_components(mol, *, cache) -> tuple[str, str] | None:
    """多片段组分名拼装：唯一配对阴阳离子走 P-77 二元名，其余按字母序以 '; ' 连接。"""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return None  # 单片段由常规管线命名（防递归）
    pair = pair_unique_ions(frags)
    if pair is not None:  # P-77.1.1 二元盐名：阳离子名 + 阴离子名
        cat, an, n_cat, n_an = pair
        c_en, c_zh = _name_component(cat, cache=cache)
        a_en, a_zh = _name_component(an, cache=cache)
        parts = _merge_repeats([(c_en, c_zh)] * n_cat + [(a_en, a_zh)] * n_an)
        return " ".join(parts[0]), " ".join(parts[1])
    names = [_name_component(f, cache=cache) for f in frags]
    names.sort(key=lambda t: component_order_key(t[0]))  # 组分名字母序（忽略位次/括号/倍增前缀）
    parts = _merge_repeats(names)
    return COMPONENT_SEP.join(parts[0]), COMPONENT_SEP.join(parts[1])


def _name_mol(
    mol,
    *,
    cache: CommonNameCache | None = None,
    t0: float | None = None,
    root_ctx: tuple | None = None,
) -> NameResult:
    """从 mol 运行 L1–L5：解盐、分析、候选装配、盐后缀。"""
    t0 = t0 if t0 is not None else time.perf_counter()
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    simple = _simple_species(mol)  # 单元素/极简单分子无母体链环，直接在入口给出保留名
    if simple is not None:
        return NameResult(en=simple[0], zh=simple[1], success=True, source="iupac",
                          time_ms=_elapsed_ms(t0), meta={"reason": "simple_species"})
    organic, salt = dissociate_salt(mol)
    if root_ctx is None:  # 顶层整分子
        root_mol, to_root = organic, list(range(organic.GetNumAtoms()))
        run_cache = CommonNameCache(max_entries=2000) if cache is not None else None  # 顶层整分子：内部 `*` 片段名仅在本分子运行内共享（同根立体一致）；跨根缓存会带入别的宿主异头立体，须禁用。
    else:
        root_mol, to_root = root_ctx  # 无盐时 organic 即 mol、索引不变；锚定碎片必为单片段不含盐，映射直接沿用。
        run_cache = cache
    info = analyze(organic) # 进入Layer1
    info["root_ctx"] = (root_mol, to_root)
    info["salt"] = salt  # 磷酸母体 producer 的盐门控与 salt_meta 来源
    result = _run_candidates(info, t0=t0, cache=run_cache)
    result = _apply_salt_suffix(result, salt)
    # 多片段体系：配对阴阳离子成盐，其余组分按字母序拼接，绝不丢弃。
    # 已配盐但有机部分组装失败时同样回落，否则整分子静默变空名。
    if not salt or not result.success:
        joined = _name_components(mol, cache=run_cache)
        if joined is not None and joined[0].strip():
            result = copy.copy(result)
            result.en, result.zh = joined
            result.success = True
    if salt and result.success:
        result.meta = {**(result.meta or {}), "salt": salt}
    return result


def _pipeline(smiles: str, t0: float, *, cache: CommonNameCache | None = None) -> tuple[NameResult, "Mol | None"]:
    """预处理 SMILES 后进入 mol 命名流程，并回传解析出的 mol。"""
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse"), None
    return _name_mol(mol, t0=t0, cache=cache), mol


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    """写缓存，容量满的 ValueError 静默忽略。"""
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


def _canonical_result(mol, result: NameResult) -> NameResult:
    """复制结果并把 meta.parent_chain 重写为规范原子排序。"""
    chain = (result.meta or {}).get("parent_chain") or []
    if not chain:
        return result
    try:
        from rdkit.Chem import CanonicalRankAtoms

        ranks = CanonicalRankAtoms(mol)
        canon = [ranks[i] for i in chain]
    except Exception:
        return result
    r = copy.copy(result)
    r.meta = {**(result.meta or {}), "parent_chain": canon}
    return r


class SMILESNNamer:
    """SMILES → IUPAC 命名的顶层命名器（含缓存）。"""

    def __init__(self, cache: CommonNameCache | None = None) -> None:
        """初始化命名器，未提供缓存则构造默认 20000 条容量的缓存。"""
        self.cache = cache if cache is not None else CommonNameCache(max_entries=20000)  # 容量留足给主分子 + 递归子结构命名（全量去重后约 8.7k 条）

    def name(self, smiles: str) -> NameResult:
        """命名单个 SMILES；先查缓存，未命中则计算并写回成功结果。"""
        t0 = time.perf_counter()
        hit = self.cache.get(smiles)
        if hit is not None:
            return hit
        memo.begin_run()  # 本次命名的中间结果记忆：不跨分子共享，见 cache/memo
        result, mol = _pipeline(smiles, t0, cache=self.cache)  # 缓存未命中：走完整管线并回传 mol 供写缓存
        if result.success and mol is not None:
            _cache_put(self.cache, smiles, _canonical_result(mol, result))
        return result
