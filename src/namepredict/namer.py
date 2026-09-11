"""顶层命名管线：SMILES → L1–L5 双语 IUPAC，带缓存、候选重试与盐拆分。"""

from __future__ import annotations

import copy
import time

from namepredict.tools import memo
from namepredict.cache.common_names import CommonNameCache
from namepredict.constants import N_PREFIX_KINDS
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer0.salt import dissociate_salt
from namepredict.layer1.analyzer import analyze
from namepredict.layer3.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import select_parent_tied
from namepredict.layer3.coverage import build_coverage_ledger
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer3.substituent_namer import SubstituentName
from namepredict.layer4.candidate_keys import prefix_locant_set, suffix_locant_set
from namepredict.layer4.numbering import number
from namepredict.layer5.assembler import assemble
from namepredict.tools.anchored_table import anchored_whole_mol
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


def _chain_meta(numbered: dict) -> dict:
    """从编号结果提取母体链、母体 kind 与 S 桥自含围栏标记（供 L3 判断前缀是否还需加括号）。"""
    parent = numbered.get("parent") or {}
    return {"parent_chain": list(parent.get("chain") or []), "parent_kind": parent.get("kind"),
            "bridge_self_enclosed": bool(numbered.get("bridge_self_enclosed"))}


def _claim_from_sub(s: dict, atoms: frozenset[int]) -> ClaimedBlock:
    """由取代基 dict 与原子集构建 ClaimedBlock 归属块。"""
    attach = s.get("attach_idx")
    return ClaimedBlock(
        slot=SideSlot.OTHER,
        attach_parent=int(attach) if attach is not None else -1,
        root=min(atoms),
        atoms=atoms,
    )


def _one_name_from_sub(s: dict) -> SubstituentName | None:
    """将单个取代基 dict 转为 SubstituentName，无原子则返回 None。"""
    atoms = frozenset(s.get("atoms") or [])
    if not atoms:
        return None
    return SubstituentName(
        claim=_claim_from_sub(s, atoms),
        en=s.get("en") or "x",
        zh=s.get("zh") or "x",
        requires_parentheses=bool(s.get("paren")),
        backend=s.get("backend") or "extract",
    )


def _names_from_subs(subs: list[dict]) -> list[SubstituentName]:
    """批量将取代基 dict 列表转换为 SubstituentName 列表（跳过无效项）。"""
    return [n for s in subs if (n := _one_name_from_sub(s)) is not None]


def _ledger_complete(mol, owned, subst: list[dict]) -> bool:
    """基于 coverage ledger 判断取代基是否覆盖全部母体原子。"""
    names = _names_from_subs(subst)
    return build_coverage_ledger(mol, owned_atoms=owned, names=names).complete


def _ok_result(numbered: dict, *, depth: int, t0: float, name_mode: str = "general") -> NameResult | None:
    """组装编号结果为 NameResult，成功且非空才返回；meta 附链元数据与母体取代基数（供 L3 判定词干是否复合）。"""
    numbered["name_mode"] = name_mode
    result = assemble(numbered, time_ms=_elapsed_ms(t0))
    if not result.success or not result.en:
        return None
    result.meta = {**(result.meta or {}), **_chain_meta(numbered),
                   "depth": depth, "coverage_complete": True,
                   "parent_substituent_count": len(numbered.get("substituents") or [])}
    return result


def _chain_set(parent: dict) -> set[int]:
    """取出母体的链原子索引集合。"""
    return set(parent.get("chain") or [])


def _remap_attach(parent: dict, s: dict) -> dict:
    """确保 attach_idx 位于母体链上，供 L4 orient 使用（环官能团连接）。"""
    if s.get("o_side"):
        return s  # ester O 侧烷基：连接点保留在酯 O 上，不做链重映射
    chain = _chain_set(parent)
    attach = s.get("attach_idx")
    if attach in chain:
        return s
    for key in ("ring_attach_idx", "amide_c_idx", "amine_c_idx", "ketone_c_idx", "oh_c_idx"):
        alt = parent.get(key)
        if alt in chain:
            return {**s, "attach_idx": alt}
    return s

_MAX_TIED_CANDIDATES = 4  # P-45.2.2 需要为每个并列候选各跑一次 L3–L5，上限防止组合爆炸（benchmark 中并列组多为 2–4 个）。


def _subs_for_numbering(parent: dict, subst: list[dict]) -> list[dict]:
    """筛选并重映射参与编号的取代基（O 侧、链上连接点与 N 端）；N- 取代基（n_* kind）位次隐含省略但仍保留进 L5 前缀组装。"""
    chain = _chain_set(parent)
    out: list[dict] = []
    for s in subst:
        s2 = _remap_attach(parent, s)
        if s2.get("o_side") or s2.get("attach_idx") in chain or s2.get("kind") in N_PREFIX_KINDS:
            out.append(s2)
    return out


def _assemble_candidate(parent, subst, *, depth: int, t0: float, name_mode: str = "general") -> NameResult | None:
    """对单个候选执行编号+组装，编号异常或失败时返回 None（meta 附 P-44.1.1 / P-45.2.2 位次键）。"""
    try:
        numbered = number(parent, _subs_for_numbering(parent, subst))
    except (ValueError, KeyError, TypeError):
        return None
    hit = _ok_result(numbered, depth=depth, t0=t0, name_mode=name_mode)
    if hit is not None:
        hit.meta = {**(hit.meta or {}), "p44_1_1_key": suffix_locant_set(numbered),
                    "p45_2_2_key": prefix_locant_set(numbered)}
    return hit


def _prepare_candidate(
    info: dict, parent: dict, *, name_mode: str = "general", cache: CommonNameCache | None = None,
    depth: int = 0,
) -> tuple[dict, list[dict], bool]:
    """完成母体归属、提取取代基并返回 (parent, subst, complete)；depth 透传给取代基命名。"""
    mol = info["mol"]
    parent = finalize_parent_ownership(parent, mol)
    if not parent.get("owned_atoms"):
        return parent, [], False
    if not parent.get("chain") and not info.get("has_ring"):
        return parent, [], False
    subst = extract_substituents(info, parent, name_mode=name_mode, cache=cache, depth=depth)
    complete = _ledger_complete(mol, parent["owned_atoms"], subst)
    return parent, subst, complete


def try_candidate(
    info: dict,
    parent: dict,
    *,
    depth: int = 0,
    t0: float | None = None,
    require_complete: bool = True,
    name_mode: str = "general",
) -> NameResult | None:
    """完成归属、提取取代基、可选要求完整 coverage ledger，然后组装。"""
    t0 = t0 if t0 is not None else time.perf_counter()
    parent, subst, complete = _prepare_candidate(info, parent, name_mode=name_mode, depth=depth)
    if require_complete and not complete:
        return None
    hit = _assemble_candidate(parent, subst, depth=depth, t0=t0)
    if hit is None:
        return None
    hit.meta = {**(hit.meta or {}), "coverage_complete": complete}
    return hit


def _candidate_key(hit: NameResult) -> tuple:
    """候选裁决键：(P-44.1.1 后缀位次集合, P-45.2.2 前缀位次集合)。"""
    meta = hit.meta or {}
    return meta.get("p44_1_1_key") or (), meta.get("p45_2_2_key") or ()


def _best_hit(hits: list[tuple]) -> NameResult | None:
    """候选裁决：P-44.1.1 后缀位次集合未决（并列）时才按 P-45.2.2 前缀位次集合取最小，否则保持候选顺序。"""
    if not hits:
        return None
    if len({suffix for suffix, _, _, _ in hits}) > 1:
        return hits[0][3]  # P-44.1.1 已决：保持候选顺序（当前即首选）
    return min(hits, key=lambda t: (t[1], t[2]))[3]


def _complete_hit(prepared, *, depth, t0, name_mode):
    """在候选集中寻找 coverage 完整且可组装的命中（并列候选取裁决键最小者）。"""
    hits = []
    for order, (parent, subst, complete) in enumerate(prepared):
        if not complete:
            continue
        hit = _assemble_candidate(parent, subst, depth=depth, t0=t0, name_mode=name_mode)
        if hit is not None:
            hit.meta = {**(hit.meta or {}), "coverage_complete": True}
            hits.append((*_candidate_key(hit), order, hit))
    return _best_hit(hits)


def _partial_hit(prepared, *, depth, t0, name_mode, attempts):
    """放宽 coverage 门控，在候选集中找可组装的命中（并列候选取裁决键最小者）。"""
    hits = []
    for order, (parent, subst, complete) in enumerate(prepared):
        hit = _assemble_candidate(parent, subst, depth=depth, t0=t0, name_mode=name_mode)
        if hit is not None:
            hit.meta = {**(hit.meta or {}), "coverage_complete": complete, "fallback": "no_coverage_gate", "attempts": attempts}
            hits.append((*_candidate_key(hit), order, hit))
    return _best_hit(hits)


def _try_phase(prepared, *, depth, t0, name_mode, attempts):
    """依次尝试完整命中与部分命中两阶段。"""
    hit = _complete_hit(prepared, depth=depth, t0=t0, name_mode=name_mode)
    return hit or _partial_hit(
        prepared, depth=depth, t0=t0, name_mode=name_mode, attempts=attempts,
    )


def _candidate_phases(info: dict, depth: int) -> list[list[dict]]:
    """选取候选母体阶段列表（P-45.2.1 并列组，上限 _MAX_TIED_CANDIDATES 个）。"""
    group = select_parent_tied(info)
    return [group[:_MAX_TIED_CANDIDATES]] if group else [[]]


def _run_candidates(
    info: dict, *, depth: int, t0: float, name_mode: str = "general", cache: CommonNameCache | None = None,
) -> NameResult:
    """仅尝试 P-44.1.1 高优先级阶段；绝不降级能力。"""
    attempts: list[dict] = []
    phase = _candidate_phases(info, depth)[0]
    prepared = [_prepare_candidate(info, cand, name_mode=name_mode, cache=cache, depth=depth) for cand in phase]
    hit = _try_phase(prepared, depth=depth, t0=t0, name_mode=name_mode, attempts=attempts)
    return hit or _fail(_elapsed_ms(t0), "no_assemblable_candidate", attempts=attempts)


def _apply_salt_suffix(result: NameResult, salt: dict) -> NameResult:
    """将盐元数据组装为名称后缀（碱金属盐/HCl 加成盐），仅成功结果生效；kind=phosphate 由 L5 phosphate worker 自行组装，此处跳过。"""
    if not result.success or not salt:
        return result
    if (result.meta or {}).get("parent_kind") == "phosphate":
        return result
    from namepredict.layer5.stems import maybe_metal_salt_names

    en, zh = maybe_metal_salt_names({"salt": salt}, result.en or "", result.zh or "")
    if en == (result.en or "") and zh == (result.zh or ""):
        return result
    out = copy.copy(result)
    out.en, out.zh = en, zh
    return out


def _name_mol(
    mol,
    *,
    depth: int = 0,
    cache: CommonNameCache | None = None,
    t0: float | None = None,
    name_mode: str = "general",
    root_ctx: tuple | None = None,
) -> NameResult:
    """从 mol 运行 L1–L5，带 coverage 门控的候选重试；root_ctx=(根分子, 原子→根索引映射) 供取代基 R/S 回根分子重算。"""
    t0 = t0 if t0 is not None else time.perf_counter()
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse")
    organic, salt = dissociate_salt(mol)
    if root_ctx is None:  # 顶层整分子
        root_mol, to_root = organic, list(range(organic.GetNumAtoms()))
        run_cache = CommonNameCache(max_entries=2000) if cache is not None else None  # 顶层整分子：内部 `*` 片段名仅在本分子运行内共享（同根立体一致）；跨根缓存会带入别的宿主异头立体，须禁用。
    else:
        root_mol, to_root = root_ctx  # 无盐时 organic 即 mol、索引不变；锚定碎片必为单片段不含盐，映射直接沿用。
        run_cache = cache
    info = analyze(organic)
    info["root_ctx"] = (root_mol, to_root)
    info["salt"] = salt  # 磷酸母体 producer 的盐门控与 salt_meta 来源
    result = _run_candidates(info, depth=depth, t0=t0, name_mode=name_mode, cache=run_cache)
    result = _apply_salt_suffix(result, salt)
    if salt and result.success:
        result.meta = {**(result.meta or {}), "salt": salt}
    return result


def _pipeline(smiles: str, t0: float, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> tuple[NameResult, "Mol | None"]:
    """预处理 SMILES 后进入 mol 命名流程，解析失败返回失败结果；整分子命中锚定表（带 * 锚点输入本身即锚定键）直接返回保留名，免经自由基母体管线。同时回传解析出的 mol，供调用方复用（免去二次解析）。"""
    mol = preprocess(smiles)
    if mol is None:
        return _fail(_elapsed_ms(t0), "parse"), None
    whole = anchored_whole_mol(mol, name_mode=name_mode)
    if whole is not None:
        en, zh, paren, kind = whole
        return NameResult(
            en=en, zh=zh, success=True, source="anchored",
            time_ms=_elapsed_ms(t0),
            meta={"parent_kind": "radical", "anchored": True},
        ), mol
    return _name_mol(mol, depth=0, t0=t0, name_mode=name_mode, cache=cache), mol


def _cache_put(cache: CommonNameCache, smiles: str, result: NameResult) -> None:
    """写缓存，容量满的 ValueError 静默忽略。"""
    try:
        cache.put(smiles, result)
    except ValueError:
        pass


def _canonical_result(mol, result: NameResult) -> NameResult:
    """复制结果并把 meta.parent_chain 重写为规范原子排序——缓存以子结构 SMILES 为键，parent_chain 不能依赖随母体变化的切分点原子顺序。"""
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


def _name_uncached(smiles: str, t0: float, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> tuple[NameResult, "Mol | None"]:
    """缓存未命中时直接走完整命名管线，回传解析出的 mol 供调用方复用。"""
    return _pipeline(smiles, t0, name_mode=name_mode, cache=cache)


class SMILESNNamer:
    """SMILES → IUPAC 命名的顶层命名器（含缓存与命名模式）。"""

    def __init__(self, cache: CommonNameCache | None = None, *, name_mode: str = "general") -> None:
        """初始化命名器，未提供缓存则构造默认 20000 条容量的缓存。"""
        self.cache = cache if cache is not None else CommonNameCache(max_entries=20000)  # 容量留足给主分子 + 递归子结构命名（全量去重后约 8.7k 条）
        self._name_mode = name_mode

    def name(self, smiles: str) -> NameResult:
        """命名单个 SMILES；先查缓存，未命中则计算并写回成功结果。"""
        t0 = time.perf_counter()
        hit = self.cache.get(smiles)
        if hit is not None:
            return hit
        memo.begin_run()  # 本次命名的中间结果记忆：不跨分子共享，见 cache/memo
        result, mol = _name_uncached(smiles, t0, name_mode=self._name_mode, cache=self.cache)
        if result.success and mol is not None:
            _cache_put(self.cache, smiles, _canonical_result(mol, result))
        return result
