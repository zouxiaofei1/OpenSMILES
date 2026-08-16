"""name_as_substituent：切割 submol → free-name 管道 → P-29 -yl 形式。
无宿主（苯甲酰胺 gate / n_block 提取）接线——纯 cut→pipeline→yl。
"""
from __future__ import annotations

from rdkit.Chem import CanonicalRankAtoms
from namepredict.constants import C
from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.submol_build import build_anchor_submol, build_cut_submol
from namepredict.tools.free_to_yl import free_to_yl as yl_form

def _radical_yl_from_sub(
    mol, atoms: frozenset, attach_old: int, *, depth: int, name_mode: str,
    cache: CommonNameCache | None,
) -> tuple[str, str, bool] | None:
    """碳连接点：锚定 * 走 radical 主基团管线，L4 权威位次 + P-22.2.4 保留名。

    取代 _is_arom_c6/_arene_yl_from_sub 的 6 元全碳芳环特判：所有碳连接点子结构
    （芳基/杂环/稠环/开链/烯）统一锚定，radical worker 直接拼 -yl 名。
    """
    from namepredict.namer import _cache_put, _canonical_result, _name_mol
    from rdkit import Chem
    anchored = build_anchor_submol(mol, atoms, attach_old)
    if anchored is None:
        return None
    smiles = Chem.MolToSmiles(anchored)
    hit = cache.get(smiles) if cache is not None else None
    if hit is None:
        hit = _name_mol(anchored, depth=depth, name_mode=name_mode, cache=cache)
        if cache is not None and hit.success and hit.en:
            hit = _canonical_result(anchored, hit)
            _cache_put(cache, smiles, hit)
    if not hit.success or not hit.en:
        return None
    if not (hit.meta or {}).get("parent_kind") == "radical":
        # 锚定分子必被 L1 radical 条目检出、principal 必选（p41=1），理论不可达，防御。
        return None
    return hit.en, hit.zh, hit.en != "phenyl"


def _locant_from_result(mol, result, attach_new: int) -> int | None:
    """parent_chain 上的连接位次（canonical rank；跨切割缓存安全）。"""
    chain = (result.meta or {}).get("parent_chain") or []
    if not chain:
        return None
    ranks = CanonicalRankAtoms(mol)
    ar = ranks[attach_new]
    if ar in chain:
        return chain.index(ar) + 1
    return None


def _locant_via_hetero(mol, attach_new: int, result) -> int | None:
    """当连接原子是不在链中的杂原子（O/S/N）时，使用其链上的邻接原子。"""
    a = mol.GetAtomWithIdx(attach_new)
    if a.GetAtomicNum() not in (8, 16, 7):
        return None
    chain = (result.meta or {}).get("parent_chain") or []
    if not chain:
        return None
    ranks = CanonicalRankAtoms(mol)
    c_nbrs = [ranks[n.GetIdx()] for n in a.GetNeighbors() if ranks[n.GetIdx()] in chain]
    return chain.index(c_nbrs[0]) + 1 if len(c_nbrs) == 1 else None


def _yl_from_sub(
    sub, *, mol, atoms, attach_old, depth: int, name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    """连接点类型分派：碳→锚定 radical 优先；非碳/锚定失败→H 封端 free-name + free_to_yl。"""
    from namepredict.namer import _cache_put, _canonical_result, _name_mol
    from rdkit import Chem

    # ① 碳连接点：锚定 * 直接走 radical 管线（信息无损——L4 权威位次、P-22.2.4 保留名、
    #    且 FG 型碳链如 HOCH2CH2- 不再被 free_to_yl 误前缀化成 ethoxy）。
    if mol.GetAtomWithIdx(attach_old).GetAtomicNum() == C:
        hit = _radical_yl_from_sub(mol, atoms, attach_old, depth=depth,
                                   name_mode=name_mode, cache=cache)
        if hit is not None:
            return hit
    # ② 非碳连接点（卤素/烷氧基/氨基/硫基）/ 锚定失败：H 封端子分子 free-name →
    #    free_to_yl FG 前缀化（甲醇→甲氧基、乙醇→乙氧基），含跨 cut 缓存。
    result = None
    smiles = None
    if cache is not None:
        # 递归子结构命名与主分子共享缓存：命中则跳过整条 L2 递归。
        smiles = Chem.MolToSmiles(sub.mol)
        result = cache.get(smiles)
    if result is None:
        # 用 canonical SMILES 重解析再命名：命名输入与缓存 key 严格一一对应，消除不同 cut 上下文（环断点/手性方向）泄漏进命名的差异。
        named_mol = sub.mol
        if smiles:
            canonical = Chem.MolFromSmiles(smiles)
            if canonical is not None:
                named_mol = canonical
        # print(smiles)
        result = _name_mol(named_mol, depth=depth, name_mode=name_mode, cache=cache)
        if cache is not None and result.success and result.en:
            # 缓存条目的 parent_chain 必须是 canonical rank，跨 cut 复用才安全。
            result = _canonical_result(named_mol, result)
            _cache_put(cache, smiles, result)
    if not result.success or not result.en:
        return None
    if result.meta and result.meta.get("parent_kind") == "radical":
        # L5 radical worker 已输出带 -yl 的完整自由基名（pentan-1-yl / pyridin-4-yl），
        # 直接消费，不再 free_to_yl 二次转换（否则成 pentan-1-yl-{loc}-yl）。
        return result.en, result.zh, result.en != "phenyl"
    loc = _locant_from_result(sub.mol, result, sub.attach_new)
    if loc is None:
        loc = _locant_via_hetero(sub.mol, sub.attach_new, result)
    return None if loc is None else yl_form(result.en, result.zh, loc)


def name_as_substituent(
    mol, attach_old: int, atoms, *, depth: int = 0, name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    """在 attach_old 处切割原子，free-name 子分子，输出 -yl 双语名称。"""
    atoms = frozenset(atoms)
    sub = build_cut_submol(mol, atoms, attach_old)
    if sub is None:
        return None
    return _yl_from_sub(sub, mol=mol, atoms=atoms, attach_old=attach_old,
                        depth=depth + 1, name_mode=name_mode, cache=cache)
