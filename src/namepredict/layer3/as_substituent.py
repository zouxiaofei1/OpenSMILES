"""name_as_substituent：切割 submol → free-name 管道 → P-29 -yl 形式。
无宿主（苯甲酰胺 gate / n_block 提取）接线——纯 cut→pipeline→yl。
"""
from __future__ import annotations

from rdkit.Chem import CanonicalRankAtoms
from namepredict.constants import C
from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.submol_build import build_anchor_submol, build_cut_submol
from namepredict.tools.free_to_yl import free_to_yl as yl_form

# 简单保留烷氧基作前缀不加括号（ethoxybenzene/phenoxybenzene；与 methoxy 一致，
# 收拢产物 ethoxy/propoxy/butoxy/phenoxy/isopropoxy 由此免括号）。
_SIMPLE_ALKOXY_NO_PAREN = frozenset({
    "methoxy", "ethoxy", "propoxy", "butoxy", "phenoxy", "isopropoxy",
})


def _radical_yl_from_sub(
    mol, atoms: frozenset, attach_old: int, *, depth: int, name_mode: str,
    cache: CommonNameCache | None,
) -> tuple[str, str, bool] | None:
    """碳连接点：锚定 * 走 radical 主基团管线，L4 权威位次 + P-22.2.4 保留名。"""
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
    return hit.en, hit.zh, hit.en not in ("phenyl", *_SIMPLE_ALKOXY_NO_PAREN)


def _yl_from_sub(
     *, mol, atoms, attach_old, depth: int, name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    """连接点类型分派：碳→锚定 radical 优先；非碳/锚定失败→H 封端 free-name + free_to_yl。"""
    from rdkit import Chem
    
    if mol.GetAtomWithIdx(attach_old).GetAtomicNum() >1 :
        hit = _radical_yl_from_sub(mol, atoms, attach_old, depth=depth,
                                    name_mode=name_mode, cache=cache)
        # print(Chem.MolToSmiles(mol),hit)
        if hit is not None:
            return hit
    return None

def name_as_substituent(
    mol, attach_old: int, atoms, *, depth: int = 0, name_mode: str = "general",
    cache: CommonNameCache | None = None,
) -> tuple[str, str, bool] | None:
    """在 attach_old 处切割原子，free-name 子分子，输出 -yl 双语名称。"""
    atoms = frozenset(atoms)
    return _yl_from_sub( mol=mol, atoms=atoms, attach_old=attach_old,
                        depth=depth + 1, name_mode=name_mode, cache=cache)
