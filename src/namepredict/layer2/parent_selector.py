"""排序、终态化并返回母体候选（P-44 选择入口）。"""
from __future__ import annotations


def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    """按评分降序排列候选母体。"""
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)


def _p45_2_prefix_count(info: dict, parent: dict) -> int:
    """P-45.2.1 键：以前缀引用的取代基团数目（=owned_atoms 边界外的 claim 个数，由 L3 iter_claims 枚举，无 name_mode/cache 依赖）。"""
    mol = info["mol"]
    from namepredict.layer3.claimable_block import iter_claims
    return len(iter_claims(mol, parent.get("owned_atoms") or frozenset()))


def _reorder_p45_2(info: dict, cands: list[dict]) -> list[dict]:
    """P-45.2 流水线：按前缀取代基团数目最多（P-45.2.1）稳定重排打平候选；P-45.2.2/2.3 的位次需 L4 编号后才可得，只先接 P-45.2.1。"""
    if len(cands) <= 1:
        return cands
    return sorted(cands, key=lambda c: _p45_2_prefix_count(info, c), reverse=True)


def _finalize_ranked(info: dict, cands: list[dict]) -> list[dict]:
    """补齐契约/词干/编号并固化 owned_atoms。"""
    from namepredict.layer2.kind_registry import pack_parent_stem
    from namepredict.layer2.parent_candidate import with_principal_group_contract
    from namepredict.layer2.parent_ownership import finalize_parent_ownership
    mol = info.get("mol")
    return [
        finalize_parent_ownership(
            pack_parent_stem(with_principal_group_contract(c), mol), mol,
        )
        for c in _rank_candidates(info, cands)
    ]

def select_parent(info: dict, *, all_candidates: bool = False) -> dict | list[dict] | None:
    """排序并终态化母体候选（P-44 评分降序 → P-45.2 重排打平 → owned_atoms 固化）：取首位为选定母体，all_candidates 时全返回。"""
    from namepredict.layer2.candidates import _collect_candidates

    cands = _finalize_ranked(info, _collect_candidates(info))
    cands = _reorder_p45_2(info, cands)
    if all_candidates:
        return cands
    return next(iter(cands), None)
