"""排序、终态化并返回母体候选（P-44 选择入口）。"""
from __future__ import annotations


def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    """按评分降序排列候选母体。"""
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)

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
    """排序后的母体候选，每个都以不可变 owned_atoms 完成最终确定。"""
    from namepredict.layer2.candidates import _collect_candidates
   
    cands = _finalize_ranked(info, _collect_candidates(info))
    if all_candidates:
        return cands
    return next(iter(cands), None)
