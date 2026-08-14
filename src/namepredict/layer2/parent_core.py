"""L2 母体组装 / 链行走 / 门控辅助（单一权威）；生产者从这里导入辅助函数，parent_selector 保留 FG try 与 select_parent 并可能为兼容性做精简再导出。"""
from __future__ import annotations

from namepredict.layer2.chain_walk import _longest_chain, _longest_from


def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, **kw}
