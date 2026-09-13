"""L3 取代基提取器：提取核心/烷基侧链并汇总 claim 命名结果。"""
from __future__ import annotations

from namepredict.cache.common_names import CommonNameCache

def extract_substituents(info: dict, parent: dict, *, cache: CommonNameCache | None = None) -> list:
    """L3 入口：提取核心取代基、烷基侧链与 claim 侧链并合并。"""
    from namepredict.layer3.claim_extract import extract_claimed_sides

    mol, chain = info["mol"], parent.get("chain") or []
    result =   extract_claimed_sides(info, parent, (), cache=cache)
    print(result,"\n\n\n")
    return result
