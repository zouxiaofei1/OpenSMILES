"""L3 取代基提取器：提取核心/烷基侧链并汇总 claim 命名结果。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.cache.common_names import CommonNameCache
from namepredict.constants import HALO_EN, HALO_ZH
from namepredict.tools.chain import carbon_neighbors
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer1.functional_group_inventory import FunctionalGroupClass, inventory_from_info

def extract_substituents(info: dict, parent: dict, *, name_mode: str = "general", cache: CommonNameCache | None = None, depth: int = 0) -> list:
    """L3 入口：提取核心取代基、烷基侧链与 claim 侧链并合并；depth 透传给 claim 侧链命名。"""
    from namepredict.layer3.claim_extract import extract_claimed_sides

    mol, chain = info["mol"], parent.get("chain") or []
    result =   extract_claimed_sides(info, parent, (), name_mode=name_mode, cache=cache, depth=depth)
    return result
