"""L3 取代基提取器：提取核心/烷基侧链并汇总 claim 命名结果。"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.cache.common_names import CommonNameCache
from namepredict.constants import HALO_EN, HALO_ZH
from namepredict.tools.chain import carbon_neighbors
from namepredict.layer1.fg_registry import FG_SPECS


_PARENT_OXO_KINDS = frozenset(sp.fg for sp in FG_SPECS if sp.oxo_parent)


def _make_oxo(attach: int) -> dict:
    """构造氧代取代基字典。"""
    return {
        "kind": "oxo", "attach_idx": attach, "atoms": [attach],
        "en": "oxo", "zh": "氧代",
    }


def _extract_oxos(info: dict, parent: dict) -> list[dict]:
    """提取链上的氧代取代基。"""
    if parent.get("kind") in _PARENT_OXO_KINDS:
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_oxo(k["c_idx"]) for k in info.get("ketones") or [] if k["c_idx"] in chain
    ]



def extract_substituents(info: dict, parent: dict, *, name_mode: str = "general", cache: CommonNameCache | None = None, depth: int = 0) -> list:
    """L3 入口：提取核心取代基、烷基侧链与 claim 侧链并合并；depth 透传给 claim 侧链命名。"""
    from namepredict.layer3.claim_extract import extract_claimed_sides

    mol, chain = info["mol"], parent.get("chain") or []
    base = (
        _extract_oxos(info, parent)
    )
    result = base + extract_claimed_sides(info, parent, base, name_mode=name_mode, cache=cache, depth=depth)
    # print(result)
    return result
