"""L4 编号入口：定向编号并组装最终 result 包。"""
from __future__ import annotations
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.layer4.locant_calc import _pack, _with_locants


def number(parent: dict, substituents: list) -> dict:
    """对 parent 定向编号，校验编号骨架事实后组装位次结果。"""
    chain = orient_numbering(parent, substituents)
    kind = parent.get("kind")
    if parent.get("numbering_scaffold_required") and not parent.get("numbering_scaffold"):
        raise ValueError("numbering_scaffold facts required for selected scaffold")
    oriented = {**parent, "chain": chain}
    plan = None
    if plan is not None: oriented["numbering"] = plan

    result = _pack(oriented, _with_locants(chain, substituents, kind, oriented.get("numbering_scaffold")))
    return result
