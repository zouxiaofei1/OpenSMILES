"""FG 位次省略规则（L4；P-14.3.4 / 环单 FG）。"""
from __future__ import annotations


def _is_cyclo(parent: dict | None) -> bool:
    """环状单环 carbocycle；有 fused_tree 的稠环不算单环环烷烃。"""
    return (parent or {}).get("scaffold_id") == "carbocycle" and not (parent or {}).get("fused_tree")


def omit_fg_locant(
    pos: int | None, n_carbons: int, parent: dict | None = None, n_subs: int = 0, *,
    single: bool = True,
) -> bool:
    """判定主官能团位次是否省略（P-14.3.4 / 环单 FG）：环状无取代省。"""
    if _is_cyclo(parent) and pos is not None and single:  # 环单 FG：无取代省位次，有取代或多官能团保留。
        return n_subs == 0
    return pos == 1 and n_carbons <= 2


def omit_unsat(
    n_carbons: int, kind: str | None = None, parent: dict | None = None, *,
    triple: bool = False,
) -> bool:
    """判定烯/炔位次是否省略（环单烯或短链）；triple 选择炔规则。"""
    if kind == "alkane" and (parent or {}).get("scaffold_id") == "carbocycle":  # 纯烃环单烯位次隐含省略；环多烯保留位次。
        if not (parent or {}).get("double_bonds"):
            return True
    # 烯 ≤C2、炔 ≤C3 时位次 '1' 省略（P-14.3.4.2(d)）。
    return n_carbons <= (3 if triple else 2)
